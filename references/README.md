# Unified Reference Registry

## Purpose

`references/registry.json` 是序列、序列集合、MSA、Profile HMM 和 annotation catalog 的权威索引。
`src/sprb_scanner/references/` 提供离线解析、类型/校验值验证、版本查找和来源依赖追溯。
本批没有 MSA/HMM 构建、自动下载、序列重切分、模型 GA 调整、分类或扫描结果重写。

## Layout / migrated input

```text
references/
  registry.json          # 20 个固定版本资源
  registry.schema.json   # 由 RegistryDocument 导出的 JSON Schema
  migration.json         # 16 个历史文件的路径/字节数/SHA-256 对照
  sprb/                  # 原 reference_data/ 的逐字节副本
    RN28.fasta
    reference_only.fasta
    reference_catalog.json / reference_catalog.tsv
    reference_manifest.json
    archive/...
  hmm/                   # 原 profiles/ 的逐字节副本
    *.profile.json
```

`RN28` 是主 reference ID，FASTA query ID 保持 `SMF_RN28__RN28`。
序列仍为 **207 aa**，位于 FJ SprB `WP_012023065.1` 的 **[681,888)**。
原始 panel 的 74 条序列、模块边界、46 个 legacy clusters 及全部注释行保持不变。
其他 73 条序列是 reference-only；没有隐式启用全 panel 搜索。

注册模型：`TIGR04131@TIGR04131.1`、`PF13585@PF13585.13`、
`TIGR03519@TIGR03519.1`、`PF11751@PF11751.15`、`PF13573@PF13573.13`。
原始 source URL、accession 和 **HMM 文件** SHA-256 不变；配置 JSON 作为单独 catalog 注册，
其文件 checksum 与模型 bytes checksum 明确区分。模型文件不提交 Git，初始状态是 missing_external。

## Schema / identity and version control

每个条目必须包含：

| 字段 | 合同 |
|---|---|
| `reference_id`、`version` | 唯一身份 `(reference_id, version)`；非空、无空白 |
| `resource_type` | sequence / sequence_set / msa / profile_hmm / annotation_catalog |
| `artifact` | `{path, format, storage}`；path 相对 Registry 文件目录，POSIX 风格，不含 `..`、绝对路径或 URL；storage 为 bundled/external |
| `sha256` | 原始 artifact **文件字节**的 64 位小写 SHA-256 |
| `source` | provider、locator，以及可选 release、notes；locator 是来源记录，不触发访问 |
| `biological_target` | description，以及可选 organism、parent accession/length、protein_interval、sequence IDs、sequence SHA-256、model ID/accession |
| `derived_from` | 可选父资源列表，每个父项固定 reference_id、version 和 sha256 |

`biological_target.protein_interval` 使用 0-based half-open。
`sequence_sha256` 用于单序列的 residue bytes，与 FASTA 文件 checksum 分开。
模型来源与自定义 RN/RC/SMF 标签都按原始含义保留，不将旧 cluster 解释为已验证的生物学分类。

同一 ID 可以有多个版本；查找时若不指定版本且有歧义则报错，不自动选择 latest。
新增数据/构建产物应新增版本和 artifact 路径，保留既有版本的文件、checksum 和 provenance。
父资源必须已注册、checksum 匹配；拒绝重复身份/父项、缺失父项、checksum 冲突与依赖环。
`lineage()` 返回传递父资源，供未来 `sequence_set → msa → profile_hmm` 工作流使用。
它不要求所有历史父文件仍在本地；这只说明 lineage 声明已固定，不等于重建验证已完成。

## Dependencies / API / CLI parameters

沿用锁定的 Pydantic、Biopython、PyHMMER 依赖；不新增下载器或扫描算法。
所有操作离线执行，默认 Registry 路径是当前目录的 `references/registry.json`；其他位置显式传 `--registry`。

```bash
uv run sprb-scanner reference-list
uv run sprb-scanner reference-validate
uv run sprb-scanner reference-resolve --reference-id RN28
```

`reference-list` 只解析身份/schema/依赖，不要求模型文件存在。
`reference-validate` 检查所有本地资源：bundled 缺失为失败；external 缺失明确标记 missing_external。
默认允许 external 缺失，整体状态 partial、exit 0；这不表示缺失模型已验证。
`--require-external` 将任何 external 缺失视为失败、exit 2。
损坏的 external 文件即使默认允许缺失也会失败，不会当作缺失而略过。

将已下载模型交给 Registry，可覆盖它的本地 artifact 路径；不会修改索引或模型：

```bash
uv run sprb-scanner reference-resolve --reference-id TIGR04131 \
  --artifact /path/to/TIGR04131.1.HMM \
  --export-profile /path/to/new-TIGR04131.profile.json
```

`--version` 指定固定版本。`--artifact` 只允许用于 external 资源。
`--export-profile` 只允许 verified profile_hmm，创建新的标准 ProfileConfig JSON，不覆盖已有文件。
适配到旧引擎时，ProfileConfig 的模型版本使用实际 HMM accession（若有版本后缀），
Registry 的资源发布版本仍保留在 Reference 身份中，不把二者混同。
之后仍用原命令：`hmm-search --prepared ... --profile /path/to/new-TIGR04131.profile.json --output ...`。
也可把模型放到 canonical `references/hmm/` 对应文件名处，继续使用那里的原始 `.profile.json`。
原始 `--profile`、`--tigr-profile`、`--pfam-profile` 调用格式及自定义 ProfileConfig 均不变。

序列搜索继续显式传 FASTA：

```bash
uv run sprb-scanner mmseqs-search --query references/sprb/RN28.fasta \
  --prepared /path/to/prepared --output /path/to/new-result
```

Python API：

```python
from sprb_scanner.references import ReferenceRegistry

registry = ReferenceRegistry.load("references/registry.json")
entry = registry.get("RN28", "sprb-legacy-20261019")
resolved = registry.resolve("RN28")  # reference、已验证的 artifact_path、解析 summary
proteins = registry.sequences("RN28")  # 保留原 FASTA IDs 的 ProteinRecord tuple
parents = registry.lineage("RN28")
profile = registry.to_profile("TIGR04131", artifact="/path/to/TIGR04131.1.HMM")
```

## Type validation / output

- sequence：严格单条 ungapped protein FASTA；不改写 lowercase/stop/gap 为合法序列。
- sequence_set：一条或多条 protein FASTA，保留顺序，拒绝重复 ID/空记录/非法残基。
- msa：支持 aligned FASTA、Stockholm、A3M；验证非空、ID 唯一、氨基酸/gap、match-column 长度。
  A3M 的 lowercase insertion 和 `.` insertion gaps 不计入 match columns。只解析，不做 MSA 构建。
- profile_hmm：必须为单个合法 amino-acid HMM，检查声明的 model ID/accession，原样读取 GA。
- annotation_catalog：支持严格 JSON object/array 和有唯一表头、矩形行的 TSV。保留历史 row 原值；
  这里只校验资源格式，不对任意 catalog 发明新的语义或重算模块边界。

先读取一次文件字节、验证 SHA-256，再从同一字节快照解析。解析不能自动修复 checksum。
CLI 输出 JSON 包含身份、版本、来源、实际本地路径、摘要及明确的状态/错误，不生成搜索命中。
Registry 不自动插入或替换历史扫描结果的字段；既有搜索 metadata 及 query ID 行为保持不变。

## Compatibility / migration policy

旧 `reference_data/` 和 `profiles/` 保留 **固定历史版本的兼容副本**，不是第二套可独立修改的资源配置。
新代码/文档以本目录的 Registry 为准。迁移表与回归测试保证两个路径下的原始文件逐字节一致；
不使用目录 symlink，避免破坏旧 manifest 所在目录的相对路径语义、复制到模型目录的用法及 CI。

`sprb/README.md`、`sprb/reference_manifest.json` 和 `hmm/README.md` 是按原 SHA-256 保存的历史文档/快照，
其中旧相对路径按旧布局解释；新 Registry 只按自己的 artifact path 定位资源，不把历史 manifest
中的 fixture 路径重新解释为 canonical 路径。本文件是新布局的使用文档。
旧例子 `--query reference_data/RN28.fasta` 与 `--profile profiles/NAME.profile.json` 仍可用。
新版本应新增 Registry entry，不追溯性重写这些历史兼容副本。

源码包携带 canonical 数据与兼容快照；wheel 包含 Python Registry 模块，但没有隐藏的资源路径发现 API。
参考/模型参数仍由用户显式提供；不提交大型 HMM、基因组、MMseqs2 数据库或批量结果。

## Error handling / test cases

错误统一为 `ReferenceRegistryError`，缺失文件为 `ReferenceUnavailableError`；CLI exit 2。
错误包含 reference ID/version、实际路径及 checksum/类型/图依赖原因。不会联网寻找缺失模型。
bundled 资源的符号链接不能逃出 Registry 目录；external 路径只通过显式本地 override 或 artifact 定位。

`tests/unit/test_reference_registry.py` 覆盖所有格式、checksum、坏文件/路径/ID、缺失 external、
DAG/版本歧义和 ProfileConfig 导出，验证原 HMM engine 可读且 GA/bytes 不变。
`tests/unit/test_reference_migration.py` 覆盖16个历史文件、74个模块/全部边界和legacy assignments、
RN28 207 aa/[681,888)、五模型来源/版本/checksum、JSON Schema 同步和完整资源清单。
原 `tests/unit/test_reference_data.py` 继续检查旧路径，保证历史命令的文件契约未破坏。

本批本地验收：315 项 pytest（含完整 FJ 回归）通过；五个本地 HMM 只读校验的来源、bytes 和 GA 均不变。
