# Batch 6A — Unified Evidence Merger

## Purpose

将既有扫描结果按 **保留命中的并集** 整合为每条蛋白的 Evidence Profile。
不重新搜索、不修改阈值，不开发功能分类、候选评分或 RN28 专属规则。
RN28-positive / CTD-negative、仅 PF13573 命中、仅 PorP/SprF 命中的蛋白都收录；
低于 GA、未达到 HMM domain reporting 条件但已由后端返回的 domain 也保留。

## Input / dependencies

- 一个完整且 checksum/质量检查通过的 prepared dataset。
- 可选：`bctd-hmm-search`、`sprf-search` 的 paired 输出、`sprb-repeat-search` 的 PF13573 输出、
  含 RN28 query 的 `mmseqs-search` 输出。省略的输入标记 `not_run`。
- RN28 成功扫描需要本地 Reference Registry；只接收其中注册的 RN28 **相同序列**，允许 query ID 别名。
  混合 query 搜索中，只导入 RN28 对应的 alignment；其他 query 不被误称为 RN28，完整原始文件仍由元数据追溯。
- 只读现有文件；使用 core/io、Reference Registry、已有 DomainEvidence、SequenceAlignment 和
  NeighborhoodIndex。不执行 MMseqs2/PyHMMER，不要求模型二进制、MMseqs2 可执行文件或网络。
- 不读取脱离来源/manifest 的单独 TSV；不依赖调用时原始 query 文件或被删除的临时 MMseqs2 数据库。

## Usage / parameters

在项目根目录运行，或者显式指定 `--registry /absolute/path/references/registry.json`：

```bash
uv run sprb-scanner merge-evidence \
  --prepared /path/genome/prepared \
  --bctd-scan /path/genome/bctd \
  --sprf-scan /path/genome/sprf \
  --repeat-scan /path/genome/repeats \
  --rn28-scan /path/genome/rn28 \
  --registry references/registry.json \
  --window 5 --output /path/new-evidence
```

`--window` 默认上下游各 5 个 CDS（非负整数）；`--no-neighborhood` 只整合扫描证据。
各扫描参数、模型版本、checksum 和阈值使用已有输出的原值。
批量工作流可以逐个调用 CLI 或 Python API，不需要重复运行搜索：

```python
from pathlib import Path
from sprb_scanner.evidence import run_merge

summary = run_merge(
    Path("/path/prepared"),
    Path("/path/new-evidence"),
    bctd_scan=Path("/path/bctd"),
    sprf_scan=Path("/path/sprf"),
    repeat_scan=Path("/path/repeats"),
    rn28_scan=Path("/path/rn28"),
    registry=Path("/path/references/registry.json"),
    window=5,
)
```

`merge(dataset, runs, hmm_hits, alignments, ...)` 提供内存整合；`run_merge` 负责磁盘输入验证与导出。
重复运行请选择新的输出目录，不覆盖既有结果。

## Evidence and status contract

每个 profile 包含 protein/source ID、全长长度、序列 SHA-256、原始 GeneRecord（含分段和 qualifiers），
以及 `sources` 中六个独立来源：TIGR04131、PF13585、PF13573、TIGR03519、PF11751、RN28。

| 层级 | 状态 | 含义 |
|---|---|---|
| genome/source | `not_run` | 没有提供此来源，命中数量为 null |
| genome/source | `failed` | 提供有效失败输出，保留错误和元数据；命中数量为 null |
| genome/source | `succeeded` | 成功输出经过验证；命中数量可以为 0 |
| protein/source | `not_run` / `failed` | 继承来源状态，不表示生物学阴性 |
| protein/source | `succeeded_no_hits` | 来源成功，但该蛋白无保留命中 |
| protein/source | `succeeded_with_hits` | 来源成功，保留该蛋白的全部记录 |

GA 在各条 HMM domain 中独立记录 `passed_ga`、`below_ga`、`unavailable`。
Protein state 不等于 GA state；一个蛋白可以有保留命中但全部低于 GA。
HMM 的 sequence E-value / bit score、independent/conditional domain E-value、domain bit score 原样保留，
不与 MMseqs2 alignment E-value / bit score 混合或比较。
每个 domain 保留 HMM/target coverage、alignment/envelope/HMM 坐标、相对位置及重叠 domain ID。
RN28 保留每段 identity、coverage、alignment length、分数和坐标，不只取最佳匹配。
PF13573 数量是 **后端保留的 domain 记录数**，不等同于不重叠物理重复单元数量。
不合并重叠区间，不将覆盖率相加，不按全长 coverage、长度或 GA 收窄并集。

## Neighborhood semantics

为证据并集中的每条有 CDS 的蛋白查询现有 NeighborhoodIndex，使用成功 PorP/SprF 来源的保留命中。
每条关系保留 neighbor gene/protein ID、原始分段、strand、坐标/转录相对方向、CDS offset、
intervening CDS、直接邻接、distance/overlap/signed gap 和 circular origin 信息。
关系通过 domain hit ID 链接到邻近蛋白自己的完整证据 profile；不把邻近蛋白的分数赋给 cargo。

| state | outcome | 含义 |
|---|---|---|
| `succeeded` | `local_hit` / `no_local_hit` | 两个 SprF 模型均成功，邻域查询完整 |
| `partial` | `local_hit` / `unknown` | 单条蛋白：一个 SprF 模型成功，另一个缺失或失败；观测到的数量是下界 |
| `not_run` | `unknown` | 禁用邻域，或没有提供任何 SprF 来源 |
| `failed` | `unknown` | SprF 没有成功来源但存在失败，或 CDS 索引/几何错误 |
| `unavailable` | `unknown` | 仅蛋白序列，无基因组 CDS 坐标 |

全局邻域状态即使候选数为零也输出。异常 CDS 的错误信息保留，不删除 CDS 或回退修改坐标。
索引构建失败仍影响所有候选的邻域状态。索引构建成功后，逐候选 `query()` 的 `InputError`
只将该蛋白的邻域标为 `failed / unknown`，保留 protein/CDS ID、原始错误及全部扫描证据；
其他候选继续查询并正常导出。失败候选不输出部分邻接关系或假零命中数量。
存在局部查询失败时，全局邻域状态汇总为 `partial`（仍有候选正常查询）或 `failed`（全部查询失败），
保留失败数量与逐蛋白错误，`status.json` 的 `evidence_complete` 为 false。
此汇总不会将其他蛋白已完成的邻域状态改为失败，也不会将模型缺失标记改成阴性。
未知程序异常仍向上传播，不作为可忽略的局部输入异常。
`no_local_hit` 不等于 orphan；任何 PorP/SprF HMM 命中都不自动代表功能性 shuttle。
RN28-only 或 SprF-only 蛋白的邻域同样只是客观关系，不会被自动称为 cargo。
SprD 的局部 domain coverage 保持独立，不自动判定为典型全长 SprF 或假阳性。

## Output

- `evidence_profiles.json`：schema 1.0.0，全局六来源状态、邻域状态和所有 profile/domain/alignment/关系。
- `evidence_profiles.tsv`：每条蛋白一行，六来源状态、保留记录数、GA domain 数、邻域状态和邻接数量。
  未运行/失败时的数量留空，成功零命中时写 0；RN28 的 GA 数量留空。
- `evidence_records.tsv`：每个 domain/alignment 一行，分数、坐标、coverage、GA、模型身份均有独立列；
  另保留完整 `evidence` JSON 对象和原始 scan directory，不丢失专有字段。
- `neighborhood_associations.tsv`：每条候选–邻近 PorP/SprF-like 蛋白关系一行，允许一对多。
- `summary.json`：候选并集数量、每来源状态/记录数、蛋白状态计数和全局/逐蛋白邻域状态。
- `status.json`：整合执行成功与否、`evidence_complete`。执行成功不代表所有来源都运行成功。
- `run_metadata.json`：版本、参数、模型/序列来源、原始 run metadata、全部读取文件 SHA-256 及输出 checksum。
  原始 raw 文件通过输入 checksum + `scan_directory` / `raw_artifact` 引用，不复制大型原始数据。

全部内部蛋白/基因组区间为 0-based half-open。JSON 含既有 computed 字段；需要 Pydantic 读回时，
使用内存结果的 `model_dump(exclude_computed_fields=True)`，与现有数据约定一致。

## Error handling

显式提供但不存在的路径、损坏 JSON/checksum、越界/重复 hit ID、不匹配蛋白 ID/序列指纹、
paired/native 冲突或不含注册 RN28 的 query 都报错（CLI exit 2），不伪装成 `not_run` 或零命中。
格式完整的 failed run 可整合为 failed source；若失败目录同时含歧义的未完成命中文件则拒绝，
不将这些文件作为完整证据。早期失败没有目标指纹时，只记录失败尝试，不附加任何命中。
所有已读取文件在导出前再次验证未变化；临时目录写完后原子改名；校验错误不产生成功输出。

## Tests / real-data validation

`tests/unit/test_evidence.py` 覆盖来源并集、三条 RN28-only 候选、弱/未报告 domain、多个 alignment、
重叠 domain 不累加、六来源状态、部分/完整/失败邻域、逐候选查询异常隔离与完整导出、
环状 CDS、无坐标输入、空输出、CLI、
SHA-256/ID/序列/模型/坐标/coverage 校验和不覆盖输出。

2026-10-20，复用 21 株现有五模型 HMM 和 RN28 输出，未执行新搜索：

- 21/21 成功整合，候选并集 **592 条蛋白**，完整保留 **2,807 条 HMM domain**。
- RN28 **17 条蛋白 / 18 段 alignment** 原样保留。
- `WP_013551150.1`、`WP_072355989.1`、`WP_229206235.1` 三条 RN28-positive / CTD-negative 蛋白仍在并集中。
- FJ 12 条参考 cargo 恢复 **9 个邻接型 + 3 个 no_local_hit**；SprF/SprD 的全长/局部 coverage 差异保留。

真实基因组及完整结果保存在仓库外，不纳入 Git；这些数值是该批特定输入的验收记录，不是新的分类阈值。
