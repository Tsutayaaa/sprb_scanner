# Batch 5A — Generic MMseqs2 Search Engine

## Purpose

通用蛋白序列/片段搜索，复用 `ProteinRecord`、`DomainHit` 与 JSON/TSV 导出体系。
`scanners.mmseqs.SequenceAlignment` 继承 `DomainHit`，额外保存 identity、alignment length 和来源 ID。
不改变共享 `DomainHit` schema、不修改 HMM Engine/NeighborhoodIndex、不推断家族功能或融合证据。

## Dependencies

Python 依赖沿用 `uv.lock`。另需 CPU 版 [MMseqs2](https://github.com/soedinglab/MMseqs2)，
可用 `brew install mmseqs2`（macOS）或 `sudo apt-get install mmseqs2`（Ubuntu）。
本地验证版本 `18-8cc5c`；CI 显式安装并运行真实后端测试，打印实际版本。
通过 `--executable /path/to/mmseqs` 选择安装。保留 launcher 符号链接路径，兼容 Debian/Ubuntu SIMD dispatcher。搜索/转换显式使用 `--threads`；建库先检查 `createdb -h`，支持线程参数时才传入。
MMseqs2 15 的 createdb 为单线程且不接受该参数，有效建库线程数和后端能力记入元数据。
不自动下载 MMseqs2、基因组或外部参考蛋白。仓库内的小型用户参考集见
[reference_data](../reference_data/README.md)，RN28 主 query 为 `reference_data/RN28.fasta`。

## Input / usage

查询为 `.faa/.fasta/.fa/.fas`（也支持 `.gz`），1 个或多个 **唯一 FASTA header 首字段 ID**。
可以混合全长蛋白和短片段；description 不作为 ID。小写转大写、一个末端 `*` 去除的操作记入元数据。
空记录、重复 query ID、内部 stop、gap 或非法残基会使整次运行失败，不丢弃后继续。
`B/X/Z/J/U/O` 沿用 ProteinRecord 的允许范围；实际打分/字母解释遵循 MMseqs2，自行保留原始序列和 alignment。
短且低复杂度的序列可能被 MMseqs2 预过滤漏掉；本模块不声称穷举所有可能的局部比对。

直接搜索一个经过校验的 prepared proteome：

```bash
uv run sprb-scanner mmseqs-search --query references.faa \
  --prepared /path/FJ-prepared --output /path/FJ-sequence-search \
  --sensitivity 7.5 --evalue 0.001 --threads 2
```

反复使用相同蛋白组时，先构建数据库：

```bash
uv run sprb-scanner mmseqs-createdb --prepared /path/FJ-prepared \
  --output /path/FJ-mmseqs-db --threads 2
uv run sprb-scanner mmseqs-search --query references.faa \
  --database /path/FJ-mmseqs-db --output /path/new-search --threads 2
```

数据库是本项目 `mmseqs-createdb` 生成的 **managed directory**，包含原生数据库、完整 ProteinRecord
映射、标准化序列 FASTA、来源、版本和 SHA-256。可以移动后复用，无需原 prepared 目录。
不直接接受缺少映射/来源的外部裸数据库前缀。更换 MMseqs2 版本后需重新构建。
同时传入 `--prepared` 和 `--database` 时，核对 ordered `(protein_id, sequence)` 指纹；不一致即失败。
内部使用 `q0/t0` 等安全 token，输出还原 query 原 ID、target 内部 ID、source ID、gene ID。
因此带 `|`、冒号等 ID 不受后端 header 解析影响，重复 target source ID 也不会合并内部 ID。

## Parameters / retention

| 参数 | 默认值 | 作用 |
|---|---|---|
| `--sensitivity` | 5.7 | MMseqs2 `-s`，范围 1–7.5 |
| `--evalue` | 0.001 | MMseqs2 alignment 报告阈值，正有限值 |
| `--threads` | 1 | 正整数 |
| `--max-seqs` | target 蛋白数 | 预过滤结果上限，避免隐式默认 300 截断 |
| `--alternative-alignments` | 最长 target 的长度 | 每对序列最佳比对之外，继续搜索 target 被掩蔽后的局部比对 |
| `--split-memory-limit` | 2G | MMseqs2 每 split 内存限制，不是整个进程 RSS 的硬上限 |
| `--timeout-seconds` | 3600 | 每个后端命令的超时；POSIX 上终止整个流程进程组 |
| `--executable` | mmseqs | PATH 名称或可执行路径 |

显式使用 `-a 1 --alignment-mode 3 --seq-id-mode 0`，获取 backtrace 和实际相同残基计数。
`--max-accept/--max-rejected` 为 INT_MAX，coverage、identity、最小 alignment length 阈值均为 0。
最大序列长度至少 65,535，按实际 query/target 最长长度增大，不故意截短长蛋白。

MMseqs2 默认 `--alt-ali 0` 只返回每对序列最佳匹配。本模块默认开启 alternative alignments，
采用 target 长度决定的有限上限；后端会在不再通过报告条件时停止。
元数据记录有效上限及 `pairs_reaching_alternative_limit`（可能触及上限的 pair 数）。
如用户显式设为 0，则只搜索最佳 alignment。所有 `convertalis` 行逐条保留，不折叠、不去重、
不将多个 alignment 的 coverage 相加。MMseqs2 的预过滤、target masking 及记录的容量参数仍限定
搜索范围；这不是所有重叠 HSP 的穷举器。默认 composition correction 和低复杂度 masking 保持后端行为，
日志保留完整的后端参数展开。进一步行为见 [官方指南](https://github.com/soedinglab/MMseqs2/wiki)。

## Alignment schema / coordinates

| 字段 | 定义 |
|---|---|
| `query_id`, `query_source_id` | query FASTA 原始首字段 ID |
| `protein_id`, `target_source_id`, `target_gene_id` | target 内部 ID、原注释 ID 和 CDS 链接 |
| `query_length`, `protein_length` | query / 完整 target 长度 |
| `query_start`, `query_end` | query 上 0-based half-open 区间 |
| `start`, `end` | target 上 0-based half-open 区间 |
| `alignment_length` | alignment 列数，包括 gap 列 |
| `identical_residues`, `sequence_identity` | 后端 nident，及 nident / alignment_length，identity 为 0–1 |
| `query_coverage` | (query_end − query_start) / query_length |
| `target_coverage` | (end − start) / **整条 target** 长度 |
| `evalue`, `evalue_kind`, `bit_score` | MMseqs2 alignment 统计；kind 为 alignment，不能与 HMM domain E-value 混用 |
| `model` | query ID、规范化 query 序列 SHA-256、来源 FASTA 路径；version 为 sequence-sha256:<hash> |
| `backend_version`, `raw_artifact` | MMseqs2 实际版本与原始 alignment TSV |
| `hmm` | null |

`convertalis` 原始坐标为 1-based inclusive，转换为 `[start-1, end)`。
Span coverage 包括插入/缺失影响的坐标跨度，不等于 identity，也不是只数匹配列。
短片段 100% query coverage / 很低 target coverage 是正常结果，**不会因 target coverage 低而过滤**。

## Output / reproducibility

新建输出目录，不覆盖已有结果；允许相同 query/数据库重复运行到不同目录。

- `alignments.tsv`、`alignments.json`：逐条标准化完整记录，唯一 hit ID；空结果保留 TSV 表头及 `[]`。
- `raw_alignments.tsv`：无表头的原始 convertalis 输出；列顺序在元数据 `raw_columns` 中。
  包括后端 fident/qcov/tcov、qaln/taln；规范化记录用整数区间/计数重算 coverage/identity，避免输出舍入误差。
- `query_mapping.json`、`target_mapping.json`：安全 token → 原始 ProteinRecord（含序列）。
- `status.json`、`run_metadata.json`、`backend.log`：成功/失败状态、软件/后端版本、完整命令、
  用户与有效参数、来源 prepared 数据/元数据/FASTA SHA-256、蛋白指纹、数据库 manifest SHA-256，
  以及输出文件 SHA-256（不对 manifest 自身作循环校验）。
- 直接 prepared 模式额外保存 `direct_database_metadata.json` / `direct_database_backend.log`，
  保留临时目标数据库创建的来源及失败诊断；原生临时 query/result DB 用毕清理。

数据库复用在搜索前后核验版本、文件 checksum、manifest 不变、token/长度/蛋白指纹。
每条 alignment 核对 query/target token、长度和坐标边界；不接受未知 ID 或越界记录。
JSON 与现有导出风格一致含 computed `ctd_distance`；Pydantic 读回时先去除 computed 字段，
或在内存记录上使用 `exclude_computed_fields=True`。

## Error handling

输入缺失、空 target、query 歧义、prepared 质量 error、损坏/不兼容数据库、未知 ID、
错误/越界 alignment、后端失败和超时均返回 CLI exit 2。
配置验证通过、输出目录建立后，正常可报告失败会保留 failed 状态、错误说明和已有日志/元数据。
配置错误、已存在输出目录或无法写入输出目录等前置错误不覆盖/强行创建结果。
成功零命中明确标记 succeeded + alignment_count=0。

## Test cases / real-data acceptance

`tests/unit/test_mmseqs.py`：gapped 区间与 identity、独立 coverage、重复行保留、错误输入、
空/坏行、配置、失败目录、超时、数据库 metadata。
`tests/integration/test_mmseqs_cli.py`：真实 MMseqs2 下 120 aa 全长自匹配、90 aa 片段对
315 aa 蛋白两个重复区段的匹配（query coverage 1，target coverage <0.3），
数据库移动/复用/重复运行、空结果、版本/校验/身份不匹配与数据库构建失败日志。
本地没有 MMseqs2 时只跳过真实后端测试；CI 设置 `SPRB_REQUIRE_MMSEQS=1`，缺依赖会失败。

FJ smoke test 独立于默认 pytest，不捆绑完整蛋白组：

```bash
uv run python scripts/smoke_mmseqs_fj.py --prepared /path/FJ-prepared \
  --output /path/new-fj-mmseqs-smoke --threads 2
```

2026-10-19：FJ UW101 `GCF_000016645.1`，5,139 条蛋白，SprB `WP_012023065.1`（6,497 aa）。
query 为全长与 C 端 218/448/1,182 aa，对应
[Kulkarni et al., 2019](https://journals.asm.org/doi/10.1128/JB.00218-19) 报道的片段长度；
这里提取原生 RefSeq 末端序列，不复刻 GFP fusion/linker，也不把它们称为 RN28。
默认 E-value 0.001、sensitivity 7.5、2 threads，共保留 122 条 alignment。

| Query | 预期 target 区间 | Identity / query coverage | Target coverage |
|---|---|---|---|
| SprB full | [0,6497) | 100% / 100% | 100% |
| C218 | [6279,6497) | 100% / 100% | 3.3554% |
| C448 | [6049,6497) | 100% / 100% | 6.8955% |
| C1182 | [5315,6497) | 100% / 100% | 18.1930% |

四个预期精确匹配全部恢复。完整结果、query/来源边界和 SHA-256 保存在仓库外的
`~/.local/share/sprb-scanner/analyses/batch5a-mmseqs-20261019/`。
E-value 数值为 0 时是后端输出的数值下溢，不代表额外的生物学结论。

本地完整验收：258 项 pytest 全部通过（含完整 FJ 注释与12-cargo回归）；Ruff、CLI help 和打包通过。
