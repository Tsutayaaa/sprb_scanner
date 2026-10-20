# RN28-first reference data

## Purpose

项目内的小型 SprB reference 数据，优先使用 **RN28.fasta** 作为序列搜索 query。
其他 73 条序列先作为参考保存，不自动参与搜索，不导入旧项目的评分或过滤规则。
RN/RC/SMF 名称与 cluster 是旧项目的自定义标签，不代表官方 Pfam 家族或已验证的功能分类。

## Input / provenance

来源为用户提供的旧项目 `sprb_genome_scan/reference_data`，四个源文件在 `archive/` 中逐字节保留。
导入时对照 FJ UW101 `GCF_000016645.1` 的 SprB `WP_012023065.1`（6,497 aa）：
74 条参考序列均与旧 domain table 指定的蛋白区间精确一致。
来源与输出 SHA-256 记录在 `reference_manifest.json`，逐条序列/坐标/cluster 记录在 catalog 中。
manifest 的实际文件路径相对此目录；parent 校验使用现有 `../tests/fixtures/fj_window.faa`。
原始 `cluster_summary.tsv` 中的旧绝对路径只作历史来源记录，使用本数据不需要访问这些路径。

## Output / layout

| 文件 | 用途 |
|---|---|
| `RN28.fasta` | 主 query，只有 `SMF_RN28__RN28`，207 aa |
| `reference_only.fasta` | 其他 73 条序列，仅供参考 |
| `archive/sprb_smf_panel.fasta` | 原始完整 74 条 panel，保留顺序/ID/字节 |
| `archive/domain_table.tsv` | 原始模块边界表 |
| `archive/cluster_assignments.tsv`、`archive/cluster_summary.tsv` | 原始 cluster 记录 |
| `reference_catalog.tsv`、`reference_catalog.json` | 角色、蛋白 ID、序列 hash、坐标、所有精确匹配位置及原始记录 |
| `reference_manifest.json` | 主 query 选择、来源、parent 和文件 SHA-256 |

RN28 对应 SprB **[681,888)**（0-based half-open），即旧表 682–888（1-based inclusive）。
它是 SprB 内部片段，与文献 C-terminal 片段 query 分开。旧表 expanded 678–892 为 215 aa，
包含两侧各 4 aa 延伸；主 query 保持原始 207 aa，不将 expanded 坐标当作其边界。
重复模块可能在 SprB 中精确匹配多处；catalog 保留所有位置，parent_start/end 使用旧表指定的区间。
`source_module_row` 中原始 start/end/expanded 坐标仍为 1-based inclusive。

## Dependencies / parameters / usage

从仓库根目录使用现有 MMseqs2 CLI（依赖及参数详见 [引擎文档](../docs/mmseqs_engine.md)）：

```bash
uv run sprb-scanner mmseqs-search \
  --query reference_data/RN28.fasta \
  --prepared /path/to/FJ-prepared \
  --output /path/to/new-RN28-search \
  --sensitivity 7.5 --evalue 0.001 --threads 2
```

也可用 `--database /path/to/managed-mmseqs-db` 复用数据库。
`RN28.fasta` 是推荐主 query；通用 CLI 仍需显式传 `--query`，不会自动搜索完整 panel。
这些文件作为仓库/source distribution 数据分发，未新增 wheel 内置资源定位 API。
query coverage 和整条 target coverage 分别保留，不因片段的 target coverage 低而自动过滤。
本次仅携带序列参考，没有 RN28 专属评分、功能分类或 HMM 融合。

## Error handling / tests

空/错误 query、重复 ID、损坏数据库和后端失败沿用 MMseqs2 Engine 的明确错误状态。
数据更新时同步重算 manifest 的 `output_sha256`（不对 manifest 自身作循环校验）。
`tests/unit/test_reference_data.py` 离线检查所有文件/序列 hash、74 条参考的原始 ID 与
分组关系、parent 坐标、207 aa 主 query 与 73 条 reference_only 的分离。

此前对 FJ 全部 5,139 条蛋白的 RN28 单 query 实测恢复 SprB 精确自匹配 [681,888)：
identity/query coverage 均为 100%，target coverage 为 3.1861%，E-value 为 3.359e-142，
bit score 为 439（MMseqs2 18-8cc5c；sensitivity 7.5；E-value 0.001；2 threads）。
该次搜索仅返回这一个 alignment。这是输入/匹配验证，不构成功能或跨物种特异性验证。
搜索结果、数据库和完整基因组保留在仓库外。
