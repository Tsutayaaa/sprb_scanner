# 第二批验收记录

验收日期：2026-10-18。本批实现基因组基础设施；不包含 BCTD 扫描或生物学分类验证。

## 本地验收

环境：macOS、Python 3.11.17、uv 0.13.1，依赖由 uv.lock 固定。

- `uv sync --locked`：通过。
- `SPRB_FJ_DATA=<本地完整组装目录> uv run --locked pytest -v`：73 passed。
- `uv run --locked ruff check .`：通过。
- `uv run --locked ruff format --check .`：36 files already formatted。
- `uv run --locked sprb-scanner --help`：通过。
- `uv build`：sdist 和 wheel 构建成功。
- 完整 FJ GBFF 经 CLI prepare 后再 validate --prepared：通过，0 个质量 error。

测试代码见 tests/unit/test_coordinates.py、test_models.py、test_parsers.py、test_ncbi.py，
tests/integration/test_cli_config.py、test_prepare.py，以及 tests/regression/test_fj.py。
下载器测试使用模拟响应，不重复联网下载。完整组装回归仅在显式提供 SPRB_FJ_DATA 时运行；
默认与 CI 使用小型离线片段，完整组装用例明确 skip，预期 72 passed、1 skipped。
GitHub CI 的实际结论以对应 PR 的检查记录为准。

## 真实参考数据

本批共下载 **1 个组装**，未超过授权上限 20 个。
F. johnsoniae UW101：GCF_000016645.1，RefSeq 注释 RS_2025_12_11（PGAP 6.10）。
校验值与预期位置记录在 tests/fixtures/fj_reference.json；下载清单与完整数据留在仓库外。

| 项目 | GBFF / GFF3+FNA+FAA 解析结果 |
|---|---|
| contig | 1 |
| CDS | 5,166 |
| 映射蛋白记录 | 5,139 |
| 假基因 | 27 |
| 部分 CDS | 15 |
| 翻译校验跳过 | 1 |
| 重复来源蛋白 ID 提示 | 1 |
| 质量 error | 0 |

FAA 有 5,138 条记录；两个 CDS 引用同一来源蛋白 ID。保留两个基因组位置，因此注释模式
生成 5,139 条内部蛋白记录。完整回归比较两种注释模式的全部坐标、链、partial/pseudo 标志、
来源 ID 和蛋白序列，并检查与 FAA 序列集合一致。

已核对 sprB（旧 Fjoh_0979）及 sprF（旧 Fjoh_0978）的位置与蛋白长度。
离线小片段为原始区间 [1094604,1115116)，变换方法与来源见 tests/fixtures/README.md。
这验证的是解析一致性，并不证明扫描器的检出率或假阳性率。

本机验收产物：`~/.local/share/sprb-scanner/results/FJ-batch2-acceptance`。
完整下载：`~/.local/share/sprb-scanner/data/GCF_000016645.1`。
输入规范见 [input_formats.md](input_formats.md)，模块合同见 [genome_infrastructure.md](genome_infrastructure.md)。
