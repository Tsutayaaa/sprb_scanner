# SprB-like Scanner V2

A modular and scalable bioinformatics framework for identifying SprB-like proteins using
sequence homology, protein domain architecture, genomic context and structural evidence.

当前包含基因组解析、数据准备、NCBI下载、HMM证据及 **Batch 4A 基因组邻域查询与 Batch 4B HMM 邻域关联**。
已增加通用 HMM Engine 与 TIGR04131 单模型检索；双模型融合及最终 SprB-like 分类尚未实现。
用法及实测结果见 [HMM Engine](docs/hmm_engine.md)。
PorP/SprF-like 家族支持 TIGR03519 + PF11751 独立双模型检索，见 [sprf-search](docs/sprf_scanner.md)。
新增 [SprB repeat 与 Type B CTD 双模型](docs/sprb_pfam_models.md)：`sprb-repeat-search`（PF13573）及
`bctd-hmm-search`（TIGR04131 + PF13585）。

新增 [`typeb-sprf-associate`](docs/typeb_sprf_association.md)：复用四模型已有扫描，保留一对多邻域关系与逐domain coverage证据，不作功能配对预测。

新增通用 [`mmseqs-search` / `mmseqs-createdb`](docs/mmseqs_engine.md)：支持全长/片段 query、数据库复用、全部报告 alignment 和独立 query/target coverage。

统一 [Reference Registry](references/README.md) 管理序列、MSA、HMM 和 annotation catalog。
RN28 主 query（207 aa）使用 `--query references/sprb/RN28.fasta`；其他73条模块作参考。
旧 `reference_data/`、`profiles/` 路径保留固定兼容副本。
可用 `reference-list`、`reference-validate`、`reference-resolve` 查看和校验资源。

新增独立的 [`neighborhood`](docs/neighborhood.md) 命令，查询任意CDS的坐标/转录上下游、间隔与重叠关系。

**首选输入是 NCBI 注释 GBFF。** 也支持 GFF3 + FNA + FAA（同一注释发布），或仅序列的 FAA。
完整支持范围、配套要求与异常策略见 [输入文件标准](docs/input_formats.md)。

## 安装与验收

需要 [uv](https://docs.astral.sh/uv/getting-started/installation/)；项目使用 Python 3.11+。
uv 会按 `.python-version` 准备 Python 3.11。仓库内 `.tools/` 是本机工具环境，不随项目分发。

```bash
uv sync --locked
uv run pytest
uv run ruff check .
uv run ruff format --check .
uv run sprb-scanner --help
uv run sprb-scanner version
uv run sprb-scanner validate --config config/default.yaml
uv build
```

## 数据准备示例

以下命令下载一个参考组装；数据和结果存放在 Git 仓库外。

```bash
uv run sprb-scanner download GCF_000016645.1 --output ~/.local/share/sprb-scanner/data/FJ
uv run sprb-scanner prepare --format genbank --genome-id GCF_000016645.1 \
  --input ~/.local/share/sprb-scanner/data/FJ/ncbi_dataset/data/GCF_000016645.1/genomic.gbff \
  --output ~/.local/share/sprb-scanner/results/FJ-prepared
uv run sprb-scanner validate --prepared ~/.local/share/sprb-scanner/results/FJ-prepared
```

GFF 模式必须传入 `--format gff3 --input annotations.gff3 --fna genome.fna --faa proteins.faa`。
仅蛋白模式用 `--format protein_fasta --input proteins.faa`；两者同样需传入新的 `--output`。
所有路径均可替换为本地文件；支持 gzip。prepare 不覆盖已有目录、不自动生成缺失翻译。

`validate --input FILE --format FORMAT` 输出完整质量报告；GFF 还需 `--fna`、`--faa`。
无输入参数时 validate 仅检查配置。`--config config/default.yaml` 可指定配置，未知字段报错。
相对配置路径以 YAML 所在目录为基准。`input.strict_quality` 默认 true：质量 error 阻止导出。
输出含标准化 FASTA、JSON、TSV、质量报告及校验值；Parquet 和候选表属于第三批。

测试默认不联网；CI 运行小型 FJ 片段回归。完整基因组回归使用现有本地下载：

```bash
SPRB_FJ_DATA=~/.local/share/sprb-scanner/data/FJ/ncbi_dataset/data/GCF_000016645.1 \
  uv run pytest tests/regression/test_fj.py
```

基准固定于 RefSeq RS_2025_12_11，重新下载若注释已更新会因 checksum 不同而失败，需显式审查基准更新。

## 开发流程

功能分支 → PR → CI → Tsutayaaa 本人确认 → Merge Commit。
不得直接提交 main，不得 squash/rebase 合并。暂不设置许可证。
仓库端保护规则是否已启用，见 [开发治理](docs/governance.md)。

```bash
uv run pre-commit install
```

## 文档

- [输入文件标准](docs/input_formats.md)
- [解析模块合同](docs/genome_infrastructure.md)
- [第二批验收](docs/batch2_validation.md)
- [架构与未来接口](docs/architecture.md)
- [数据模式](docs/data_schema.md)
- [坐标规范](docs/coordinate_conventions.md)
- [开发路线与验收](docs/roadmap.md)
- [统一资源与模型使用约定](references/README.md)

依赖和锁文件管理遵循 [uv 项目文档](https://docs.astral.sh/uv/guides/projects/)。
模型验证使用 [Pydantic validators](https://docs.pydantic.dev/latest/concepts/validators/)。
