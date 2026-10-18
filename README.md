# SprB-like Scanner V2

A modular and scalable bioinformatics framework for identifying SprB-like proteins using
sequence homology, protein domain architecture, genomic context and structural evidence.

当前为 **Phase 1 / 第一批 Project Foundation**：可安装的 Python 包、统一数据模型、
坐标工具、严格配置校验和 CI。尚未实现基因组解析或 BCTD 检索，不能用于生物学结论。

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

`validate` 目前只检查配置格式和参数，不验证基因组、模型文件或扫描结果。
默认数据和结果目录位于 `~/.local/share/sprb-scanner/`，校验命令不会创建目录。
相对路径以配置文件所在目录为基准。未知配置项直接报错，避免参数拼写错误被忽略。
`bctd` 中阈值仅为待验证的探索性参数；第一批不执行扫描，也不据此判定 SprB-like。

## 开发流程

功能分支 → PR → CI → Tsutayaaa 本人确认 → Merge Commit。
不得直接提交 main，不得 squash/rebase 合并。暂不设置许可证。
仓库端保护规则是否已启用，见 [开发治理](docs/governance.md)。

```bash
uv run pre-commit install
```

## 文档

- [架构与未来接口](docs/architecture.md)
- [数据模式](docs/data_schema.md)
- [坐标规范](docs/coordinate_conventions.md)
- [开发路线与验收](docs/roadmap.md)
- [模型数据约定](profiles/README.md)

依赖和锁文件管理遵循 [uv 项目文档](https://docs.astral.sh/uv/guides/projects/)。
模型验证使用 [Pydantic validators](https://docs.pydantic.dev/latest/concepts/validators/)。
