# Batch 6B — Single-Genome Pipeline Runner

## Purpose

以一个已经 prepare/validate 的基因组为单位，离线编排现有科学模块：

```text
prepared + selected Registry resources
    ├── typeb_ctd  (TIGR04131 + PF13585)
    ├── porp_sprf  (TIGR03519 + PF11751)
    ├── sprb_repeat (PF13573)
    └── rn28       (existing MMseqs2 service)
             └── evidence (Batch 6A merge-evidence + existing NeighborhoodIndex)
```

不修改 HMM/MMseqs2 算法、GA、报告阈值、参考序列、domain/coverage 或邻域语义。
不引入 SLURM、系统发育、HMM 构建、功能分类、评分或可视化。

## Input / dependencies

- 完整、checksum 与质量检查通过、非空的 prepared dataset。
- 严格的版本化 YAML，必填 `schema_version: "1.0.0"` 与 `registry`。
- 本地 Reference Registry，以及启用阶段需要的 HMM/sequence artifact。
  外部 HMM 使用显式 `artifact` 路径；不下载缺失资源、不自动替换版本。
- 现有 PyHMMER、Pydantic、Biopython、PyYAML；启用 RN28 才要求本地 MMseqs2。
- 本地 POSIX 文件系统，使用非阻塞进程文件锁；不需要调度器或网络。
- 启用 neighborhood 时要求 CDS 坐标且 `porp_sprf` 启用。仅蛋白输入可禁用 neighborhood。

## Usage / parameters

先检查并按需编辑 [完整配置示例](../config/pipeline.yaml) 中的本地模型路径：

```bash
uv run sprb-scanner run \
  --prepared /path/FJ-prepared \
  --config config/pipeline.yaml \
  --output /path/FJ-pipeline

# 同一工作区继续：验证并跳过成功阶段；失效、失败或未完成阶段写入新 attempt。
uv run sprb-scanner run \
  --prepared /path/FJ-prepared \
  --config config/pipeline.yaml \
  --output /path/FJ-pipeline --resume
```

Python API：

```python
from sprb_scanner.pipeline import load_pipeline_config, run_pipeline

config = load_pipeline_config("config/pipeline.yaml")
report = run_pipeline("/path/prepared", config, "/path/workspace", resume=True)
print(report["state"])
```

`run_pipeline` 同样接受 YAML 路径。所有 YAML 内的相对路径以 YAML 目录为基准；
Python Config 对象的相对路径以调用 cwd 为基准。MMseqs2 可执行文件名按 PATH 查找，
包含路径分隔符的值按 YAML 路径规则解析，保留 launcher 路径。

未知字段、重复 YAML key、非字符串 key、错误版本、字符串代替数值/布尔值、非法线程/窗口等报错。
所有 scanner 参数直接沿用现有模型和默认值；本模块不设新的科学默认阈值。
示例 sensitivity 默认为现有 MMseqs2 的 5.7；重现此前 RN28 pilot 请明确设为 7.5。
`max_seqs` / `alternative_alignments: null` 沿用引擎的完整 target/最长序列容量策略。

Reference binding 的 key 固定为五个模型及 RN28，可选 `version` 与 `artifact`：

```yaml
references:
  TIGR04131:
    version: TIGR04131.1
    artifact: /local/models/TIGR04131.1.HMM
  RN28:
    version: sprb-legacy-20261019
```

版本是 Registry resource version，不推断最新版；多个版本时必须显式选定。
artifact override 仅用于 Registry 声明为 external 的资源，仍校验预期 checksum。
现有 bundled RN28 继续用其原始 207 aa 和 `[681,888)` 来源边界。
本批不会编辑 Registry 或资源；资源更新必须由用户另行登记版本和 checksum。

## Execution / strict and partial modes

- `strict`：启用阶段失败或仅部分完成后停止，后续阶段明确 `blocked`；不导出伪完整证据。
- `partial`：独立阶段继续，随后尝试 Merger。可验证的原生部分输出保留成功子模型和失败子模型；
  无有效原生输出的失败阶段以明确 `failed` 来源传入 Merger，不降级为 `not_run`。
- 配置禁用阶段标记 `not_run / disabled`，Merger 同步保留缺失来源状态。
- 成功零命中是 `succeeded`，记录数为 0，可以复用；失败计数保持未知。
- Merger 成功执行与证据完整性不同。部分模式下可以生成完整格式的证据文件，
  但 pipeline state 仍为 `partial`。Neighborhood 查询失败同样影响其阶段状态。
- Pipeline `succeeded` 表示所有**启用**的阶段完成；主动禁用来源时，Merger 的
  `evidence_complete` 仍可能为 false。两者均可查看，不能把禁用证据解释为阴性。

CLI 退出码：0 为所有启用阶段成功，3 为有明确部分结果，2 为失败/配置或输入错误。
Python API 对正常完成的 strict/partial 流程返回带状态的报告；无效输入、锁冲突和中断向上传播。

## Output / checkpoints

```text
workspace/
  workspace.json                    # protocol 与 genome ID 绑定
  .pipeline.lock                    # 仅持锁进程可写，不自动删除他人的锁
  latest.json                       # 最近 invocation 路径、checksum 与最终状态
  genomes/<escaped-genome-id>/
    runs/run-0001/run_manifest.json  # 本次完整配置、来源、阶段动作与状态
    runs/run-0001/configuration.yaml # 原始 YAML 的字节快照（YAML 调用时）
    runs/run-0001/effective_config.json # 解析并补全默认值后的绝对路径配置
    runs/run-0002/run_manifest.json
    stages/<stage>/<fingerprint>/attempt-0001/
      execution_manifest.json       # fingerprint payload、依赖、状态与所有 artifact checksum
      stage.log                     # 阶段诊断、异常堆栈
      result/                       # 未改写的原生 scanner/Merger 输出目录
    stages/<stage>/<fingerprint>/attempt-0002/
```

每次调用保存独立 run manifest，阶段前后写 checkpoint；执行中的 manifest 使用 `running`。
通过临时文件、flush/fsync 和原子 replace 写 checkpoint；原生 scanner 继续使用自己的原子导出。
成功 attempt 之后不再写入。失败或 KeyboardInterrupt 保留诊断；进程被终止时遗留的 running attempt
不会在恢复时被当作成功。恢复将建立新 attempt，不覆盖或删除旧目录。

阶段单位是上述四个扫描服务及最终 Merger；双模型阶段内不另行缓存单个模型。
即使一侧成功、另一侧失败，重试也会调用原有双模型服务并写入新的完整 attempt。

## Fingerprints / validated reuse

Fingerprint 为排序 canonical JSON 的 SHA-256，包括 protocol、相关实现、prepared 文件内容、
选中的 Registry entry 与 pinned lineage、artifact SHA-256、有效参数和后端身份。
后端身份包含版本、平台/Python/相关依赖以及 MMseqs2 executable 或 PyHMMER library checksum。
恢复仍会执行廉价版本/能力探测与 checksum 检查，不执行被复用阶段的搜索。

- HMM 阶段只包含其相关 HMM 引用和 HMM 参数，不包含 RN28 或 MMseqs2 参数。
- RN28 参数或引用变化，仅使 RN28 与 evidence 失效。
- prepared 文件内容或元数据变化，使全部启用阶段失效。
- HMM 引用/参数/后端变化，使相应 HMM 阶段及 evidence 失效。
- neighborhood 窗口变化只使 evidence 失效。
- Registry 中无关 entry 变化不强迫重跑扫描；最终 Merger 可因完整 Registry 来源 checksum 改变而重建。

复用同时验证 fingerprint payload、依赖绑定、完整输出文件集合/checksum，以及原生加载器的
状态、蛋白 ID/序列指纹、模型/query 身份、实际参数和 backend。损坏或不兼容输出会记录 cache rejection，
随后执行新 attempt，旧结果保留。

Evidence 的确定性 fingerprint 使用依赖阶段 identity/state；另行绑定具体 dependency attempt 的
manifest checksum。即使相同输入导致一次扫描重建，最终证据也会针对新的已验证来源重新生成，
不会继续指向旧损坏输出。不同工作区的原始路径、时间戳和日志不要求字节相等；
相同有效输入/设置的科学证据应语义等价。

## Error handling / limits

配置/Registry 结构或 prepared 校验失败在执行前报错，不生成成功阶段。
启用资源缺失、checksum 不符、backend 缺失等是阶段失败，有 attempt、日志和来源信息。
并发修改输入或依赖时拒绝继续，不将不同版本的输入拼接成一个结果。
同一工作区只允许一个进程；不同工作区独立。`--resume` 不允许更换 genome ID 或协议不兼容的 workspace。
不导入缺少执行 identity 的任意历史扫描目录作为“已完成 pipeline 阶段”；原始目录仍可直接用于 Merger。
不自动重试无限次，每次 resume 对失败阶段创建一次新 attempt。

## Tests

- `tests/unit/test_pipeline_config.py`：严格 YAML、路径、默认参数、canonical fingerprint、失败状态 adapter。
- `tests/integration/test_pipeline.py`：真实小型 HMM/MMseqs2 的完整流程、resume、语义等价、
  RN28 参数/版本/内容、prepared/backend 变化、损坏缓存、失败与部分子模型、未完成/中断恢复、
  成功零命中、禁用来源、CLI 退出码及独占锁。
- FJ 完整流程与此前五模型/RN28 结果逐条比较，验证恢复和选择性失效；
  21 株的既有 evidence 输入通过 Merger 兼容性复核，不为此重新运行全部搜索。

真实数据、模型文件、attempt 和批量结果均保存在 Git 仓库外。

### 本地验收记录（2026-10-20）

完整 pytest：405 passed（包含本地完整 FJ 回归与真实 MMseqs2）；Ruff、CLI、wheel/sdist 和包内容检查通过。

FJ 验收使用 sequence/domain report E-value 10、seed 42、HMM 2 threads、RN28 sensitivity 7.5、
E-value 0.001、2 threads；五套 HMM 原始 domain 记录与此前归档逐条完全相等：
TIGR04131 25、PF13585 62、TIGR03519 11、PF11751 11、PF13573 120。
RN28 原始 alignment 同样逐条完全相等（1 段）；最终 39 条 Evidence Profile。
9 个邻接型 / 3 个 no_local_hit 参考 cargo 保持一致。

FJ 再次运行时五个阶段均 verified reuse；单独把 sensitivity 改为 6.5 后，
三个 HMM 阶段复用，仅 RN28 与 evidence 执行新 attempt。
21 株既有结果复核：592 条 Profile、2,807 条 HMM domain、17 条 RN28-positive 蛋白 / 18 段 alignment；
3 条 RN28-positive / CTD-negative 蛋白全部保留，未为此重新扫描 21 个蛋白组。
