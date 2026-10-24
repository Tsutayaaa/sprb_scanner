# Batch 6C MVP — SLURM arrays

## Purpose

在每个数组任务内先严格 prepare 一个注释基因组，再调用 Batch 6B `run_pipeline(..., resume=True)`。
科学算法、GA、原始参考序列、证据并集及邻域语义不变。不实现分类、评分、清理或归档。

## Input / dependencies

- NCBI `ncbi_dataset/data/`：每个版本化 GCF/GCA 目录恰好一个 `.gff`、`.fna`、`.faa`。
- 现有严格版本化 Pipeline YAML、Registry、本地五个 HMM，MMseqs2 与 PyHMMER。
- Python 3.11+，共享 POSIX 文件系统及 Slurm `sbatch`/`sacct`；不要求管理员权限。
- 独立资源 JSON：partition、account、qos、cpus、memory、time、concurrency、max_array_tasks、
  max_array_size、绝对 Python 路径均必填；modules 是可选的显式模块列表。
- 原始组装目录只读，生成文件和临时文件全部位于 batch 工作区。

`config/slurm.resources.example.json` 是 XJTLU 现场检查后给出的**站点示例**，不是通用默认值。
本次查得 `MaxArraySize=1001`，`52cores` QOS 上限为每用户 52 核、52 个提交作业。
示例一波至多 40 个数组项、并发 6，每任务 4 核；控制进程也占用一个作业和 CPU。
其他账户中的作业同样可能影响额度，提交被拒会保留诊断；不自动提高配额或换账户。

## Preflight / deployment

先检查 `sinfo`、`sacctmgr show assoc/qos`、`scontrol show config`、`mmlsquota`、`module avail`。
不能只根据登录节点的软件推断计算节点环境。通过 `salloc`/`srun` 验证导入、资源 checksum、
MMseqs2 版本及服务接口，小批量实际执行后再扩大运行。登录节点不运行搜索。

建议从 Conda 已提供的合适 Python 建立项目目录内 venv，或使用项目内 Conda prefix。
不要改动既有 base 环境。用 `uv export --locked --no-dev --no-emit-project` 导出带哈希依赖；
从国内镜像安装时仍使用 `pip --require-hashes --only-binary=:all:`。轮子可直接传输，
用 `pip install --no-deps /path/project.whl` 安装。Registry/FASTA 与外部 HMM 另行传输并校验。
记录 lockfile、wheel、模型、Registry、实际安装版本和 MMseqs2 backend 的 checksum/版本。
不自动下载、替换模型，也不在作业执行期间更新共享环境。

## Parameters / commands

所有路径最好使用绝对路径。Pipeline YAML 仍使用 Batch 6B 的路径解析和参数语义。
RN28 旧 pilot 使用 sensitivity 7.5；复现时显式设置，示例科学配置的原有默认仍是 5.7。

```bash
sprb-scanner batch-plan \
  --dataset /input/ncbi_dataset/data \
  --config /project/pipeline-hpc.yaml \
  --resources /project/slurm.resources.json \
  --output /project/batches/run-001

sprb-scanner batch-submit --manifest /project/batches/run-001/manifest.json --dry-run
sprb-scanner batch-submit --manifest /project/batches/run-001/manifest.json --prepare-only
sprb-scanner batch-submit --manifest /project/batches/run-001/manifest.json
sprb-scanner batch-status --manifest /project/batches/run-001/manifest.json
sprb-scanner batch-retry --manifest /project/batches/run-001/manifest.json --indexes 2,5-7
sprb-scanner batch-disk-usage --manifest /project/batches/run-001/manifest.json
sprb-scanner batch-collect --manifest /project/batches/run-001/manifest.json --output /project/tables/snapshot-001
```

`batch-plan` 读取所有原始文件计算 SHA-256，建议也在计算节点运行。`--accession` 可重复用于小样本。
它只做稳定输入清单，**ready 不表示注释已验证**。每个 task 才严格调用现有 GFF3/FNA/FAA
prepare API 并验证 prepared checksums。可用 `--prepare-only` 单独运行准备阶段，也可直接运行完整任务。

资源请求不会覆盖 YAML 中 scanner threads；规划时检查 threads 不超过请求的 CPU。
脚本设置项目内 TMPDIR、OMP_NUM_THREADS 和 OPENBLAS_NUM_THREADS，并加载显式模块。
`batch-task --manifest ... --index N` 必须位于 Slurm allocation；`--local` 仅用于明确的本地测试。
Python API 与各命令同名，见 `sprb_scanner.batch`。

## Immutable planning / resume / retry

- 版本化 assembly accession 字典序排序，索引为固定的 `0..N-1`；重试使用原始稀疏索引。
- `manifest.json` + `manifest.sha256` 不再修改；里面保留三件套大小/哈希、来源注释对象、
  计划 prepared 路径和初始 `preparation_status=not_run`。实际准备状态在 task attempt/status 表中。
- 计划内 YAML 是补齐默认值、解析为绝对路径的快照；Registry 和启用资源固定预期 checksum。
- 输入、资源或实现不兼容会明确失败，需建立新计划。不会默默复用旧数据或修复注释。
- 每次 task 使用新的 `tasks/attempt-NNNN/`，保存状态、错误与 Slurm identity；每基因组独立锁。
- 已验证 prepared 可复用。损坏 prepared 拒绝覆盖；保留现场，使用新的计划重新准备。
- 每基因组 Pipeline 工作区独立，成功 stage 的原生输出不可覆盖；重试委托 Batch 6B 的
  fingerprint/checksum/backend/dependency 校验，不另写搜索或 stage-cache 逻辑。
- `batch-submit` 只选择未运行任务，完整执行时也选择已 prepared 的任务；一次最多一波。
  `batch-retry` 只选择失败、中断、部分或损坏结果，不重新派发成功任务。
- 提交锁防止同时选择同一任务；sbatch 超时/模糊响应可能已被接受，标记 `submission_unknown`，
  保守禁止自动重复提交，需人工核对 Slurm，不能把它解释为未提交。
- `batch-status` 默认轻量检查状态；`--verify` 重新校验输入、prepared 和所有阶段文件，可能较慢。
- 调度使用 `JobID` 保留 `arrayjob_index`，不能用内部数字型 JobIDRaw。Slurm OOM/超时/启动失败
  即使没有 task.json 也记录为失败；调度状态查不到则明确 unknown，不认作成功或零命中。
- 计算中断的 running 记录没有调度信息时保持 unknown，不冒险并发重试。

可在一个小型计算 allocation 中运行 `scripts/hpc/drive_batches.py --manifest ...`，依次提交有界波次，
等待每波结束后继续。加 `--prepare-only` 完成全量准备；不加则准备并运行，自动跳过失败准备数据，
但不自动反复重试错误。重启前先确认旧控制器已退出；它不要求常驻服务或网络。
控制器日志中的失败数不等于生物学阴性数。

## Output / scientific semantics

```text
batch/
  manifest.json + manifest.sha256     # immutable identities and original indexes
  pipeline.yaml                      # effective config snapshot, pinned checksum
  submissions/submission-NNNN/        # job.sh, submission.json, scheduler diagnostics
  logs/                              # stdout/stderr per array job + index
  genomes/GCF_x.y/
    prepared/                        # unchanged native prepare output
    tasks/attempt-NNNN/               # task.json, tmp/, optional error.log
    latest_task.json                  # pointer with checksum
    pipeline/                        # unchanged Batch 6B workspace and native scanners
```

`batch-collect` 总是创建新快照，验证后逐基因组流式导出；不重新搜索、不筛掉 below-GA 域、不合并
重叠 alignment。原生 TSV 字段逐行保留，另外附 assembly/genome/index、source protein ID 和 cargo gene。
`gene`（邻居）与 `cargo_gene` 必须区分，两者均保留分段坐标和 qualifiers。

- `genome_summary.tsv`：全部基因组状态、完整性标记、逐模型命中/GA 蛋白数、Type B/邻域统计与磁盘占用。
- `protein_summary.tsv`：六来源证据并集；`cargo_summary.tsv` 只是其中有任一 retained CTD 域的子集。
- `evidence_records.tsv`：所有原始 HMM domain 和 RN28 alignment，含原始 JSON，不压缩成最佳匹配。
- `domain_summary.tsv` / `rn28_alignments.tsv`：上述记录按证据种类分表，坐标为 0-based half-open。
- `architecture_summary.tsv`：PorP/SprF 两个模型的每条域证据，含 coverage、坐标、GA 与原始域位置。
- `association_summary.tsv`：已有 Merger 输出的邻域关系，完整保留关系和候选 neighbor IDs。
- `source_status.tsv`：逐基因组逐来源状态；无可验证 Merger 输出时用 `unavailable`，不填零命中。
- `provenance.jsonl` / `collection_metadata.json` / `status.json`：注释来源、实际参数/版本/模型身份、
  输入输出 SHA-256、纳入与排除数量、失败详情。

`typeb_counts_complete` / `association_counts_complete` 为 false 时，非空计数只是观测到的证据下限。
缺失来源计数为空；成功且无命中才是 0。`no_local_hit` 不等于 orphan，HMM 阳性不意味着功能性 shuttle。
原始分数不同统计类型不混用，不添加覆盖率 cutoff 或综合评分。

## Disk measurement / error handling

`batch-disk-usage` 和 task.json 保留当前文件数、apparent bytes、allocated bytes (`st_blocks * 512`)；
按 inode 去除硬链接重复，遇到非预期 symlink 明确报错。两类容量可能因 GPFS 分配/压缩而不同，
不声称是历史峰值。collection 自身容量可额外用 `du` 测量。
不新增结果清理、归档或删除机制；既有扫描器内部的临时数据库生命周期保持不变。

严格保留 GFF 翻译不一致、无法支持的重叠 CDS、缺失蛋白等错误，不通过降低 strict_quality 绕过。
CLI task 成功/准备完成返回 0，partial 返回 3，失败返回 2。collect 完成只表示快照导出成功，
不表示全部组装已成功；必须查看 included/excluded 和 complete_evidence_genomes。
无法读取调度状态不会伪造运行结果。无效清单或变更资源全局拒绝，单株错误局部隔离。

## Test cases

`tests/unit/test_batch.py` 覆盖严格资源/不可变清单、固定索引、输入变更/注释错误、准备结果损坏、
进程锁、生成脚本真实执行、调度失败和未知派发防重、缺失证据计数。
`tests/integration/test_batch_pipeline.py` 用真实小型 HMM/MMseqs2 验证六来源、全部 alignment、
原生 TSV 逐行保留、成功阶段复用及坏输出排除。HPC 验收必须记录 allocation/job ID、节点、
全部组装成功/失败清单、实际磁盘占用、环境/模型身份和小批次续跑情况。
