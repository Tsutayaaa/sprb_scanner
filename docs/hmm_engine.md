# 通用 HMM Engine 与 TIGR04131

## 范围与模块合同

- Purpose：单个蛋白 HMM 检索，保留引擎返回的全部 domain 证据；不进行最终分类。
- Input：`ProteinRecord` 序列集、单模型 HMM 文件、含来源/版本/SHA-256 的 ProfileConfig。
  CLI 接受已通过 prepare 的数据目录，首先校验完整性及质量错误。
- Output：hits.json、hits.tsv、raw_domains.json、targets.tbl、domains.domtbl、status.json、
  run_metadata.json。成功无命中输出空列表及带表头 TSV；空蛋白集也可成功返回零命中。
- Dependencies：PyHMMER 0.12.3（uv.lock）、core 记录及公共导出工具。
- Parameters：sequence/domain 报告 E-value 默认各 10；线程默认 1；随机种子 42；
  `--max-sensitivity` 关闭 bias filter 并将 F1/F2/F3 设为 1，保留 null2 修正。
- Error Handling：错误路径、校验值/身份/版本不符、错误/空/多模型 HMM、非蛋白 HMM
  均报错；引擎失败在新输出目录记录 failed 状态并使 CLI 非零退出。输入准备目录或
  参数校验失败在启动前报错；输出不可写时无法保证生成状态文件。拒绝覆盖已有目录。
- Test Cases：tests/unit/test_hmm.py 覆盖上述失败、零命中、多 domain、未报告 domain 保留、
  与原生 PyHMMER 逐字段对照、输出校验值、重复内部 ID 和后端故障。

## 统计量与坐标

`DomainHit.evalue` 固定为独立 domain E-value，`evalue_kind=independent_domain`；
`bit_score` 为 domain score。`hmm` 字段另存 sequence_evalue、sequence_bit_score、
independent_domain_evalue、conditional_domain_evalue、envelope，以及 sequence/domain
的 reported/included 标志。不得把 sequence E-value 当作 domain E-value。
旧 DomainHit 可省略 hmm 字段；本引擎生成的每条记录均包含它。

HMMER 的 alignment 与 envelope 坐标为 1-based inclusive；规范化记录使用 0-based
half-open，即 start=from-1，end=to。query 表示 HMM 模型，target 表示蛋白序列，
不将蛋白区间冒充基因组坐标。query_coverage 是模型跨度比例，target_coverage 是蛋白
跨度比例；两者仅记录，不过滤。ctd_distance=protein_length-end。

“全部命中”指本次后端返回的全部 domain，包含未报告的 domain，而非所有可能的弱相似性。
HMMER 的 sequence 报告阈值和启发式过滤仍可能排除目标；不会进一步按长度、覆盖率、
末端距离或最佳 domain 删选。默认探索性阈值不是生物学分类标准，也不自动使用 GA/TC。
原生 tbl/domtbl 只包含 HMMER 报告的记录，可能少于 hits.json/raw_domains.json。
原生表使用安全内部 token（p0 等），run_metadata.json 提供 token 到蛋白 ID 的映射。
run_metadata 还记录 Z/domZ、完整蛋白 ID/序列摘要、模型来源、配置与输出校验值。

## 实际使用

在项目目录运行。以下模型/数据已经存在于本机，无需再次下载：

```bash
.tools/bin/uv run --locked sprb-scanner hmm-search \
  --prepared "$HOME/.local/share/sprb-scanner/results/FJ-batch2-acceptance" \
  --profile "$HOME/.local/share/sprb-scanner/models/TIGR04131.1/profile.json" \
  --output "$HOME/.local/share/sprb-scanner/results/FJ-hmm-$(date +%Y%m%d-%H%M%S)"
```

其他机器先通过 prepare 生成数据目录，再从 manifest 的 source URL 下载模型，将
profiles/TIGR04131.profile.json 与 TIGR04131.1.HMM 放在同一目录，用 --profile 指向 manifest。
模型校验失败时需核对版本，不应直接修改 checksum 使未知模型通过检查。

## 已完成的真实数据验收

FJ GCF_000016645.1、RS_2025_12_11 的完整 GBFF 蛋白集：5,139 条记录；模型
TIGR04131.1，NCBI PGAP HMM release 20.0，长度 85，描述 T9SS type B sorting domain。
模型 SHA-256：`501e3ad25dace893fdd99d68c2dde1643ca9b2f40a051c6a8bbff36fbcc564fe`。
在 PyHMMER 0.12.3、默认设置下得到 15 条蛋白上的 25 条 domain，其中 23 条 domain
被报告，2 条未被报告但保留。完整结果在仓库外 `FJ-TIGR04131-engine-v1` 目录。

已知 SprB（WP_012023065.1，6,497 aa）被识别：

| 指标 | 结果 |
|---|---|
| 蛋白 alignment，0-based half-open | [6413, 6497) |
| 对应 1-based inclusive | 6414–6497 aa |
| HMM 区间，0-based half-open | [3, 85) |
| Domain bit score | 87.21968078613281 |
| Sequence E-value | 5.752534390697332e-26 |
| Independent domain E-value | 1.6160409413097677e-25 |
| Conditional domain E-value | 4.716990488353088e-28 |
| C 端距离 | 0 aa |

这是一次真实数据验收记录，不新增生物学回归测试，也不代表完成检出率/假阳性率评估。
本批不实现 MMseqs2、自定义 CTD 模型、证据融合或最终 SprB-like 分类。

依据：[PyHMMER 统计量](https://pyhmmer.readthedocs.io/en/latest/guide/significance.html)、
[结果字段](https://pyhmmer.readthedocs.io/en/stable/api/plan7/results.html)、
[NCBI 模型目录](https://ftp.ncbi.nlm.nih.gov/hmm/current/)。


## GA 证据分层

每条新命中在 `hmm.ga` 保存 `status`、模型的 sequence/domain GA cutoff、各自的
score-cutoff 差值及通过标志。阈值从校验后的模型读取，不固定为35，也不替代原始分数。

- `passed_ga`：sequence score 和该 domain score 均 ≥对应 GA（包含等于边界）。
- `below_ga`：已有后端 domain 命中，但至少一项低于 GA。弱命中和未报告命中同样保留；
  此状态仅表示比较结果，不自动确认为 Type B CTD，也不证明有分泌功能。
- `unavailable`：模型没有 GA，cutoff/差值/通过标志为 null，不推导为失败或通过。

旧记录的 `hmm.ga` 可以缺省（null），表示尚未计算；这与新记录的 unavailable 不同。
原始 raw_domains.json、分数、三种 E-value、坐标、reported/included 标志均不改变。

新增 `hit_evidence.tsv` 将 GA 字段、原始分数、E-value、坐标、末端距离、覆盖率和模型
来源标识平铺成列；`hits.json` 与 `hits.tsv` 也包含完整嵌套 GA 证据。
新增 `evidence_summary.json` 给出 domain 和蛋白两个层面的三类计数。蛋白类别互斥：
任何 domain 通过 GA 即计入 passed_ga，否则计入 below_ga；无 GA 则为 unavailable。
零命中时各类计数均为0，ga_available 仍准确反映模型是否有阈值。失败不输出成功汇总。

这些字段用于后续组合文献证据、末端距离、覆盖率和邻域关系；不在当前引擎引入新的
独立 E-value、蛋白长度或末端位置过滤，也不自动降低 GA 阈值。

测试新增双阈值边界、单项不通过、无GA、旧记录读取、低于GA命中保留、互斥蛋白汇总、
阈值来源和错误差值校验。临时 FJ 实测用于核对新输出，不新增生物学回归测试。
