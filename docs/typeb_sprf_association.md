# Batch 4B：Type B CTD / PorP-SprF 邻域关联

## 模块合同

| 项目 | 约定 |
|---|---|
| Purpose | 将既有HMM证据与客观CDS邻域关联，允许一个cargo对应多个家族命中蛋白 |
| Input | 同一基因组prepared目录、完整bctd-hmm-search与sprf-search输出目录 |
| Output | associations.json、associations.tsv、domain_evidence.tsv、summary.json、run_metadata.json |
| Dependencies | core记录、NeighborhoodIndex、Pydantic、标准库；不调用HMM搜索 |
| Parameters | window默认5；候选为任一CTD模型的所有保留命中，不加GA或coverage筛选 |
| Error Handling | 身份/序列摘要不一致、文件损坏、缺模型/失败扫描、缺来源摘要、非法坐标/覆盖率明确失败；输出不覆盖 |
| Tests | tests/unit/test_association.py；tests/regression/test_fj_association.py（完整FJ可选本地回归） |

## 命令

```bash
uv run sprb-scanner typeb-sprf-associate \
  --prepared /path/to/FJ-prepared \
  --bctd-scan /path/to/FJ-bctd \
  --sprf-scan /path/to/FJ-sprf \
  --window 5 --output /path/to/new-association
```

两个扫描目录分别来自TIGR04131/PF13585与TIGR03519/PF11751独立双模型命令。
四个模型必须成功完成（成功但0命中允许）。不接收单独hits.tsv、缺模型或partial扫描，
避免把失败/未执行当成无局部证据。prepared必须含CDS坐标，protein-only输入不适用。

校验每个扫描目录与子目录的输出清单、SHA-256、模型版本/来源/模型SHA-256、子模型与
汇总hits的一致性，以及完整蛋白ID集合、长度、序列SHA-256。摘要按扫描记录的蛋白顺序
重建，prepared中的存储顺序变化不影响比对；不接受不同序列、部分蛋白组或缺失历史摘要。
清单路径不能越过输入根目录。读取后及导出前再次检查输入哈希，临时目录原子导出。
这些是输入一致性检查，不提供对人为重写全部清单的密码学身份认证。

旧扫描若没有input_proteins_sha256会明确拒绝：不得猜测身份或伪造摘要，需使用带摘要的
已有完整扫描，或另外显式运行现有扫描命令；关联模块不会自动重扫。
HMM摘要证明的是蛋白ID与序列匹配；输出坐标来自本次prepared，历史HMM导出没有记录
完整基因组坐标摘要，不能凭HMM摘要证明历史注释坐标逐字节相同。

## 候选与关系语义

- CTD候选=任一TIGR04131/PF13585保留命中的蛋白并集，包括仅below_ga的蛋白。
  hit_models/ga_models与每条domain的原始GA状态独立保留；这不是最终Type B或SprB分类。
- 对每个cargo调用NeighborhoodIndex查询±window个CDS排名，不跨contig，按既有规则处理
  负链、compound、重叠和显式环状回绕。伪基因/partial/无翻译CDS仍计入CDS排名。
- 邻域中任一SprF模型保留命中的蛋白都关联，不限于直接邻接或同向，不作一对一选择。
  一个邻居上的多个domain/model通过domain_hit_ids引用，不合并分数或覆盖率。
- status=local_hit表示存在局部保留证据；no_local_hit仅表示当前窗口和已有扫描没有
  局部保留命中。它不是orphan判定，也不能排除远端shuttle或未被原扫描保留的弱同源。
- 关联只是位置关系，不产生cognate配对、典型全长shuttle、假阳性、operon或评分标签。
  没有MMseqs2、RN28或搜索阈值变更。

所有坐标为0-based half-open。CDS排名距离、intervening_cds、distance_bp（非负）、
overlap_bp、signed_gap_bp、双方strand、相对cargo转录方向、same_strand、直接邻接及
跨原点标志沿用[邻域模块合同](neighborhood.md)。有符号gap不表示上下游。

## Domain位置与coverage

完整DomainHit原样嵌入证据记录，分别保留sequence E-value/bit score、independent/
conditional domain E-value、domain bit score、GA cutoff/margin/状态和模型版本/来源/哈希。

- hmm_coverage=query_coverage=(query_end-query_start)/query_length。
- protein_coverage=target_coverage=(end-start)/protein_length。
- alignment坐标、HMM模型坐标、envelope坐标分别保留；覆盖率是alignment跨度比例，
  不是sequence identity或剔除gap后的逐残基覆盖率。
- domain_start/end_fraction及N/C端剩余aa给出蛋白内位置；domain_order_by_protein_start
  仅表示同蛋白、同模型按alignment起点排序的序号，不冒充HMMER原始编号。
- 多domain逐条保留并列出overlapping_same_model_hit_ids；不简单累加重叠区间，也不将
  两个模型的覆盖率相加。无任何coverage cutoff或N端/内部/C端的任意硬分类阈值。

FJ SprF（Fjoh_0978）TIGR03519/PF11751的HMM coverage为100%/99.63%，
protein coverage为95.06%/88.89%；SprD（Fjoh_0980）分别为100%/96.70%及18.95%/16.56%。
SprD匹配接近完整的家族模型，但alignment局限于1588 aa蛋白N端[0,301)或[36,299)。
该证据不能自动等同于典型全长SprF，也不自动判定为假阳性。两条均通过GA且均保留。

## 输出与复现

- associations.json：版本化结果，cargo及一对多邻居、完整基因注释/分段和关系；domains
  是全部CTD与SprF原始证据库，包括没有参与任何局部关系的远端命中，通过hit_id引用。
- associations.tsv：每条cargo–邻居蛋白一行；no_local_hit cargo也单独一行，邻居字段空白。
  两侧domain证据以JSON列表列完整保存，不做两个模型间domain的笛卡尔积。
- domain_evidence.tsv：每个domain一行，模型/GA/raw_hit等复合列为JSON；包含所有四模型
  保留命中、蛋白ID、完整基因组分段、coverage和位置，避免因局部关联而丢失原证据。
- summary.json：全部候选数、local_hit/no_local_hit数、关系数和domain数。
- run_metadata.json：版本、参数、输入文件绝对路径/SHA-256、原始prepared输入来源及
  扫描元数据（含模型来源/版本/模型SHA-256）、输出SHA-256。raw_artifact路径以每条
  domain的scan_directory为基准；原始文件保留在输入目录，不复制大型数据。

JSON不按参考蛋白名单筛选，FJ候选可多于12个：12个为文献参考集合，额外below-GA候选
如实保留。回归测试只对12个参考cargo核对9个有邻接证据、3个无局部命中的分组，
不把额外候选当成已验证阳性。小型参考清单在tests/fixtures/fj_association_reference.json。

默认测试不联网、不下载模型。完整FJ回归复用已有数据（变量都必须设置）：

```bash
SPRB_FJ_PREPARED=/path/to/FJ-prepared \
SPRB_FJ_BCTD_SCAN=/path/to/FJ-bctd \
SPRB_FJ_SPRF_SCAN=/path/to/FJ-sprf \
uv run pytest tests/regression/test_fj_association.py
```

完整FJ/扫描文件不入库；未提供本地路径时该测试skip。CI离线专项覆盖一对多、below-GA、
重叠多domain、coverage、成功空扫描、失败/损坏/错基因组输入、CLI和跨contig/环状/compound。
算法构建一次NeighborhoodIndex，候选查询复用索引；同蛋白/模型m个domain的交叠标注
为O(m²)，不是全基因组两两比较。没有扩大NeighborhoodIndex的复杂度承诺。

## 本批本地验收（2026-10-19）

从最新main e22766b创建codex/typeb-sprf-association。完整本地FJ及四模型既有结果下，
137项pytest全部通过，Ruff lint/format及sdist/wheel构建通过。
真实CLI导出：25个全保留CTD候选、10个有局部命中、15个no_local_hit、12条关联、
109条原始domain。12个文献参考cargo子集恢复9个邻接型和3个无局部命中；
其余13个候选不是文献验证阳性。SprB–SprF距离43 bp，Fjoh_3952–Fjoh_3951为61 bp。
SprD与SprF同时保留，coverage和位置显著不同；输出不产生功能或真假阳性标签。
完整临时输出位于/private/tmp/sprb-association-batch4b-FJ-20261019，不随Git提交。
