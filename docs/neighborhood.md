# Batch 4A：Genomic Neighborhood Analyzer

## 模块合同

| 项目 | 约定 |
|---|---|
| Purpose | 任意CDS的客观基因组邻域，不推断功能、operon、cargo–shuttle配对或候选评分 |
| Input | 已有GenomeDataset、GeneRecord、ContigRecord；CLI读取prepared目录和内部gene_id |
| Output | NeighborhoodResult，含完整target、contig和NeighborRecord；JSON/TSV/运行元数据 |
| Dependencies | core记录、标准库、Pydantic与公共导出工具；不导入HMM/scanners |
| Parameters | window默认5，非负整数；窗口0返回空邻居列表但仍包含target |
| Error Handling | 未知gene_id、负/非整数窗口、损坏输入、异常分段顺序及已有输出目录明确报错 |
| Test Cases | tests/unit/test_neighborhood.py；tests/integration/test_neighborhood_cli.py；仓库外完整FJ验收 |

入口：`sprb_scanner.neighborhood.NeighborhoodIndex(dataset)`，建立一次后重复调用
`index.query(gene_id, window=5)`。全部CDS均进入索引，包括假基因、partial和无translation
记录。RNA或只有gene feature而没有CDS的记录不属于本模块的CDS序列。

## 顺序、方向及环状窗口

按contig分组，以坐标递增方向的CDS起点排序，同起点依次按最大分段终点和gene_id
排序。负链location_parts先按生物学提取顺序反转为坐标前进顺序。跨原点的CDS起点在
原点前的高坐标部分，不使用简单min/max包络把它误放到整个contig开头。
保留原始location_parts及qualifiers，不改变核心数据模式。

每个contig都有独立排序数组和gene_id→(contig,position)字典。查询逐级取左右各window
个CDS排名，输出顺序为左1、右1、左2、右2……，在边界截断或环状回绕，不跨contig。
环状短contig每个CDS最多返回一次，永不返回目标自身；同一CDS可从两侧到达时优先较少
CDS步数，等步数时坐标上游优先。该去重规则意味着小环上不一定两侧各返回window个。

- genomic_position：在坐标递增的CDS排序中为upstream（左）或downstream（右）。
- relative_position：相对目标CDS的转录方向；负链上下游翻转，strand=0时为null。
- genomic_cds_offset：左负右正；transcription_cds_offset对负链翻转，未知方向为null。
- intervening_cds：沿选定索引路径的步数减1；directly_adjacent只表示CDS排序中的相邻。
- same_strand：两条链均已知时比较，否则null。邻居自身方向不改变目标的上下游定义。

仅显式标为circular的contig回绕；topology=unknown按线性边界处理，不自行猜测环状。
多分段按已有生物学顺序解释；未知链时使用提供的分段顺序作为坐标前进顺序。
非环状contig中的倒序分段、无明确slippage依据的同一CDS内部重叠、跨原点多次或超过
一圈的跨度明确报错。带注释的slippage按下面规则建立独立几何并集，不修改原始记录。

## 带ribosomal slippage注释的CDS

Batch 4C的8株失败来自48个CDS，原始GBFF均有/ribosomal_slippage及programmed
frameshift注释：44个重叠1 bp、4个重叠2 bp；正负链各24个。注释支持这种表示，
并不代表本模块对生物学机制或组装正确性做了独立实验验证。

允许内部重叠须具有明确证据：无值的ribosomal_slippage qualifier（键存在，值为空），
或exception精确为ribosomal slippage（忽略大小写及首尾空格，对应NCBI GFF3表示）。
仅note含frameshift、普通exception、疑似组装/测序差错本身均不够；非空flag值或冲突
exception仍报错。本次不支持把其他translation discrepancy自动解释成slippage。

验证原始生物学顺序转换后的坐标前进路径；环状最多展开一次、总跨度不超过一圈。
重复、包含/嵌套而不继续前进的分段即使带slippage标记仍报错。不凭观察到的1/2 bp
设硬性重叠长度阈值；其他重叠长度也需相同明确注释与前进路径。这是位置表示支持，
不是判定任意长度的重叠在生物学上合理。

每个合法CDS构建一次排序、去重的genomic interval union（接触/重叠片段合并，真实
内部间隙保留），只缓存于IndexedCDS.parts。几何计算使用并集，重叠碱基只计算一次；
排序anchor仍来自原始坐标前进路径，CDS排名/直接邻接按基因而不是按片段计算。
GeneRecord.location_parts、原始分段顺序、qualifiers及蛋白翻译均原样保留，导出JSON/TSV
中的基因坐标仍是原始分段。两个CDS各自可有slippage，交集使用它们的去重并集。

geometry仍接收非空、排序且内部不重叠的区间元组；索引负责验证注释并建立该并集，
不能把未经验证的重叠分段直接传给geometry。距离及关联元数据的geometry_policy记录
这条规则。原有无内部重叠的线性、compound、环状行为保持不变。

来源：[INSDC feature table](https://www.insdc.org/submitting-standards/feature-table/)定义
无值flag及join表示；[NCBI GFF3文档](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/)
区分slippage与组装差错造成的相似重叠，列出exception=ribosomal slippage。

## 距离与重叠

所有位置均0-based half-open，字段保留完整分段。geometry对各CDS分段并集做精确计算：

- distance_bp始终非负：不相交时为最近分段边缘间隔，接触或编码区重叠时为0。
  可直接用于后续 `distance_bp <= max_distance` 筛选；如需排除重叠，额外检查overlaps。
- overlap_bp始终非负：两个CDS分段并集实际共有的编码碱基总数。接触不算重叠；
  overlaps严格等于overlap_bp>0，不用包络跨度代替实际编码区交集。
- signed_gap_bp = distance_bp - overlap_bp，独立保留带符号间隔：不相交为正间隔，
  接触为0，重叠为负的交集长度。该符号不表示上下游或链方向，不参与默认邻域选择。
- geometry(a,b,L)返回(distance_bp,overlap_bp,distance_crosses_origin)，三项在交换a/b后
  均不变。输入为按坐标排序、内部不重叠、非空的Interval元组；原始生物学分段顺序由
  索引负责转换。距离字段与同名JSON/TSV输出、元数据规则保持一致。
- 复合CDS的内部间隙不算编码区；另一个CDS落在该间隙时不会被误报为编码区重叠。
- 环状物理距离取两个方向中的较短分段边缘间隔，等距时优先不跨原点的路径。
  该最短物理路径与选取邻居的CDS排名路径可以不同，两者分开标记。

index_crosses_origin表示CDS索引路径回绕；distance_crosses_origin表示最短物理间隔跨原点；
gene_crosses_origin与target_crosses_origin分别表示CDS本身跨原点。crosses_origin是这些
标志的并集，不把线性化表示中的负间隔当作真实重叠。
对于普通单段CDS，非重叠距离就是两段最近边缘间隔；负链仍使用同一物理距离。
对称性测试覆盖普通/接触/重叠/嵌套、compound间隙与多段交集、环状最短间隔、跨原点
compound重叠及等距路径；同时校验非负字段和独立有符号字段之间的一致性。

## 复杂度与复用

N为CDS数，分组/排序与ID索引构建为O(N log N)，常规单段CDS的查询为O(k)，不重新
扫描整个基因组。复合位置还需一次预排序分段；若允许任意多分段，额外构建成本为
Σ O(p_i log p_i)，精确查询为O(k + k*p_target + Σp_neighbor)。这是精确处理编码分段
的成本，不是对全基因组的遍历；不能在无限分段输入上宣称严格O(k)。分段数有界的
常见原核CDS仍为O(k)。JSON/TSV序列化同样受返回的分段和注释大小影响。

CLI每次启动需要加载/校验prepared并建立索引；批量查询应复用Python API中的同一个
NeighborhoodIndex实例。索引引用现有冻结记录；与core约定一致，不应修改共享qualifiers。

## CLI及输出

```bash
.tools/bin/uv run --locked sprb-scanner neighborhood \
  --prepared "$HOME/.local/share/sprb-scanner/results/FJ-batch2-acceptance" \
  --gene 'GCF_000016645.1::NC_009441.1::CDS_000999' \
  --window 5 \
  --output "$HOME/.local/share/sprb-scanner/results/FJ-neighborhood-$(date +%Y%m%d-%H%M%S)"
```

必须使用内部gene_id；可在prepared/genes.tsv查找locus_tag或gene_name获得它。

- neighborhood.json：完整target/contig与邻居GeneRecord、方向/排名/距离/原点关系。
- neighborhood.tsv：平铺常用注释与全部关系字段，location_parts为JSON列，0邻居仍有表头。
- run_metadata.json：参数、软件版本、prepared文件SHA-256、原始来源清单、质量问题数量、
  距离约定和输出校验值。输出使用新目录、暂存后原子重命名，不覆盖现有结果。

输入先校验prepared完整性；quality错误计入元数据，不因无效蛋白而丢弃结构有效的CDS。
损坏的基因链接、越界等结构性错误仍被core验证拒绝。protein-only数据无法查询CDS。
终端显示目标、contig、strand、邻居数和每条邻居的简洁关系。

## FJ UW101验收

本地完整GCF_000016645.1、RS_2025_12_11，索引包含全部5,166条CDS；默认窗口返回10条
不重复邻居。已核对：

| 属性 | 实测 |
|---|---|
| Contig | NC_009441.1 |
| SprB旧locus | Fjoh_0979 |
| SprF旧locus | Fjoh_0978 |
| SprB / SprF strand | -1 / -1 |
| SprB坐标 | [1095622,1115116) |
| SprF坐标 | [1094604,1095579) |
| 直接邻接 / intervening CDS | true / 0 |
| 编码区边缘间隔 | 43 bp |
| SprF相对SprB | 转录下游、基因组坐标上游 |

真实结果在仓库外/private/tmp/sprb-neighborhood-validation/FJ，另有verify_fj.py核查
目标/邻居身份、坐标关系、10邻居去重与输出校验值。CI使用已有小型FJ片段作离线验证。
新增tests/unit/test_slippage_neighborhood.py覆盖48条原始位置及合成边界情况；
tests/integration/test_neighborhood_cli.py验证GBFF→prepared→CLI不改变翻译/分段。
本批不添加HMM或MMseqs2扫描、不做功能配对；Batch 4B再使用这些客观关系关联候选。


## Slippage修复验收（追加于PR #5）

完整本地pytest 214项通过，Ruff lint/format与sdist/wheel通过。48条真实位置及合成
1/2/5/19 bp重叠、双向union交集、compound内部间隙、正负链、环状跨原点、异常拒绝
均有回归；GBFF→prepared→CLI验证原始分段、qualifiers与translation不变。

仅用原有prepared和HMM结果重跑Batch 4C，8/8原失败株恢复，21/21成功；此前13株完整
associations.json逐字段不变，FJ仍为9个邻接型+3个no_local_hit，43/61 bp重点位点不变。
全部源文件/模型/扫描的SHA-256不变；没有重扫HMM或更改任何分数、阈值或cargo规则。

原Batch 4C结果保持不变。修复后的21株中278个GA cargo有119个直接邻接GA SprF-like
候选、136个在±5 CDS有GA SprF-like候选。完整结果及注释核查位于仓库外：
`~/.local/share/sprb-scanner/analyses/batch4c-slippage-fix-20261019/`。
其中validation_report.md、rerun_summary.tsv及cds_audit.tsv说明恢复情况；不把旧13株
分母的比例与新21株比例混用。结果不代表功能性cargo–shuttle配对。

## Batch 6D 兼容扩展

除原有显式 slippage 依据外，接受由严格 GFF 解析器生成、绑定原始分段和 strand 的
FAA 完整翻译核对凭据。它仅表示坐标表示可验证，不作 slippage 或组装错误诊断。
区间并集与距离/邻接算法保持不变；详见 [NCBI GFF 兼容合同](ncbi_gff_compatibility.md)。
