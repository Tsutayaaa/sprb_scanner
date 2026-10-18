# 输入文件标准 v1.0（Phase 1）

适用范围：单个原核生物组装的 CDS/蛋白数据准备。多 contig 可用；每次处理一个组装。
本阶段不做基因预测，不宣称支持所有 GFF 方言或真核转录本模型。
必须显式指定 `--format`；扩展名用于检查，不用于猜测生物学内容。

| 优先级 | 模式 | 文件及扩展名 | 必需内容 | 可用于后续功能 |
|---|---|---|---|---|
| **首选** | `genbank` | `.gbff`、`.gbk`、`.gb` | 核酸序列、CDS 坐标/方向、尽量完整的 `/translation` | 序列扫描、基因组位置及未来邻域分析 |
| 支持 | `gff3` | `.gff3` 或 `.gff` **加** FNA **加** FAA | 同一组装、同一注释发布的三个文件；明确 ID 映射 | 同上 |
| 仅序列 | `protein_fasta` | `.faa`（推荐）、`.fasta`、`.fa`、`.fas` | 氨基酸 FASTA | 后续序列扫描；无基因座/邻域能力 |

以上文件均支持单文件 gzip，例如 `.gbff.gz`、`.gff3.gz`、`.fna.gz`、`.faa.gz`。
配套核酸 FASTA 扩展名接受 `.fna`（推荐）、`.fasta`、`.fa`、`.fas`。
普通文本采用 UTF-8；FASTA 第一行必须为 `>` 标题，序列记录不得为空。
不接受 ZIP 直接作为 prepare 输入；download 命令负责解压 NCBI Datasets 包。

## 推荐数据来源及版本

优先使用 NCBI RefSeq 的 versioned assembly accession，如 `GCF_000016645.1`。
GenBank `GCA_` 也接受，但必须有完整注释包。下载一个组装时同时获取 GBFF/GFF/FNA/FAA，
保存 NCBI assembly report、注释发布日期和每个文件 SHA-256。
组装版本固定不代表注释版本永久固定；复现应同时固定注释发布和文件校验值。
**不要混用不同发布的 GFF、FNA 和 FAA。** 三个文件应来自同一次下载。
默认推荐直接使用其中的 GBFF：单文件保存序列、CDS 和翻译，减少配对错误。

## GenBank 合同

- 至少一个 CDS；contig ID 不可重复；一个输入只能包含一个物种/组装。
- 保留 join、负链、环状跨原点及 fuzzy location 原始表示；内部每段为 `[start,end)`。
- mixed-strand、remote location、非 join 的复合位置不支持，明确报错。
- `/translation` 是蛋白来源；缺失时保留 CDS 并报告 warning，**不自动推导替代蛋白**。
- 假基因保留在 genes 中，记录原因，不放入可扫描 proteome。
- genes 表在本阶段表示 CDS 记录，不包括仅有 gene feature 而无 CDS 的记录、RNA 基因。

## GFF3 配套文件合同

- 第一行为 `##gff-version 3`；9 列 tab 分隔；属性采用 `key=value`，支持百分号转义。
- GFF seqid 必须精确匹配 FNA 标题的首个空白前 ID；FNA ID 必须唯一。
- CDS 必须有 ID、合法 strand 和 phase（0/1/2）。重复 CDS ID 表示同一 CDS 的分段，
  必须位于同一 contig/strand，关联信息一致；不支持多 Parent 或循环/缺失 Parent。
- Parent 的基因注释会继承，CDS 自身字段优先。保存各分段属性、行号和 phase。
- FAA 匹配顺序：显式 `protein_id` 精确匹配；没有该属性时，CDS ID 或去除 `cds-` 前缀后精确匹配。
  不用描述文本、locus tag 相似性或模糊匹配。歧义直接报错；FAA ID 必须唯一。
- 多个 CDS 可以引用同一 FAA ID：保留独立内部蛋白 ID，并报告重复来源 ID。
- 普通分段按生物学方向排序；显式 `part=X/Y` 时按 part 排序。环状越界坐标要求
  `region` 上有 `Is_circular=true`，跨原点转换为两段。拆成多个原点两侧片段时要求提供
  正确 `part` 顺序，程序不猜测环状连接顺序；重叠片段/复杂移码模型不支持。
- 含 `##FASTA` 的 GFF（如部分注释软件输出）不直接接受，须先拆分；不支持 GTF/GFF2。
- 没有可靠物种名称时 genome.organism 保持 null，不由文件名猜测。

## 序列与质量标准

核酸接受无 gap 的 IUPAC DNA。蛋白接受常规氨基酸及 B/X/Z/J/U/O。
小写转大写、删除一个末端 `*` 都报告 warning；B/X/Z/J 报告 ambiguity warning。
空序列、gap、内部 stop 或非法字符产生 error；原始文件保持不变。
不能仅凭字符完全可靠地区分 DNA 与蛋白，因此用户必须正确选择显式输入模式。

完整且无 exception 的 CDS 按 `transl_table` 校验翻译；未指定时默认表 11（可配置）。
部分 CDS、假基因、未知方向、transl_except/exception 不做简单翻译一致性判断，报告跳过原因。
GFF phase 不一致、缺失 FAA 对应记录、孤立 FAA 记录、翻译不一致均为 error。

`prepare` 默认 `input.strict_quality: true`：有 error 就不导出。
先用 `validate --input ...` 获取完整质量报告。显式设置 `strict_quality: false` 可导出带问题的
数据与报告；不能关闭结构性错误（无效格式、歧义 ID、越界等）。无效蛋白没有可扫描记录，
对应原因始终留在质量报告中。`validate` 遇到质量 error 返回非零状态，即使允许带问题导出。

## ID、坐标与输出

提供 `--genome-id` 可使用 versioned accession；缺省 ID 是输入模式及原始文件 SHA-256 的
完整摘要。相同字节改名不改 ID；压缩方式或注释版本变化会改变缺省 ID。
CDS 按 contig ID、最小坐标、最大坐标及原始顺序确定编号；重复源 protein_id 不合并。
只有相同输入和约定能保证相同 ID；跨注释版本以原始 ID/校验值关联，不能只依赖编号。

输出新目录包含 dataset.json、genome_metadata.json、contigs.tsv、genes.tsv、proteins.tsv、
proteome.faa、quality_report.json、summary.json、run_metadata.json。
TSV 分段坐标列为 JSON；内部/导出均 0-based half-open。FASTA 标题为内部 protein_id。
运行记录包含输入 SHA-256、配置、软件版本、时间及输出校验值；拒绝覆盖已有目录。
Parquet 和 BCTD 命中/候选导出留待第三批，本批不宣称已实现。

## 明确不支持

裸 FNA 直接扫描、FASTQ、BAM/SAM、GTF/GFF2、多物种混合注释、真核多 Parent 转录本、
隐式格式猜测、自动基因预测、缺失翻译自动补全、自动修复错配注释。

参考：[NCBI GFF3](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/)、
[NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/)、
[BioPython SeqFeature](https://biopython.org/docs/latest/api/Bio.SeqFeature.html)。
