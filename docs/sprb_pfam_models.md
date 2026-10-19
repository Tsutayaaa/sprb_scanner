# SprB repeat 与 CHU_C HMM 模型

## 模型与来源

| 模型 | 固定版本 | 长度 | sequence GA | domain GA |
|---|---|---:|---:|---:|
| PF13573 SprB repeat | PF13573.13 | 37 | 28 | 20 |
| PF13585 CHU_C | PF13585.13 | 86 | 27 | 27 |

从 InterPro/Pfam 官方接口取得模型，文件内容的 accession、GA 与 SHA-256 均已核对。
模型字节在仓库外，仓库只保存 profiles/PF13573.profile.json 与 PF13585.profile.json。
接口返回 gzip；清单中的 SHA-256 对应解压后的HMM。source URL是可变来源，版本与
校验值固定本次实际使用的文件；不能随意替换checksum而沿用旧模型身份。

来源：[PF13573 HMM](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF13573/?annotation=hmm)、
[PF13585 HMM](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF13585/?annotation=hmm)。

## PF13573：内部重复证据

- Purpose：识别蛋白内部 SprB repeat 的局部匹配，为结构域架构判断提供证据。
- Input：prepare输出目录、PF13573 ProfileConfig、输出新目录。
- Output：现有HMM引擎的所有输出，加 repeat_counts.json/tsv，逐命中蛋白记录全部domain、
  passed_ga/below_ga/unavailable数量和reported数量。不将多个重复合为一个命中。
  无命中蛋白不出现在计数表；空结果为[]和带表头TSV。
- Dependencies：现有PyHMMER Engine，无新搜索算法或依赖。
- Parameters：sequence/domain报告E-value默认各10；threads=1；max-sensitivity=false，
  含义与hmm-search一致。GA从模型读取，sequence与domain阈值可不同。
- Error Handling：错误模型ID、坏HMM、无效路径、失败状态与不覆盖规则沿用引擎。
- Test Cases：tests/unit/test_sprb_models.py 中的重复命中保留、计数/校验值一致、
  空结果、错误模型身份；模型解析和统计量细节复用已有引擎测试。

本机已有数据使用示例：

```bash
.tools/bin/uv run --locked sprb-scanner sprb-repeat-search \
  --prepared "$HOME/.local/share/sprb-scanner/results/FJ-batch2-acceptance" \
  --profile "$HOME/.local/share/sprb-scanner/models/sprb-pfam/PF13573.profile.json" \
  --output "$HOME/.local/share/sprb-scanner/results/FJ-repeats-$(date +%Y%m%d-%H%M%S)"
```

每个HMM domain alignment是一条证据。原始计数包含弱匹配、局部/重叠片段，因此
all_repeat_domains不代表已确认的物理重复单元数；passed_ga_domains也只是模型支持计数。
不添加全蛋白长度或全长覆盖率硬过滤，不根据重复数量自动分类为SprB。

## PF13585 + TIGR04131：Type B CTD独立双模型

- Purpose：用两个既有HMM独立评价Type B CTD，保留互补与不一致证据。
- Input：prepare目录、TIGR04131与PF13585各一个清单、新输出目录。
- Output：每模型的完整输出，加合并 hits.json、protein_evidence.json/tsv、summary.json、
  run_metadata.json；输出约定与sprf-search相同。
- Dependencies：复用scanners.pair调度与scanners.hmm搜索。PorP/SprF也使用相同调度。
- Parameters：与sprf-search一致，两个模型使用相同报告参数，分别使用自身GA。
- Error Handling：单模型失败保留另一模型结果，整体partial/failed时CLI退出2；
  两个成功（含零命中）退出0。未运行/失败不可记成阴性。
- Test Cases：精简测试正确模型对与原始证据路径；调度的部分失败、空结果和CLI行为
  由已有tests/unit/test_sprf.py复核。

```bash
.tools/bin/uv run --locked sprb-scanner bctd-hmm-search \
  --prepared "$HOME/.local/share/sprb-scanner/results/FJ-batch2-acceptance" \
  --tigr-profile "$HOME/.local/share/sprb-scanner/models/TIGR04131.1/profile.json" \
  --pfam-profile "$HOME/.local/share/sprb-scanner/models/sprb-pfam/PF13585.profile.json" \
  --output "$HOME/.local/share/sprb-scanner/results/FJ-bctd-$(date +%Y%m%d-%H%M%S)"
```

其他机器把下载模型与对应仓库清单放在同一目录，或在清单中指定本地绝对路径。
两个模型的原始分数、三种E-value、坐标与GA状态分别保留，不平均分数、不要求区间完全
一致。hit_models/ga_models分别代表已有命中/至少一个domain通过GA；model_states区分
阴性与失败。蛋白级双模型支持不等于已融合的独立生物学证据，也不构成最终分类。

## FJ实战验收

输入：GCF_000016645.1、RS_2025_12_11，完整5,139条蛋白，默认引擎参数。

- PF13573：5条蛋白、120条原始domain，41条通过GA。已知SprB WP_012023065.1
  有70条原始匹配、19条通过GA；弱匹配全部保留，不将70解释为真实重复数。
- PF13585：62条原始domain，12条domain/12条蛋白通过GA。
- TIGR04131：25条原始domain，10条domain/10条蛋白通过GA。
- 两模型合并：25条蛋白有任一后端命中、13条有双模型后端命中；12条有至少一个模型
  GA支持，其中10条有双模型GA支持。原始命中共87条，未删选。
- PF13585通过GA的12条蛋白与2019文献Table S2映射到当前GBFF的集合完全一致。
  Fjoh_1645与Fjoh_3952在PF13585下passed_ga，在TIGR04131下below_ga，差异证据均保留。

临时结果：/private/tmp/sprb-pfam-validation/FJ-repeats 与 FJ-bctd。已核对全部输出校验值、
命中坐标、E-value、重复计数与文献12个映射集合。不新增生物学回归测试。
文献映射来源：[2019 Table S2](https://journals.asm.org/doi/suppl/10.1128/jb.00218-19/suppl_file/jb.00218-19-s0001.pdf)。
此验收不代表跨菌株检出率/假阳性率已确立。
