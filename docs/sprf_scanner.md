# PorP/SprF-like 家族双模型扫描

## 模块合同

- Purpose：提供 PorP/SprF-like 家族证据，复用 `scanners.hmm.run_hmm`，不重新实现搜索算法。
- Input：已通过 prepare/validate 的单组装蛋白数据；TIGR03519 与 PF11751 各一个
  ProfileConfig JSON。HMM 文件路径相对于各自 manifest，或使用绝对路径。
- Output：每模型完整 HMM 输出、合并 hits.json、protein_evidence.json/tsv、summary.json、
  run_metadata.json。保留每个 domain 的 ID、坐标、独立来源和全部原始统计量。
- Dependencies：现有 PyHMMER Engine、core 模型和导出工具，无新增依赖。
- Parameters：`--threads` 默认1，`--sequence-evalue`/`--domain-evalue` 默认10，
  `--max-sensitivity` 默认false；含义与 hmm-search 相同，随机种子由引擎固定为42。
- Error Handling：重复/错误模型ID、重复蛋白ID、输入损坏在运行前报错。单模型加载/搜索
  失败仍运行另一模型，保留各模型 failed 状态和原因；整体 partial/failed 时 CLI 退出2。
  两模型都成功（包括零命中）为 succeeded、退出0。拒绝覆盖输出目录。
  IO故障等无法可靠导出时直接报错，不能保证留下完整状态目录。
- Test Cases：tests/unit/test_sprf.py 使用临时生成HMM，验证证据独立、GA差异、多domain、
  部分/全部失败、空结果、重复ID、相对路径与CLI状态；不联网。

## 官方模型与来源

| 模型 | 本次固定版本 | 模型长度 | sequence/domain GA |
|---|---|---:|---|
| TIGR03519 | TIGR03519.1 | 292 | 90.55 / 90.55 bits |
| PF11751 | PF11751.15 | 273 | 22.3 / 22.3 bits |

NCBI模型描述为 PorP/SprF family type IX secretion system membrane protein；Pfam模型
描述为 Type IX secretion system membrane protein PorP/SprF。每个模型的实际版本和
SHA-256 固定在 profiles/TIGR03519.profile.json、profiles/PF11751.profile.json。
模型文件放在仓库外。NCBI current 和 InterPro API 都是可变来源；更换文件需明确审查
版本与校验值。Pfam API 返回 gzip，manifest 的校验值对应解压后的HMM文件。
本模块不重新分发模型文件，也不把NCBI CDD的PSSM替代为Pfam HMM。

来源：[NCBI模型](https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/TIGR03519.1.HMM)、
[InterPro/Pfam模型](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF11751/?annotation=hmm)、
[Pfam下载说明](https://pfam-docs.readthedocs.io/en/latest/pfam-entry.html)。

## 使用方法

```bash
.tools/bin/uv run --locked sprb-scanner sprf-search \
  --prepared "$HOME/.local/share/sprb-scanner/results/FJ-batch2-acceptance" \
  --tigr-profile "$HOME/.local/share/sprb-scanner/models/porp-sprf/TIGR03519.profile.json" \
  --pfam-profile "$HOME/.local/share/sprb-scanner/models/porp-sprf/PF11751.profile.json" \
  --output "$HOME/.local/share/sprb-scanner/results/FJ-sprf-$(date +%Y%m%d-%H%M%S)"
```

其他机器将两个 HMM 文件与对应 profiles/*.profile.json 放在同一目录；也可修改 manifest
中的 path，不能修改来源身份和校验值来掩盖模型不匹配。输入可以是基因组注释准备结果，
也可以是 protein_fasta 准备结果；后者没有基因组邻域信息。

单模型仍可使用已有 hmm-search 加任一 manifest，无需双模型入口。

## 证据语义

- `hit_models` / `support`：哪些模型产生了后端命中，不要求通过GA。
- `ga_models` / `ga_support`：哪些模型至少有一个 domain 通过自身 sequence/domain 双GA。
- `ga_by_model`：该蛋白在命中模型上的 passed_ga/below_ga/unavailable 状态。
- `model_states`：两个模型的成功/失败状态。一个模型失败时 support 会明确使用
  observed_other_failed，不能解释成另一模型的阴性。ga_support=none仅表示没有观察到
  通过GA的证据；须结合模型状态和原始命中判断。
- 不合并/平均分数，不要求两个模型区间一致，不删除低于GA、局部或多domain命中。
  combined hits.json 中 raw_artifact 路径指向对应模型子目录。

PorP/SprF-like 是家族归属证据；不是特定 Type B CTD 底物的配对结论。模型可能识别
较长蛋白中的局部家族区域（如SprD），因此不加全长/覆盖率硬过滤。邻域、膜拓扑和
具体底物关联留给后续模块。

## FJ真实数据验收

本地完整 FJ GCF_000016645.1（RS_2025_12_11）的5,139条蛋白：两个模型各有11条
蛋白/11个domain，全部通过各自GA。合并保留11条蛋白、22条独立domain证据。

- SprF：Fjoh_0978 → WP_044048232.1，324 aa，两模型通过GA。
- orphan shuttle：Fjoh_1677 → WP_012023753.1，329 aa，两模型通过GA。
- SprD：Fjoh_0980 → WP_012023066.1，1,588 aa，两模型通过GA；仅作家族区域证据。

核查了蛋白ID、区间、独立domain E-value、原始文件引用与全套输出SHA-256。
结果在仓库外 `/private/tmp/sprb-sprf-validation/FJ`。此为实测记录，不作为跨物种准确率证明。
