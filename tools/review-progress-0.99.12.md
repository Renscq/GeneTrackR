# GeneTrackR 0.99.12 验收记录

2026-11-21。本地基线为 devel 上的 0.99.11 提交 `a86804e`。
[初审逐条回复草稿](review-response-0.99.12.md)

| 验收 | 实测结果 |
| --- | --- |
| Windows R CMD check | R 4.6.1，`Status: OK`。六篇教程重建及 PDF 手册通过。1,640 个断言通过，测试日志有 4 条 Bioconductor 依赖初始化警告。 |
| BiocCheck 1.49.45 | 同一归档的新包检查为 0 ERROR、0 WARNING、6 NOTE，Git clone 检查为 0 ERROR、0 WARNING、0 NOTE。 |
| Linux 定向回归 | 同一源码包在 Bioconductor 3.23 的 Linux 容器中安装。LD、图形、信号和索引相关 8 个测试文件的 447 个断言通过，0 失败、0 警告、0 跳过。 |
| Tests-only 行覆盖率 | Windows R 4.6.1 和 covr 4.0.0。9,895/11,812 → 10,400/11,816，83.77% → 88.02%。[汇总](coverage/0.99.12-summary.csv) · [逐文件对比](coverage/0.99.11-to-0.99.12-files.csv) |

验收归档的 SHA256 为
`FB33858A6F4124955F71C7FB845429BB3BF562BC79ACF6B0AB80D2BC6BC58351`。

## 修改范围

- `R/bwg_ops.R` 使重命名统一应用到样本、信号数据与染色体元数据。测试包括重名内存轨道、惰性轨道、生成名与已有 ID 冲突，以及两种存储模式混合时的明确报错。`man/merge_bwg.Rd` 记录了适用条件。
- `R/plot_signal.R` 保留具名分组和颜色向量的名称，并用绘图视窗保留与图窗相交的柱形和读码框图元。
- `tests/testthat/test-ld-block.R`、`test-signal-ops-contracts.R`、`test-signal-plots.R` 和 `test-feature-plot-contracts.R` 验证 LD 数值与样本对齐、信号归一化和合并、负链坐标、热图分箱、要素标签及空区间。
- `DESCRIPTION`、`NEWS.md` 与 `vignettes/annotation.Rmd` 同步版本、变更和真实 I/O 证据。`tools/README.md` 说明覆盖率口径，`.github/workflows/test-coverage.yaml` 使用 0.99.11 的逐文件比较基线。版本化结果保存在 `tools/coverage/`。

## 待处理的边界

macOS 的真实 VCF/GFF3 对照仍待实跑，Windows 和 Linux 的真实数据证据来自 0.99.11。[复跑说明](io-evaluation/0.99.11-assessment.md)

在 0.99.12 上另行复现了三个超出本轮修复范围的行为。

- 两个同名的纯惰性信号对象按 `sum` 聚合时，合并报缺少 `sample_id`。
- 显式零库大小用于 `norm_bwg(method = "RPM")` 时，正的原始信号得到 `Inf`。
- 手工构造 `LDTrack()` 并传入未排序的位点与同顺序矩阵时，构造器排序位点元数据，矩阵仍按原顺序保留。通过 `compute_ld_block()` 构造的对象已通过行列对齐回归。
