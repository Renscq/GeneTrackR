# GeneTrackR 0.99.10 P0 阶段验收

2026-11-21。基线是 `devel` 的 0.99.9 提交 `808ab14`。本轮完成可选的多父级注释视图。 [初审意见](https://github.com/Bioconductor/BiocContributions/issues/177#issuecomment-6853763197) · [0.99.9 进度](review-progress-0.99.9.md)

| 项目 | 结果 |
|---|---|
| 转录本视图 | `as_granges(x, level = "transcript_features")` 返回命名 `GRangesList`，包含 exon 和 CDS。`S4Vectors::mcols(view)$source_rows` 将每组成员对应到当前输入的平面 Feature 行。 |
| 去重与保真 | 同一来源行在每个转录本内只关联一次。保留同坐标的独立记录、分段子特征及原始元数据。同一 ID 的分段转录本归入同一组。 |
| GFF3 和区域 | 支持常见 RNA 类型、核心 SO 编号和百分号转义的大小写匹配。未知父级、缺失分组信息及歧义标识有明确报错。区间裁剪与合并重命名已通过回归。 |
| 旧对象兼容 | 用实际 0.99.7 生成并保存的两份 GFF 对象验证读取和分组，包括原始 `parent_id` 与另一转录本展示 ID 同名的情况。基因、转录本与外显子统计与原版本结果一致。 |
| TxDb 对照 | 8 条输入记录均保留。按显式 `Name` 与 `ID` 映射，3 条转录本声明及 4 条成员关系全部吻合。新视图精确匹配 3 条 exon 与 1 条 CDS 关系，缺失和多余均为 0。 [关系表](io-evaluation/0.99.10-txdb-rows.csv) · [数量](io-evaluation/0.99.10-txdb-counts.csv) |
| 测试 | 1362 条包断言与 121 条工具断言通过，0 失败、0 跳过。包含本轮 65 条新断言。 |
| 覆盖率 | Windows 上同口径 tests-only 覆盖率由 82.59% 升至 82.77%（9769/11803）。新视图实现文件覆盖 79/80 行。 [汇总](coverage/0.99.10-summary.csv) · [逐文件对照](coverage/0.99.9-to-0.99.10-files.csv) |
| `R CMD check` | `Status: OK`，6 篇教程重建及 PDF 手册通过。归档与当前源码、帮助、示例和测试的 175 个文件逐项校验一致。 |
| BiocCheck | 发布验收用根目录的同一归档完整重跑，得到 0 ERROR、0 WARNING、6 NOTE，邮箱注册和 `genetrackr` 关注标签均已核验。此前完整检查曾因 Support Site 邮箱查询超时或 HTTP 429 产生 WARNING。Git clone 检查为 0 ERROR、0 WARNING、0 NOTE。 |

根目录 `GeneTrackR_0.99.10.tar.gz` 与检查归档 SHA256 一致，均为 `3BF352B0EF26B5E1778D96A23A596AED11E74ECD1AD4C3761D5049C164951576`。

## 修改文件

- `R/as_granges.R` 和 `R/transcript_granges.R` 提供选择入口与分组实现。
- `tests/testthat/test-transcript-granges.R` 覆盖共享、重复、分段、错误输入、区域、重命名与旧对象。
- `man/as_granges.Rd` 和 `vignettes/annotation.Rmd` 提供可执行示例及输出说明。
- `DESCRIPTION` 与 `NEWS.md` 将包版本与变更说明同步到 0.99.10。
- `tools/assess_standard_backends.R` 和 `tools/io-evaluation/0.99.10-txdb-*` 保存关系对照与软件版本。`tools/coverage/` 保存新基线、来源校验值及与 0.99.9 的比较。
- [初审逐条回复草稿](review-response-0.99.10.md) 汇总修改与验收证据，`fnd` 记为不适用。

验证在 Windows、R 4.6.1 上完成。多平台评估及 0.99.10 的远程检查待后续执行。
