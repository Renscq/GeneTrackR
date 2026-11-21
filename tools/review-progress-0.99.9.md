# GeneTrackR 0.99.9 review 进度表

2026-11-21。基于个人仓库 `devel` 的 0.99.8 提交 `ea074ad` 完成本地优化。
[评审意见](https://github.com/Bioconductor/BiocContributions/issues/177#issuecomment-6853763197) · [0.99.8 进度](review-progress-0.99.8.md)

| 顺序 | 评审项 | 0.99.9 进度与验收 |
|---|---|---|
| 已交付 | DESCRIPTION、入口整理 | `roxygen2` 已移除。推荐 `Feature()` 和 `read_vcf()`，兼容别名继续可用。别名策略待评审核认，资助角色和 ORCID 待维护者核实。 |
| 已交付 | 教程、工具对比与 GRanges 入口 | 六篇编号 BiocStyle 教程完成，包含 Gviz、rtracklayer 和 VariantAnnotation 对比，以及 `GRanges ↔ Feature` 元数据示例。 |
| 已交付 | GFF3、BED12 元数据 | 保留 GFF3 全部 Parent 和 BED12 块、厚区间及 RGB，旧主层级分析通过回归。 |
| P0 | 标准 VCF 双向适配 | `as_variant_track()`、`as_vcf()` 可选复用 `CollapsedVCF`。重排、筛选、重复 ID、RDS 跨进程恢复和样本选择已验证。 |
| P0 | 基因型导出、示例及兼容性 | 教程演示标准写入及重读，旧八列写出和无可选依赖运行通过。标准写入器可舍入浮点 FORMAT、更新 `fileDate` 并补齐缺失数组项，精确型化值可保存在 RDS。 |
| P0 | 测试、CI 与覆盖率 | 1297 条包断言和 121 条工具断言通过，0 失败、0 跳过。65 条适配断言纳入 CI 强制执行。同口径覆盖率 82.43% → 82.59%（9680/11721）。 |
| P1 | TxDb 关系映射 | 8 条 GFF3 原始记录完成计数。显式映射 `Name` 与 `ID` 后，7 条转录本、外显子和 CDS 关系全部吻合。主层级少 1 条共享外显子关系。 |
| P1 | 索引 VCF 对照 | 50,000 位点、8 样本的 Windows 合成评估中，两条路径按同一 POS 条件匹配 50/50。记录了依赖、时间和进程峰值内存。Linux/macOS 与大型真实数据仍待评估。 |
| P2 | 多父级分析图、BED12 基因模型、替代后端 | 各自的迁移验收标准已记录，等待更广的映射和跨平台证据。 |
| P3 | 检查工具路径校验 | 创建目录前解析路径并拒绝落在源码树内。Windows 相对路径、大小写和旁路目录回归通过，Unix 符号链接用例待 CI 执行。 |
| 验收 | 版本、源码包与本地检查 | `R CMD check` 为 `Status: OK`，6 篇教程重建和 PDF 手册通过。BiocCheck 0 ERROR、0 WARNING、6 NOTE，Git clone 检查 0/0/0。 |
| 发布 | 远程 CI 与 BiocStaging | 0.99.9 本地候选待发布。0.99.8 的远程 CI 证据已纳入逐条回复草稿。 |

## 本轮范围与文件

- 版本、命名空间与站点入口为 `DESCRIPTION`、`NEWS.md`、`NAMESPACE`、`_pkgdown.yml`。
- VCF 源码为 `R/variant_vcf.R`、`R/retrieve_vcf.R`、`R/merge_vcf.R`、`R/hap_variant.R`、`R/as_granges.R`、`R/write_vcf.R`。适配拒绝缺少 `fileformat`、存在歧义空 ALT 元素、来源无效或投影字段被修改的对象。
- 帮助与示例为 `man/as_variant_track.Rd`、`man/write_vcf.Rd`，以及 `vignettes/annotation.Rmd`、`vignettes/haplotype.Rmd`、`vignettes/export-results.Rmd`。
- 回归与 CI 为 `tests/testthat/test-vcf-adapter.R`、`tests/testthat/fixtures/standard-variants.vcf`、`tests/testthat/test-bioconductor-interoperability.R`，以及 `.github/workflows/check-package.yaml`、`.github/workflows/test-coverage.yaml`。
- 开发工具及记录为 `tools/check_package.R`、`tools/tests/test-check-package.R`、`tools/assess_standard_backends.R`、`tools/README.md`、`tools/io-evaluation/README.md`，以及 `tools/coverage/`、`tools/io-evaluation/` 下的 0.99.9 评估证据和逐条回复草稿。
- 根目录 `GeneTrackR_0.99.9.tar.gz` 与最终完整检查的源码包逐字节一致。SHA256 为 `176C25352EDD87E8FAA98219EE7C69CAD8E8B5F299AC5CFD51490DAA98C4BAD9`。

## 本地检查证据

BiocCheck 的 6 条 NOTE 涉及 ORCID、`suppressWarnings` 使用、函数长度、行长、缩进和 Bioc-Devel 邮件列表确认。测试所用 Bioconductor 依赖在 Windows 加载时报告命名空间提示。R CMD check 最终状态为 OK。

[覆盖率汇总](coverage/0.99.9-summary.csv) · [逐文件覆盖率对照](coverage/0.99.8-to-0.99.9-files.csv) · [后端评估](io-evaluation/0.99.9-assessment.md) · [评审回复草稿](review-response-0.99.9.md)
