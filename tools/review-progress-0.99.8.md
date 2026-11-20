# GeneTrackR 0.99.8 review 进度表

本表记录 0.99.8 的工程修改和本地验收。审查意见及验收进度见投稿 issue。
对照原 preliminary review：
https://github.com/Bioconductor/BiocContributions/issues/177#issuecomment-6853763197

| 项目 | 当前状态 | 0.99.8 的变化 / 剩余内容 |
|---|---|---|
| roxygen2 依赖清理 | 已完成 | 沿用既有修改 |
| README/API 指引 | 已完成兼容方案 | 推荐入口明确，保留别名；待 reviewer 确认 |
| GRanges 与 Feature 双向适配 | 已完成 | 新元数据继续参与往返与筛选 |
| 六篇编号 BiocStyle 教程 | 已完成 | annotation 增加可执行保真示例 |
| Gviz/rtracklayer/VariantAnnotation 对比 | 已完成定位说明 | 不作性能优劣结论 |
| GFF3 多 Parent 元数据 | 本轮修复 | parent_ids 保留全列表，parent_id/原主层级不变，重建属性保真 |
| GFF3 多父级分析图 | 尚未实现 | 不扩增原记录或修改分析；有损主模型提取明确拒绝，可走 GRanges 筛选 |
| BED12 块/厚区间/RGB | 本轮修复 | 结构化保留及读写；重复名称不合并，非法或矛盾编辑报错 |
| BED12 转完整基因模型 | 尚未实施 | 保持扁平 Feature，未迁移类或算法 |
| 回归与覆盖率 | 已完成本轮验证 | 1232包断言通过，新增88；115工具断言通过；覆盖率82.05%→82.43% |
| 标准读取器对照 | 本轮更新 | 同一组原始数据，GFF3/BED12 差异消除；其余范围保持透明 |
| 标准 VCF 双向适配/基因型导出 | 尚未实施 | 现有 VCF 行为和站点级写出不变 |
| 后端复用、GenomicFeatures 映射 | 后续阶段 | 不替换默认后端；大数据、区域索引、峰值内存及跨平台待评估 |
| 包版本、NEWS、帮助、示例 | 本轮更新 | 0.99.8；兼容性边界已说明 |
| R CMD check / BiocCheck / 根目录 tar.gz | 已完成 | Status OK；BiocCheck 0 ERROR/0 WARNING/6 NOTE；根目录包与最终检查包SHA256一致 |
| review 逐条回复 | 草稿已准备 | 已附本地检查结果，待发布 |
| 远程 CI | 见 GitHub Actions | 按发布提交核对 [devel 的运行记录](https://github.com/Renscq/GeneTrackR/actions?query=branch%3Adevel) |
| BiocStaging | 待同步 | 投稿仓库更新及官方构建仍待维护者授权 |
| 原工具路径 Minor | 低优先级未改 | 拒绝非法路径前可能留下空目录 |

新增数据列是兼容性变化；公开函数签名、别名、既有分析模型、VCF/信号后端保留。
GFF3 原始属性默认仍优先保留；重新生成属性时，编辑两份父级字段需要一致。
合并重命名同步全部父级 token；统计、绘图和单倍型定位使用原主模型的内部投影，
原对象不变，也不将临时投影当作完整元数据的替代物。
旧覆盖率、I/O 证据、0.99.7 源码快照和根目录旧包均保留。
