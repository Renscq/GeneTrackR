# GeneTrackR 0.99.7 review 进度表

本表区分工程进度与 reviewer 验收；当前版本为本地候选，未提交、推送或同步
BiocStaging。原始 preliminary review：
https://github.com/Bioconductor/BiocContributions/issues/177#issuecomment-6853763197

| 项目 | 工程状态 | 已完成 | 剩余内容 |
|---|---|---|---|
| 依赖清理 | 完成 | 移除 Suggests 中的 roxygen2 | 无 |
| README/API 指引 | 部分完成 | 推荐入口明确，保留兼容别名 | reviewer 是否接受兼容方案 |
| GRanges 与 Feature | 完成 | 双向适配、元数据保真、示例和回归 | reviewer 验收 |
| 教程整合 | 完成 | 六篇编号 BiocStyle 教程，旧网址映射 | reviewer 验收 |
| 工具对比 | 完成 | Gviz/rtracklayer/VariantAnnotation 定位说明 | 性能优劣未作结论 |
| 重复汇总修复 | 完成 | 每转录本/特征一条汇总，双链回归 | 无 |
| 测试与覆盖率 | 本轮完成 | 1144 断言通过；80.17% 升至 82.05% | 信号模块等遗漏分支继续补测 |
| 标准读取器对照 | 首轮完成 | 八组比较、完整记录核对、字段差异、耗时/对象大小 | 大文件峰值内存、区域索引、其他平台 |
| GFF3 多 Parent | 限制已记录 | 原始 attribute 保留；标准化列差异可复现 | 结构化多父级适配尚未修复 |
| BED12 块信息 | 限制已记录 | 公共帮助和教程说明当前扁平读取 | block/thick 元数据读取尚未实现 |
| 标准 VCF 互通 | 暂缓 | INFO/FORMAT 样例对照，保留现有原始字段 | 完整头信息、标准类双向适配及保留基因型导出 |
| 类与 I/O 体系 | 保留现状 | 取舍说明、GRanges 边界和对照证据 | GenomicFeatures 映射及后端决策 |
| 工具路径提示 | 低优先级未改 | 非法源码内输出路径会被拒绝 | 拒绝前可能遗留空目录 |
| review 逐条回复 | 草稿完成 | tools/review-response-0.99.7.md | 待授权发布及 reviewer 后续意见 |
| 本地打包检查 | 完成 | R CMD check: Status OK；BiocCheck: 0 ERROR / 0 WARNING / 6 NOTE；源码检查 0/0/0 | 保留并解释既有 NOTE |
| 根目录安装包 | 完成 | GeneTrackR_0.99.7.tar.gz 来自完整检查通过的归档 | 仅本地交付，未提交 |
| 远程 CI / BiocStaging | 本轮未执行 | 上一已推送 0.99.6 的 CI 已通过 | 当前候选未推送；Staging 未同步 |

本轮不删除别名、不修改分析算法、不迁移类、不更换默认读取后端。
后续优先解决已确认的元数据保真缺口，再实施标准 VCF/类适配和大型数据后端评估。
ORCID、可选 fnd 及邮件列表订阅仅依据维护者确认的信息维护。
