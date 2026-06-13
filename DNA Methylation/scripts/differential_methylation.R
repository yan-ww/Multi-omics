#!/usr/bin/env Rscript
"""
差异甲基化区域分析 (methylKit)
"""

library(methylKit)
library(dplyr)

# 读取配置
sample_csv <- "data/samples.csv"
samples_df <- read.csv(sample_csv, stringsAsFactors = FALSE)

# 按条件分组
treatment_samples <- samples_df$sample_name[samples_df$condition == "treatment"]
control_samples <- samples_df$sample_name[samples_df$condition == "control"]

# 读取 CX 报告文件
file.list <- list()
for (sample in samples_df$sample_name) {
    file.list[[sample]] <- paste0("results/03_methylation/", sample, ".CX_report.txt")
}

# 创建 methylKit 对象
myobj <- methRead(
    location = as.character(file.list),
    sample.id = as.list(samples_df$sample_name),
    assembly = "hg19",
    treatment = as.numeric(samples_df$condition == "treatment"),
    context = "CpG",
    min.cov = 10  # 最小覆盖度
)

# 合并样本
meth <- unite(myobj, destrand = FALSE)

# 差异甲基化分析
myDiff <- calculateDiffMeth(meth, 
                            overdispersion = "MN",
                            test = "Chisq",
                            adjust = "BH")

# 筛选显著差异位点
myDiffSig <- getMethylDiff(myDiff, 
                           difference = 25,  # 25% 差异
                           qvalue = 0.05)

# 注释差异位点（需要 TxDb 包）
# 这里简化：直接保存结果
diff_df <- as.data.frame(myDiffSig)
write.csv(diff_df, "results/04_dmr/dmr_summary.csv", row.names = FALSE)

# 保存 R 会话
save.image("results/04_dmr/methylkit_results.RData")

# 生成差异甲基化热图
pdf("results/04_dmr/diff_methylation_heatmap.pdf", width = 8, height = 6)
heatmaps <- function() {
    # 聚类热图
    clusterSamples(meth, dist = "correlation", method = "ward", plot = TRUE)
}
dev.off()

cat("DMR analysis completed. Found", nrow(diff_df), "significant DMRs\n")