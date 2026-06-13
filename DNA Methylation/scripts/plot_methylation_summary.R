#!/usr/bin/env Rscript
"""
生成甲基化分析汇总图
"""

library(ggplot2)
library(dplyr)
library(tidyr)

# 读取统计数据
stats <- read.csv("results/05_reports/methylation_stats.csv")
dmr <- read.csv("results/04_dmr/dmr_summary.csv")

# 1. 比对率柱状图
p1 <- ggplot(stats, aes(x = sample, y = mapping_rate, fill = sample)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Bismark Mapping Efficiency", y = "Mapping Rate (%)") +
    ylim(0, 100)

ggsave("results/05_reports/mapping_efficiency.pdf", p1, width = 6, height = 4)

# 2. 甲基化水平
if ("methylation_level" %in% colnames(stats)) {
    p2 <- ggplot(stats, aes(x = sample, y = methylation_level, fill = sample)) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Global CpG Methylation Level", y = "Methylation (%)") +
        ylim(0, 100)
    
    ggsave("results/05_reports/global_methylation.pdf", p2, width = 6, height = 4)
}

# 3. DMR 火山图
if (nrow(dmr) > 0 && all(c("meth.diff", "qvalue") %in% colnames(dmr))) {
    dmr$log10q <- -log10(dmr$qvalue)
    dmr$significant <- ifelse(dmr$meth.diff > 25 & dmr$qvalue < 0.05, "Hypermethylated",
                              ifelse(dmr$meth.diff < -25 & dmr$qvalue < 0.05, "Hypomethylated", "Not Significant"))
    
    p3 <- ggplot(dmr, aes(x = meth.diff, y = log10q, color = significant)) +
        geom_point(alpha = 0.6, size = 1) +
        theme_minimal() +
        scale_color_manual(values = c("Hypermethylated" = "red", 
                                      "Hypomethylated" = "blue", 
                                      "Not Significant" = "grey")) +
        labs(title = "DMR Volcano Plot", 
             x = "Methylation Difference (%)", 
             y = "-log10(q-value)") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
        geom_vline(xintercept = c(-25, 25), linetype = "dashed", color = "black")
    
    ggsave("results/05_reports/dmr_volcano.pdf", p3, width = 8, height = 6)
}

# 4. 生成组合报告
pdf("results/05_reports/methylation_summary.pdf", width = 10, height = 8)
print(p1)
if (exists("p2")) print(p2)
if (exists("p3")) print(p3)
dev.off()

cat("Visualization completed. Summary plots saved to results/05_reports/\n")