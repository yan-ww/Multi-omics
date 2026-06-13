#!/usr/bin/env Rscript
# scripts/volcano_plot.R

library(ggplot2)
library(ggrepel)
library(tidyverse)

# 读取DESeq2结果
de_results <- read_csv(snakemake@input[[1]])

# 设置阈值
padj_cutoff <- 0.05
log2fc_cutoff <- 1

# 添加显著性标记
de_results <- de_results %>%
  mutate(
    significance = case_when(
      is.na(padj) ~ "NS",
      padj < padj_cutoff & abs(log2FoldChange) > log2fc_cutoff ~ "Significant",
      padj < padj_cutoff & abs(log2FoldChange) <= log2fc_cutoff ~ "Significant (low FC)",
      TRUE ~ "Not significant"
    ),
    # 添加-log10(padj)用于绘图
    neg_log10_padj = -log10(padj)
  )

# 选择top基因用于标记（按padj排序）
top_genes <- de_results %>%
  filter(significance == "Significant") %>%
  arrange(padj) %>%
  head(20)

# 创建火山图
volcano_plot <- ggplot(de_results, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Significant" = "red",
    "Significant (low FC)" = "orange",
    "Not significant" = "gray",
    "NS" = "gray"
  )) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "blue") +
  geom_text_repel(data = top_genes,
                  aes(label = gene_name),
                  max.overlaps = 15,
                  size = 3,
                  box.padding = 0.5) +
  labs(
    title = paste("Volcano Plot -", snakemake@params[["contrast"]]),
    x = expression(log[2] ~ "Fold Change"),
    y = expression(-log[10] ~ "Adjusted P-value"),
    color = "Significance"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  xlim(c(-6, 6))  # 限制x轴范围

# 保存图像
ggsave(snakemake@output[[1]], 
       plot = volcano_plot,
       width = 10, 
       height = 8, 
       dpi = 300)

cat("Volcano plot saved to", snakemake@output[[1]], "\n")