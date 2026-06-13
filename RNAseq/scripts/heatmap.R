#!/usr/bin/env Rscript
# scripts/heatmap.R

library(pheatmap)
library(tidyverse)

# 读取数据
norm_counts <- read_csv(snakemake@input[["norm_counts"]])
de_results <- read_csv(snakemake@input[["de_results"]])

# 获取显著差异基因
sig_genes <- de_results %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(50)  # Top 50显著基因

# 提取这些基因的表达量
expression_matrix <- norm_counts %>%
  filter(gene_id %in% sig_genes$gene_id) %>%
  column_to_rownames("gene_id") %>%
  select(-gene_name)

# 标准化表达量（z-score）
expression_matrix_scaled <- t(scale(t(expression_matrix)))

# 准备样本注释
sample_annotation <- data.frame(
  row.names = colnames(expression_matrix),
  Condition = sub("_.*", "", colnames(expression_matrix))
)

# 颜色设置
annotation_colors <- list(
  Condition = c("control" = "blue", "treatment" = "red")
)

# 绘制热图
heatmap_plot <- pheatmap(
  expression_matrix_scaled,
  main = paste("Top 50 DEGs Heatmap -", snakemake@params[["contrast"]]),
  clustering_method = "complete",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  annotation_col = sample_annotation,
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  border_color = NA,
  filename = snakemake@output[[1]],
  width = 10,
  height = 12
)

cat("Heatmap saved to", snakemake@output[[1]], "\n")