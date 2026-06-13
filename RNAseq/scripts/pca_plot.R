#!/usr/bin/env Rscript
# scripts/pca_plot.R

library(ggplot2)
library(tidyverse)

# 读取归一化counts和样本信息
norm_counts <- read_csv("results/deseq2/all_samples_normalized_counts.csv")
samples_df <- read_tsv("config/samples.tsv")

# 准备PCA数据
# 使用variance-stabilizing transformation或rlog变换
# 这里简单使用log2变换（实际DESeq2对象中可以做更好）
expression_matrix <- norm_counts %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

# 对counts进行log2变换（添加伪计数）
log_expression <- log2(expression_matrix + 1)

# 计算PCA
pca_result <- prcomp(t(log_expression), scale = TRUE)
pca_df <- as.data.frame(pca_result$x[, 1:4])
pca_df$sample <- rownames(pca_df)

# 合并样本信息
pca_df <- pca_df %>%
  left_join(samples_df, by = "sample")

# 计算方差解释比例
var_explained <- summary(pca_result)$importance[2, 1:4] * 100

# 绘制PCA图
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, shape = batch)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = condition), level = 0.95, linetype = "dashed") +
  labs(
    title = "PCA Plot of Samples",
    x = paste0("PC1: ", round(var_explained[1], 1), "% variance"),
    y = paste0("PC2: ", round(var_explained[2], 1), "% variance"),
    color = "Condition",
    shape = "Batch"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

# 保存图像
ggsave(snakemake@output[[1]], 
       plot = pca_plot,
       width = 8, 
       height = 6, 
       dpi = 300)

cat("PCA plot saved to", snakemake@output[[1]], "\n")