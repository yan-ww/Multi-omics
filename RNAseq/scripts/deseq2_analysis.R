#!/usr/bin/env Rscript
# scripts/deseq2_analysis.R

library(DESeq2)
library(tidyverse)
library(apeglm)

# 加载配置
samples_df <- read_tsv(snakemake@input[["samples_tsv"]])

# 读取counts数据
counts_df <- read_csv(snakemake@input[["counts"]])
counts_matrix <- counts_df %>%
  column_to_rownames("gene_id") %>%
  select(-gene_name) %>%
  as.matrix()

# 准备DESeq2对象
coldata <- samples_df %>%
  select(sample, condition, batch) %>%
  column_to_rownames("sample")

# 确保顺序一致
counts_matrix <- counts_matrix[, rownames(coldata)]

# 创建DESeq2数据集
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = coldata,
  design = ~ batch + condition  # 考虑批次效应
)

# 预过滤（至少10个reads在至少3个样本中）
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# 运行DESeq2
dds <- DESeq(dds, parallel = TRUE)

# 保存DESeq2对象
save(dds, file = snakemake@output[["rdata"]])

# 对每个对比组进行分析
contrasts <- snakemake@params[["contrasts"]]
for (contrast in contrasts) {
  # 解析对比组（例如 "treatment_vs_control"）
  contrast_parts <- strsplit(contrast, "_vs_")[[1]]
  numerator <- contrast_parts[1]
  denominator <- contrast_parts[2]
  
  # 提取结果
  res <- results(dds, 
                 contrast = c("condition", numerator, denominator),
                 alpha = 0.05)
  
  # LFC收缩（使用apeglm）
  coef_name <- paste0("condition_", numerator, "_vs_", denominator)
  res_shrunk <- if (coef_name %in% resultsNames(dds)) {
    lfcShrink(dds, coef = coef_name, type = "apeglm")
  } else {
    res
  }
  
  # 转换为数据框并添加基因名称
  res_df <- as.data.frame(res_shrunk) %>%
    rownames_to_column("gene_id") %>%
    left_join(counts_df %>% select(gene_id, gene_name), by = "gene_id") %>%
    relocate(gene_name, .after = gene_id)
  
  # 输出差异表达结果
  de_output <- str_glue("results/deseq2/{contrast}_DE_results.csv")
  write_csv(res_df, de_output)
  
  # 输出归一化counts
  norm_counts <- counts(dds, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    left_join(counts_df %>% select(gene_id, gene_name), by = "gene_id")
  
  norm_output <- str_glue("results/deseq2/{contrast}_normalized_counts.csv")
  write_csv(norm_counts, norm_output)
  
  # 输出显著基因列表
  sig_genes <- res_df %>%
    filter(padj < 0.05, abs(log2FoldChange) > 1)
  
  write_csv(sig_genes, str_glue("results/deseq2/{contrast}_significant_genes.csv"))
  
  cat("Completed analysis for", contrast, "\n")
  cat("  Total significant genes:", nrow(sig_genes), "\n")
}

# 保存所有样本的归一化counts（用于PCA等）
all_norm_counts <- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("gene_id")
write_csv(all_norm_counts, snakemake@output[["all_norm_counts"]])

cat("DESeq2 analysis completed successfully!\n")
