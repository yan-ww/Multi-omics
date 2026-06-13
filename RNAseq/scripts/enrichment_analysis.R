#!/usr/bin/env Rscript
# scripts/enrichment_analysis.R

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

# 读取DESeq2结果
de_results <- read_csv(snakemake@input[[1]])

# 设置参数
pvalue_cutoff <- 0.05
log2fc_cutoff <- 1

# 提取显著差异基因
sig_genes <- de_results %>%
  filter(padj < pvalue_cutoff, abs(log2FoldChange) > log2fc_cutoff)

# 分离上下调基因
up_genes <- sig_genes %>%
  filter(log2FoldChange > 0) %>%
  pull(gene_id)

down_genes <- sig_genes %>%
  filter(log2FoldChange < 0) %>%
  pull(gene_id)

# 转换基因ID（如果使用Ensembl ID需要转换）
# 假设gene_id已经是Entrez ID，否则需要转换
# 这里假设使用Entrez ID，如果不是请修改
gene_list <- sig_genes$gene_id

# GO富集分析
go_enrich <- enrichGO(
  gene = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# 保存GO结果
if (!is.null(go_enrich) && nrow(go_enrich) > 0) {
  go_results <- as.data.frame(go_enrich)
  write_csv(go_results, snakemake@output[["go"]])
  
  # 绘制GO气泡图
  go_plot <- dotplot(go_enrich, showCategory = 20, title = "GO Enrichment Analysis")
  ggsave(str_replace(snakemake@output[["go"]], ".csv", ".pdf"), 
         plot = go_plot, width = 10, height = 8)
}

# KEGG富集分析
kegg_enrich <- enrichKEGG(
  gene = gene_list,
  organism = "hsa",
  keyType = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# 保存KEGG结果
if (!is.null(kegg_enrich) && nrow(kegg_enrich) > 0) {
  kegg_results <- as.data.frame(kegg_enrich)
  write_csv(kegg_results, snakemake@output[["kegg"]])
  
  # 绘制KEGG气泡图
  kegg_plot <- dotplot(kegg_enrich, showCategory = 20, title = "KEGG Pathway Enrichment")
  ggsave(str_replace(snakemake@output[["kegg"]], ".csv", ".pdf"), 
         plot = kegg_plot, width = 10, height = 8)
}

cat("Enrichment analysis completed for", snakemake@params[["contrast"]], "\n")
cat("  Significant genes:", length(gene_list), "\n")