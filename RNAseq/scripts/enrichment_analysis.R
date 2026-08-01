#!/usr/bin/env Rscript

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

de_results <- read_csv(snakemake@input[[1]])
sig_genes <- de_results %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

gene_ids <- unique(sub("\\..*$", "", sig_genes$gene_id))
gene_map <- if (length(gene_ids)) {
  bitr(gene_ids, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
} else {
  tibble(ENTREZID = character())
}
gene_list <- unique(gene_map$ENTREZID)

write_empty <- function(path) write_csv(tibble(), path)

if (length(gene_list) == 0) {
  write_empty(snakemake@output[["go"]])
  write_empty(snakemake@output[["kegg"]])
} else {
  go_enrich <- enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )

  if (!is.null(go_enrich) && nrow(as.data.frame(go_enrich)) > 0) {
    write_csv(as.data.frame(go_enrich), snakemake@output[["go"]])
    ggsave(sub("\\.csv$", ".pdf", snakemake@output[["go"]]),
           dotplot(go_enrich, showCategory = 20, title = "GO Enrichment Analysis"),
           width = 10, height = 8)
  } else {
    write_empty(snakemake@output[["go"]])
  }

  kegg_enrich <- enrichKEGG(
    gene = gene_list,
    organism = "hsa",
    keyType = "kegg",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )

  if (!is.null(kegg_enrich) && nrow(as.data.frame(kegg_enrich)) > 0) {
    write_csv(as.data.frame(kegg_enrich), snakemake@output[["kegg"]])
    ggsave(sub("\\.csv$", ".pdf", snakemake@output[["kegg"]]),
           dotplot(kegg_enrich, showCategory = 20, title = "KEGG Pathway Enrichment"),
           width = 10, height = 8)
  } else {
    write_empty(snakemake@output[["kegg"]])
  }
}

cat("Enrichment analysis completed for", snakemake@params[["contrast"]], "\n")
cat("  Mappable significant genes:", length(gene_list), "\n")
