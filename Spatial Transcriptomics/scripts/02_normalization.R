#!/usr/bin/env Rscript
# 标准化和降维

library(Seurat)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
input_rds <- args[1]
output_rds <- args[2]
output_pca_plot <- args[3]
method <- args[4]

# 读取数据
data <- readRDS(input_rds)

# 标准化
message(paste("Normalization using", method))
if (method == "SCT") {
    data <- SCTransform(data, assay = "Spatial", verbose = FALSE)
} else {
    data <- NormalizeData(data)
    data <- FindVariableFeatures(data)
    data <- ScaleData(data)
}

# PCA降维
data <- RunPCA(data, assay = "SCT", npcs = 50)

# PCA可视化
pca_plot <- ElbowPlot(data, ndims = 50)
ggsave(output_pca_plot, pca_plot, width = 8, height = 6)

# 保存
saveRDS(data, output_rds)
message("Normalization completed")