#!/usr/bin/env Rscript
# 空间聚类

library(Seurat)
library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
input_rds <- args[1]
output_rds <- args[2]
output_plot <- args[3]
resolution <- as.numeric(args[4])

# 读取数据
data <- readRDS(input_rds)

# 聚类
data <- FindNeighbors(data, dims = 1:30)
data <- FindClusters(data, resolution = resolution)

# 可视化
p1 <- DimPlot(data, reduction = "umap", label = TRUE) + 
      ggtitle("UMAP Clustering")
p2 <- SpatialDimPlot(data, label = TRUE, label.size = 3) + 
      ggtitle("Spatial Clusters")

combined_plot <- p1 + p2
ggsave(output_plot, combined_plot, width = 14, height = 6)

# 保存
saveRDS(data, output_rds)
message(paste("Number of clusters:", length(unique(data$seurat_clusters))))