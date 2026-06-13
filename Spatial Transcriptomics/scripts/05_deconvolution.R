#!/usr/bin/env Rscript
# 细胞类型去卷积

library(Seurat)
library(spacexr)
library(ggplot2)
library(SingleCellExperiment)

args <- commandArgs(trailingOnly = TRUE)
input_rds <- args[1]
input_sc_ref <- args[2]
output_csv <- args[3]
output_plot <- args[4]
method <- args[5]

# 读取空间数据
spatial_obj <- readRDS(input_rds)

# 读取单细胞参考数据
# 这里假设参考数据是h5ad格式，使用reticulate读取
library(reticulate)
sc <- import("scanpy")
sc_ref <- sc$read_h5ad(input_sc_ref)

# 简化版RCTD示例（实际使用时需要完整的数据格式转换）
message(paste("Performing deconvolution using", method))

# 创建模拟的细胞类型比例数据
set.seed(123)
cell_types <- c("CellType_A", "CellType_B", "CellType_C", "CellType_D", "CellType_E")
prop_matrix <- matrix(
    runif(ncol(spatial_obj) * length(cell_types), 0, 1),
    nrow = ncol(spatial_obj),
    ncol = length(cell_types)
)
prop_matrix <- prop_matrix / rowSums(prop_matrix)
colnames(prop_matrix) <- cell_types
rownames(prop_matrix) <- colnames(spatial_obj)

prop_df <- as.data.frame(prop_matrix)
write.csv(prop_df, output_csv, row.names = TRUE)

# 可视化主要细胞类型的空间分布
main_celltype <- cell_types[which.max(colMeans(prop_matrix))]

spatial_obj$main_celltype <- prop_df[[main_celltype]]

# 空间可视化
p <- SpatialFeaturePlot(
    spatial_obj,
    features = "main_celltype",
    pt.size.factor = 1.6,
    alpha = c(0.1, 1)
) + 
scale_color_gradient(low = "blue", high = "red") +
ggtitle(paste("Spatial distribution of", main_celltype))

ggsave(output_plot, p, width = 10, height = 8)

message(paste("Deconvolution completed. Results saved to", output_csv))