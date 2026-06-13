#!/usr/bin/env Rscript
# 数据预处理和质量控制

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(DropletUtils)

# 读取参数
args <- commandArgs(trailingOnly = TRUE)
input_h5 <- args[1]
input_tissue <- args[2]
input_scalefactors <- args[3]
output_rds <- args[4]
output_qc_plot <- args[5]
min_genes <- as.numeric(args[6])
max_mito <- as.numeric(args[7])

# 读取10X空间数据
message("Loading 10X Visium data...")
data <- Load10X_Spatial(
    data.dir = dirname(input_h5),
    filename = basename(input_h5),
    assay = "Spatial"
)

# 添加图像信息
image <- Read10X_Image(
    image.dir = dirname(input_tissue),
    image.name = "tissue_hires_image.png"
)
data <- AddMetaData(data, image@coordinates)

# 质量控制
message("Performing QC...")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

# 筛选spot
data <- subset(
    data,
    subset = nFeature_Spatial > min_genes &
             percent.mt < max_mito
)

# QC可视化
qc_plot <- VlnPlot(
    data,
    features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"),
    ncol = 3,
    pt.size = 0.1
) + theme(legend.position = "none")

ggsave(output_qc_plot, qc_plot, width = 12, height = 4)

# 保存
saveRDS(data, output_rds)
message(paste("Saved", output_rds))
message(paste("Remaining spots:", ncol(data)))