rm(list=ls())
options(stringsAsFactors = FALSE)
library("Matrix")
library("dplyr")
library("Seurat")
library("patchwork")
data<-Seurat::Read10X_h5(filename = "../miniproject2/filtered_feature_bc_matrix.h5",use.names = T)
count<-CreateSeuratObject(counts = data,project = "count",min.cells = 3,min.features = 200)

##计算线粒体比例
count[["percent_mito"]]<-PercentageFeatureSet(object = count,pattern = "^mt-")
head(count@meta.data,5)
feats<-c("nFeature_RNA","nCount_RNA","percent_mito")
VlnPlot(count,features = feats,ncol = 3)
plot1<-FeatureScatter(count,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
plot2<-FeatureScatter(count,feature1 = "nCount_RNA",feature2 = "percent_mito")
plot1+plot2

##过滤数据
count<-subset(count,subset = nFeature_RNA > 200 & nFeature_RNA < 2000)

##数据标准化
hist(colSums(count$RNA@data),breaks = 100,main = "Total expression before normalisation",
     xlab = "Sum of expression")
count<-NormalizeData(count,normalization.method = "LogNormalize",scale.factor = 10000)
hist(colSums(count$RNA@data),breaks = 100,main = "Total expression after normalisation",
     xlab = "Sum of expression")

##变化基因鉴定
count<-FindVariableFeatures(count,selection.method = "vst",nfeatures = 2000)
top10<-head(VariableFeatures(count),10)
plot3<-VariableFeaturePlot(count)
plot4<-LabelPoints(plot = plot3,points = top10,repel = TRUE)
plot4

##降维
count<-ScaleData(count)
count<-RunPCA(count,features = VariableFeatures(object = count))
ElbowPlot(count)
count<-FindNeighbors(count,dims = 1:10)
count<-FindClusters(count,resolution = 0.5)
count<-RunUMAP(count,dims = 1:10)
DimPlot(count,reduction = "umap")


count.markers<-FindAllMarkers(count,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
count.markers%>%group_by(cluster)%>%slice_max(n=2,order_by = avg_log2FC)
FeaturePlot(count,features = c("Tmx4","Lag3","Cd83","Dnajb1","Trbc2","Xcl1","Ly6d","H2-Aa","Lyz2","S100a9"))
VlnPlot(count,features = c("Sparc"))

new.cluster.ids<-c(rep("Microglia",5),"T cells","NK cells","B cells","cDC2","Monocytes","neutophils")
names(new.cluster.ids)<-levels(count)
count<-RenameIdents(count,new.cluster.ids)
DimPlot(count,reduction = "umap",label = TRUE,pt.size = 0.5)+NoLegend()
