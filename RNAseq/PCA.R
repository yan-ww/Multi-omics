setwd(args[1])
library(DESeq2)
library(factoextra)
library(FactoMineR)
library(vegan)

count<-read.csv(args[2],header=T,row.names=1)
coldata=read.table(args[3],headerT)
datacount<-count
summary(is.numeric(datacount[,1]))
dds <- DESeqDataSetFromMatrix(datacount, coldata, design = ~group)
dds <- DESeq(dds) 
vsd<-varianceStabilizingTransformation(dds,blind = FALSE)
#这里是用DEseq的标准化方法，也可以直接用fpkm
vstmat<-assay(vsd)
pca<-PCA(t(vstmat),scale.unit = FALSE,ncp = 5,graph = FALSE)
fviz_pca_ind(pca, geom.ind  = c("point","text"),
            col.ind = condition, # color by groups
            pointsize = 1.5, pointshape = 16,
            addEllipses = T, 
            repel = F ,
            legend.title = "Tissue",
            title="PCA analysis")+
theme(plot.title = element_text(hjust = 0.5))

