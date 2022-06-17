setwd("E:/BIO/EP/data/")
library(DESeq2)
library(factoextra)
library(FactoMineR)
library(vegan)

count<-read.csv("DEgenecount.csv",header=T,row.names=1)
summary(is.numeric(count[,1]))
group<-read.csv("group.CSV")
condition<-group$GROUP
coldata<-data.frame(row.names=colnames(count),condition)
dds <- DESeqDataSetFromMatrix(count, coldata, design = ~condition)
dds <- DESeq(dds) 
vsd<-varianceStabilizingTransformation(dds,blind = FALSE)
#这里是用DEseq的标准化方法，也可以直接用fpkm
vstmat<-assay(vsd)
pca<-PCA(t(vstmat),scale.unit = FALSE,ncp = 5,graph = FALSE)
group<-read.csv("group.CSV")
fviz_pca_ind(pca, geom.ind  = c("point","text"),
            col.ind = group$GROUP, # color by groups
            pointsize = 1.5, pointshape = 16,
            addEllipses = T, 
            repel = F ,
            legend.title = "Tissue",
            title="youth vs adults"
)+theme(plot.title = element_text(hjust = 0.5))



