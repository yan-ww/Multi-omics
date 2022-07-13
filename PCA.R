setwd("E:/")
library(DESeq2)
library(factoextra)
library(FactoMineR)
library(vegan)

count<-read.csv("DEgenecount.csv",header=T,row.names=1)
datacount<-count
summary(is.numeric(datacount[,1]))
b<-colnames(datacount)
condition<-substring(b,nchar(b))
people<-substr(b,3,nchar(b)-1)
coldata<-data.frame(row.names=colnames(datacount),condition,people)
dds <- DESeqDataSetFromMatrix(datacount, coldata, design = ~people+condition)
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
            title="FvsN")+
theme(plot.title = element_text(hjust = 0.5))



