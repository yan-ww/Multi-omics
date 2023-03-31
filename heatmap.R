setwd("E:/BIO/EP/data/")
library(pheatmap)


#在做热图前最重要的一步是标准化
#方法一：DEseq标准化方法
library(DESeq2)
DEgenecount<-read.csv("DE")
rownames(DEgenecount)<-DEgenecount[,1]
DEgenecount<-DEgenecount[,-1]
group<-read.csv("group.CSV")
condition<-group$GROUP
coldata<-data.frame(row.names=colnames(DEgenecount),condition)
dds <- DESeqDataSetFromMatrix(DEgenecount, coldata, design = ~condition)
dds <- DESeq(dds) 
vsd<-varianceStabilizingTransformation(dds,blind = FALSE)
vstmat<-assay(vsd)
vstp<-t(scale(t(vstmat)))
table(abs(vstp)>2)
vstp[vstp>=2]=2
vstp[vstp<=-2]=-2
pheatmap(vstp,show_rownames = F,annotation_col = coldata)



#方法二：vegen包的标准化函数
library(vegan)
d<-decostand(DE,"standardize",MARGIN = 1)#画图
pheatmap(vstmat,show_rownames = FALSE,show_colnames = FALSE)

dev.off()
rm(list = ls())
