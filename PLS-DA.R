setwd("D:/")
library(mixOmics)
library("DESeq2")
count<-read.csv("count.csv",header = T,row.names = 1)
count[1:8,1:8]
group<-read.csv("Group.CSV",header = T)
group[1:8,1:8]
c<-colnames(count)
condition<-substring(c,nchar(c))
coldata<-data.frame(row.names=colnames(count),condition)
dds<-DESeqDataSetFromMatrix(count,colData = coldata,design=~condition)
dds<-DESeq(dds)
vsd<-varianceStabilizingTransformation(dds,blind = FALSE)
vstmat<-assay(vsd)
pls_ana <- plsda(t(vstmat),condition,ncomp = 2)
plotIndiv(pls_ana,comp = c(1,2),group = condition,ind.names = F, 
          ellipse = T, legend = TRUE,style = 'ggplot2',pch =16,
          cex =5,title="FCD3",guide="none")
ggsave("temporalFCD1_scaled.png")
dev.off()
rm(list=ls())
