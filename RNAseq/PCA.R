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


###PCA三维图画法
#在pca分析完得到pca文件之后：
F_color<-rep("red",count(condition=="F"))
names(F_color)<-rownames(coldata)[which(coldata$condition=="F")]
N_color<-rep("blue",count(condition=="N"))
names(N_color)<-rownames(coldata)[which(coldata$condition=="N")]
groups<-c(F_color,N_color)
plot3d(pca$ind$coord[,1:3], # 取前三个主成分
       xlab="Comp.1", ylab="Comp.2", zlab="Comp.3", 
       col=groups, # 按groups填充颜色
       type="s", # 画球，'p' for points, 's' for spheres, 'l' for lines, 'h' for line segments 
       size=1, #球的大小
       lwd=2, box=T)+
  text3d(pca$ind$coord[,1:3],texts = rownames(pca$ind$coord),col="black")
rgl.snapshot("E:/lncRNA_analysis/PCA_new/3dPCAoftemfcd3.png")
