args = commandArgs(T)
setwd(args[1])
suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)
datExpr0<-read.table(args[2],sep = "\t",header=T,row.names=1)
dim(datExpr0)
head(datExpr0)
gsg<-goodSamplesGenes(datExpr0,verbose=3)##检查样本有无缺失值
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
datExpr<-t(datExpr0)
sampleTree<-hclust(dist(datExpr),method = "average")
par(cex=1.0)##设置文字大小
plot(sampleTree,sub = "",xlab = "",cex.lab=1.5,cex.axis=1.5,cex.main=2)
clust<-cutreeStatic(sampleTree,cutHeight = 15,minSize = 10)
nSamples<-nrow(datExpr)
nGenes<-ncol(datExpr)
traitData<-read.csv(args[3])
Samples<-rownames(datExpr)
traitRows<-match(Samples,traitData$SAMPLE)
datTraits<-traitData[traitRows,-1]
rownames(datTraits)<-Samples
collectGarbage()
traitColors<-numbers2colors(datTraits,signed = FALSE)
#png(file="Sample dendrograms and trait heatmap.png",width=400,height=350,res=72)
pdf(file="Sample dendrograms and trait heatmap.pdf",width=5,height=5)
plotDendroAndColors(sampleTree,traitColors,groupLabels = names(datTraits),
                    main="Sample dendrograms and trait heatmap")


powers<-c(c(1:10),seq(from=12,to=30,by=2))
sft<-pickSoftThreshold(datExpr,powerVector = powers,verbose = 1)
sft$powerEstimate

##生成阈值和网络的特征之间的关系
#png(file="Scale_indenpence.png",width=400,height=350,res=72)
pdf(file="Scale_indenpence.pdf",width=5,height=5)
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold(power)",ylab = "Scale Free Topology",type = "n",
     main = paste("Scale independence"))+text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers,cex=1.0,col = "red")+abline(h=0.85,col="red")


##对连接度的均值进行可视化
#png(file="Mean connectivity.png",width=400,height=350,res=72)
pdf(file="Mean connectivity.pdf",width=5,height=5)
plot(sft$fitIndices[,1],sft$fitIndices[,5],
     xlab = "Soft Threshold (power",ylab = "Mean Connectivity",type = "n",
     main = paste("Mean connectivity"))+text(sft$fitIndices[,1],sft$fitIndices[,5],labels = powers,cex = 0.9,col = "red")

##构建基因网络
net<-blockwiseModules(datExpr,power = sft$powerEstimate,TOMType = "unsigned",
                      minModuleSize = 30,reassignThreshold = 0,mergeCutHeight = 0.25,
                      numericLabels = TRUE,pamRespectsDendro = FALSE,saveTOMs = T,verbose = 3)
save(net,file = "WGCNAnet.RData")
table(net$colors)
mergedColors<-labels2colors(net$colors)
#png(file="dendrograms.png",width=400,height=350,res=72)
pdf(file="dendrograms.pdf",width=5,height=5)
plotDendroAndColors(net$dendrograms[[1]],mergedColors[net$blockGenes[[1]]],"Module colors",
                    dendroLabels = FALSE,hang=0.3,addGuide = TRUE,guideHang = 0.5)


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)#对给定的（特征）向量进行重新排序，以使相似的向量（通过相关性度量）彼此相邻
moduleTraitCor = cor(MEs,datTraits, use = "p")#计算module的ME值与表型的相关系数
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
write.csv(moduleTraitCor,"moduleTraitCor.csv")
write.csv(moduleTraitPvalue,"moduleTraitPvalue.csv")
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
#png(file="moduleTraitRelationshipWithPvalue.png",width=800,height=700,res=72)
pdf(file="moduleTraitRelationshipWithPvalue.pdf",width=12,height=8)
par(mar = c(6, 10, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(datTraits),yLabels = names(MEs),ySymbols = names(MEs),
               colorLabels = FALSE,colors = blueWhiteRed(50),textMatrix = textMatrix,setStdMargins = FALSE,
               cex.text = 1.0,zlim = c(-1,1),main = paste("Module-trait relationships"))
dev.off()

g<-substring(colnames(datExpr0),nchar(colnames(datExpr0)))
datF<-datExpr0[,which(g=="F")]
datN<-datExpr0[,which(g=="N")]
multiExpr<-list(A1=list(data=t(datF)),A2=list(data=t(datN)))
multiColor<-list(A1=moduleColors)
mp<-modulePreservation(multiExpr,multiColor,referenceNetworks = 1,verbose = 3,
                       networkType="signed",nPermutations = 30,maxGoldModuleSize = 1000,
                       maxModuleSize = 4000)
stats<-mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2
stats[order(-stats[,2]),c(1:2)]
stats$module = rownames(stats)
library(ggplot2)
library(ggrepel)
ggplot(data=stats,aes(x=moduleSize,y=Zsummary.pres,col=module))+geom_point(alpha=0.8, size=5) +
  theme_bw(base_size=15)+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("ModuleSize") + ylab("Zsummary.pres") +
  ggtitle( "Preservation Zsummary" ) +
  theme(text=element_text(size=20),plot.title = element_text(size=20,hjust = 0.5))+
  scale_colour_manual(values = c(stats$module))+
  ## 去掉图注
  theme(legend.position='none')+
  ## 添加阈值线
  geom_hline(yintercept = c(2,10),lty=4,lwd=1,col=c("blue","red"))+
  ## 添加文本信息
  geom_text_repel(aes(label=module),color="black",alpha = 0.8)
ggsave(filename="Zsummary.pdf")
plotData<-mp$preservation$observed$ref.A1$inColumnsAlsoPresentIn.A2
write.table(plotData,"PreservationRank.txt",sep="\t")
plotData = plotData[which(rownames(plotData) != "grey"),]
ggplot(data=plotData,aes(x=moduleSize,y=medianRank.pres,col=rownames(plotData)))+geom_point(alpha=0.8, size=5) +
  theme_bw(base_size=15)+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  xlab("ModuleSize") + ylab("Preservation Median rank") +
  ggtitle( "Preservation Median rank" ) +
  theme(text=element_text(size=20),plot.title = element_text(size=20,hjust = 0.5))+
  scale_colour_manual(values = c(rownames(plotData)))+
  theme(legend.position='none')+
  geom_text_repel(aes(label=rownames(plotData)),color="black",alpha = 0.8)
ggsave(filename="preservationRank.pdf")
dev.off()
rm(list=ls())


