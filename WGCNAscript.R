##WGCNA
setwd("E:/BIO/WGCNA/")
library(WGCNA)
options(stringsAsFactors = FALSE)
femData<-read.csv("LiverFemale3600.csv")
dim(femData)
head(femData)
datExpr0<-as.data.frame(t(femData[,-c(1:8)]))
names(datExpr0)<-femData$substanceBXH
row.names(datExpr0)<-names(femData)[-c(1:8)]
datExpr0[1:8,1:8]
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
sampleTree<-hclust(dist(datExpr0),method = "average")
par(cex=0.5)##设置文字大小
plot(sampleTree,sub = "",xlab = "",cex.lab=1.5,cex.axis=1.5,cex.main=2)
abline(h=15,col="red")##在图上画线
clust<-cutreeStatic(sampleTree,cutHeight = 15,minSize = 10)
table(clust)
keepSamples<-(clust==1)
datExpr<-datExpr0[keepSamples,]
nSamples<-nrow(datExpr)
nGenes<-ncol(datExpr)
#加载性状数据
traitData<-read.csv("ClinicalTraits.csv")
dim(traitData)
allTraits<-traitData[,-c(31,16)]
allTraits<-allTraits[,c(2,11:36)]
head(allTraits)
femaleSamples<-rownames(datExpr)
traitRows<-match(femaleSamples,allTraits$Mice)
datTraits<-allTraits[traitRows,-1]
rownames(datTraits)<-allTraits[traitRows,1]
collectGarbage()

## Re-cluster samples
sampleTree2<-hclust(dist(datExpr),method = "average")
traitColors<-numbers2colors(datTraits,signed = FALSE)
plotDendroAndColors(sampleTree2,traitColors,groupLabels = names(datTraits),
                    main="Sample dendrogran and trait heatmap")
powers<-c(c(1:10),seq(from=12,to=20,by=2))
powers
sft<-pickSoftThreshold(datExpr,powerVector = powers,verbose = 0)
str(sft)
cex1=0.9

##生成阈值和网络的特征之间的关系
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold(power)",ylab = "Scale Free Topology",type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers,cex = cex1,col = "red")
abline(h=0,9,col="red")

##对连接度的均值进行可视化
plot(sft$fitIndices[,1],sft$fitIndices[,5],
     xlab = "Soft Threshold (power",ylab = "Mean Connectivity",type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels = powers,cex = cex1,col = "red")

##构建基因网络
net<-blockwiseModules(datExpr,power = sft$powerEstimate,TOMType = "unsigned",
                      minModuleSize = 30,reassignThreshold = 0,mergeCutHeight = 0.25,
                      numericLabels = TRUE,pamRespectsDendro = FALSE,saveTOMs = TRUE,
                      saveTOMFileBase = "femaleMouseTOM",verbose = 3)
#如果这里出错了见https://blog.csdn.net/liyunfan00/article/details/91686840
##cor <- WGCNA::cor
##cor<-stats::cor
table(net$colors)
mergedColors<-labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]],mergedColors[net$blockGenes[[1]]],"Module colors",
                    dendroLabels = FALSE,hang=0.03,addGuide = TRUE,guideHang = 0.05)

##表型与模块之间的相关系数
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
##用彩色标签重新计算MEs，在给定的单个数据集中计算模块的模块本征基因
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)#对给定的（特征）向量进行重新排序，以使相似的向量（通过相关性度量）彼此相邻
moduleTraitCor = cor(MEs, traitExpr, use = "p")#计算module的ME值与表型的相关系数
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# sizeGrWindow(10,6)
# 显示相关性及其p值
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 10, 3, 3));
# Display the correlation values within a heatmap plot\
# ySymbols 当ylabels使用时所使用的其他标签； colorLabels 应该使用颜色标签吗
# colors 颜色； textMatrix 单元格名字
labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(traitExpr),yLabels = names(MEs),ySymbols = names(MEs),
               colorLabels = FALSE,colors = blueWhiteRed(50),textMatrix = textMatrix,setStdMargins = FALSE,
               cex.text = 0.4,zlim = c(-1,1),
               main = paste("Module-trait relationships"))

##绘制有表型参与的聚类树
# Recalculate module eigengenes 重新计算基因特征值
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
group = as.data.frame(traitExpr$GROUP);
names(group) = "group"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, group))
# Plot the relationships among the eigengenes and the trait
#sizeGrWindow(5,7.5);
par(cex = 0.9)
# 画树形图 marDendro给出树状图的边距设置，marHeatmap热图边距设置
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2),
                      cex.lab = 0.8, xLabelsAngle= 90)


##量化阵列上所有基因与每个模块的相似性寻找重要模块
# Define variable weight containing the weight column of datTrait定义包含数据特征权重列的变量权重
weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
# 基因和模块的相关系数
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="")

##在棕色模块中绘制了基因重要性与模块成员关系的散点图
module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

##网络分析结果总结，将此统计信息与基因注释合并
names(datExpr)# 提取表带数据样本名称
names(datExpr)[moduleColors=="brown"]# 指定颜色数据名称
# 基因注释数据
annot = read.csv(file = "../data/FemaleLiver-Data/GeneAnnotation.csv");
dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.
# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = probes,
                       geneSymbol = annot$gene_symbol[probes2annot],
                       LocusLinkID = annot$LocusLinkID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for(mod in 1:ncol(geneModuleMembership)){
        oldNames = names(geneInfo0)
        geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], MMPvalue[, modOrder[mod]]);
        names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),paste("p.MM.", modNames[modOrder[mod]], sep=""))
        print(mod)  
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "../result/geneInfo.csv")

##