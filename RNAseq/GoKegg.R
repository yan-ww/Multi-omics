#! /public/home/anyu_lab/wengyp/anaconda3/envs/r/bin/Rscript
arg=commandArgs(T)
setwd(arg[1])
suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(DOSE))

GOplot = function(go.res,color) {
  goBP <- subset(go.res, subset = (ONTOLOGY == "BP"))
  if(nrow(goBP)>10){
    goBP = goBP[1:10,]
  }
  goCC <- subset(go.res, subset = (ONTOLOGY == "CC"))
  if(nrow(goCC)>10){
    goCC = goCC[1:10,]
  }
  goMF <- subset(go.res, subset = (ONTOLOGY == "MF"))
  if(nrow(goMF)>10){
    goMF = goMF[1:10,]
  }
  term <- rbind(goBP, goCC, goMF)
  term$Description=factor(term$Description, levels = term$Description)
  gobar = ggplot(data=term, aes(x=Description, y=Count,color=ONTOLOGY)) + geom_bar(stat="identity",width=0.8,aes(fill=ONTOLOGY)) + 
    coord_flip()+ xlab("GO term") + ylab("Num of Genes") + theme_bw() + labs(title = sprintf("GO Enrichment of %s",color)) + 
    theme(text=element_text(size=20))+scale_x_discrete(labels = function(x) str_wrap(x,width=50))
  #ggsave(gobar,filename = sprintf("GOafWGCNA/GOplotOf%s.png",color),width=500,height=400,units="mm")
  ggsave(gobar,filename = sprintf("GOafWGCNA/GOplotOf%s.pdf",color),width=500,height=400,units="mm")
}

KEGGplot = function(kegg,color){
  if(nrow(kegg)>15){
    kegg = kegg[1:15,]
  }
  kegg$Description = factor(kegg$Description,levels = rev(kegg$Description))
  kegg <- mutate(kegg, ratio = parse_ratio(GeneRatio))
  p = ggplot(data = kegg,aes(x = ratio, y = reorder(Description,Count)))+
    geom_point(aes(size = Count,color = -log10(p.adjust)))+
    theme_bw()+
    scale_colour_gradient(low = "green",high = "red")+
    scale_y_discrete(labels = function(x) str_wrap(x,width = 40))+
    labs(x = "GeneRatio",y = "",title = sprintf("KEGG pathway of %s",color),
       color = expression(-log10(p.adjust)),size = "Count")+
    theme(text=element_text(size=20))
  ggsave(p, filename = sprintf("GOafWGCNA/KEGGplotOf%s.pdf",color),width=500,height=400,units="mm")
}


load("WGCNAnet.RData")
options(stringsAsFactors = FALSE)
datExpr0<-read.table(file=arg[2],sep = "\t",header = T,row.names = 1)
#options(clusterProfiler.download.method="wget")
options(download.file.method="wget")
moduleColors = labels2colors(net$colors)
a<-table(moduleColors)
n<-length(a)
for (i in 1:n) {
  color <- names(a[i])
  genes <- rownames(datExpr0)[moduleColors==color]
  entrezid<-bitr(genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db") 
  enriGo<-enrichGO(gene = entrezid$ENTREZID,OrgDb = "org.Hs.eg.db", ont = "ALL" ,pAdjustMethod = "fdr",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T,keyType = "ENTREZID")
  #enriGo<-enrichGO(gene = genes,OrgDb = "org.Hs.eg.db", ont = "ALL" ,pAdjustMethod = "fdr",pvalueCutoff = 0.05,qvalueCutoff = 0.2,readable = T,keyType = "SYMBOL")
  go.res = data.frame(enriGo)
  if(nrow(go.res)>=1){
    GOplot(go.res,color)
  }
  write.csv(enriGo,file = sprintf("GOafWGCNA/enrichGO_%s.csv",color),quote = FALSE)
  if(nrow(go.res)>5){
    edoxGO<-pairwise_termsim(enriGo)
    #png(filename = sprintf("GOafWGCNA/GOtree_%s.png",color),width = 700, height = 480, units = "px", pointsize = 12,bg="white")
    emapplot(edoxGO)
    #ggsave(filename = sprintf("GOafWGCNA/GOcnetplotOf%s.png",color),width=500,height=400,units="mm")
    ggsave(filename = sprintf("GOafWGCNA/GOcnetplotOf%s.pdf",color))
  }
  
  enriKegg<- enrichKEGG(gene = entrezid$ENTREZID,keyType = "kegg",organism  = "hsa",pvalueCutoff  = 0.05,pAdjustMethod  = "BH",qvalueCutoff  = 0.05)
  kegg = data.frame(enriKegg)
  kk = setReadable(enriKegg,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
  if(nrow(kegg)>=1){
    KEGGplot(kegg,color) 
  }
  write.csv(kk,file = sprintf("GOafWGCNA/enrichKEGG_%s.csv",color),quote = FALSE)
  if(nrow(kegg)>5){
    edoxKEGG<-pairwise_termsim(enriKegg)
    #png(filename = sprintf("GOafWGCNA/KEGGtree_%s.png",color),width = 700, height = 480, units = "px", pointsize = 12,bg="white")
    emapplot(edoxKEGG)
    #ggsave(filename = sprintf("GOafWGCNA/KEGGcnetplotOf%s.png",color),width=500,height=400,units="mm")
    ggsave(filename = sprintf("GOafWGCNA/KEGGcnetplotOf%s.pdf",color))
  }
}




