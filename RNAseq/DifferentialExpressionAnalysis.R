args = commandArgs(T)
options(stringsAsFactors = FALSE)
setwd(args[1])
count = read.table(args[2],sep="\t",row.names=1,header=T)

suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
options(download.file.method="wget")
DE<-function(a){
    b<-colnames(a)
    condition<-substring(b,nchar(b))
    people<-substr(b,3,nchar(b)-1)
    coldata<-data.frame(row.names=colnames(a),people,condition)
    dds<-DESeqDataSetFromMatrix(a,coldata,design=~people+condition)
    dds<-DESeq(dds)
    res<-results(dds,contrast=c("condition","F","N"))
    res<-res[order(res$pvalue),]
    res$change<-as.factor(ifelse((res$pvalue) < 0.05 & (res$log2FoldChange) > log2(1.5),"up",ifelse(res$pvalue < 0.05 & res$log2FoldChange < (-log2(1.5)),"down","not")))
    return(res)
}


volcano<-function(de.res){
    c<-de.res[c("pvalue","log2FoldChange","change")]
    c$gene<-rownames(de.res)
    c$label<-""
    upgenes<-head(c$gene[which(c$change=="up")],10)
    downgenes<-head(c$gene[which(c$change=="down")],10)
    top10genes<-c(as.character(upgenes),as.character(downgenes))
    c$label[match(top10genes,c$gene)]<-top10genes
    c<-as.data.frame(c)
    volcanoplot<-ggplot(data=c, aes(x=log2FoldChange, y=-log10 (pvalue),color=change)) + geom_point(alpha=0.5, size=1.75) + 
        theme_set(theme_bw(base_size=20))+  scale_color_manual(values=c('#2f5688','#BBBBBB','#CC0000'))+ 
        xlab("log[2] fold change(F/N)") + ylab("-log[10] pvalue") + theme(text=element_text(size=20),plot.title = element_text(size=20,hjust = 0.5),
        legend.position="right",legend.title = element_blank()) + geom_hline(yintercept = 1.30,linetype="dashed")+ 
        geom_vline(xintercept = c(-0.6,0.6),linetype="dashed")+ geom_text(aes(label = label),size = 5,color = "black",show.legend = FALSE )
        #改动如下：y的值从pvalue改成padj，geom_vline的xintercept从(-0.6,0.6)改成c(-1,1)
    return(volcanoplot)
}

double_enrich_plot = function(df,cate){
  df$pl = ifelse(df$category == "up",-log10(df$p.adjust),log10(df$p.adjust))
  df = arrange(df,category,pl)
  df$Description = factor(df$Description,levels = unique(df$Description),ordered = TRUE)
  tmp = with(df, labeling::extended(range(pl)[1], range(pl)[2], m = 10))
  lm = tmp[c(1,length(tmp))]
  lm = c(floor(min(df$pl)),ceiling(max(df$pl)))
  ggplot(df, aes(x=Description, y= pl)) +labs(y="p.adjust")+
    geom_bar(stat='identity', aes(fill=category), width=.7)+
    scale_fill_manual(values = c("#2874C5", "#f87669"))+
    coord_flip()+theme_light() +ylim(lm)+
    scale_x_discrete(labels=function(x) str_wrap(x, width=30))+
    scale_y_continuous(breaks = tmp,labels = abs(tmp))+
    theme(text=element_text(size=20),panel.border = element_blank())
  ggsave(filename = sprintf("%sofDEG.pdf",cate))
  ggsave(filename = sprintf("%sofDEG.png",cate),width=500,height=400,units="mm")
}

utils::globalVariables(c("category","pl","Description"))

########## DifferentialExpressionAnalysis ##########
# 三种方法：limma,DESeq2,EdgeR
#limma
#limma算法原理：核心在于lmfit()和eBays()，参考https://www.jianshu.com/p/716cb81d1f61
library(limma)
count<-read.csv("count-FvN.csv",row.names=1)
# count的行名是基因名，列名是样本名
condition = c(rep("F",3),rep("N",3)) %>% factor(.,levels = c("F","N"),ordered = F) 
design = model.matrix(~0+condition)
colnames(design) = c("F","N")
df.fit = lmFit(count,design) #数据与list进行匹配
df.matrix = makeContrasts(F-N,levels=design)
fit = contrasts.fit(df.fit,df.matrix)
fit = eBayes(fit)
res = topTable(fit,number=Inf)
res = res[order(res$P.Value),]

#DEseq2 methods
DEresult = DE(count)
DE.res = data.frame(DEresult)
rownames(DE.res) = trimws(rownames(DE.res),whitespace = " ",which = "l")
write.table(DE.res,file="DEGresult.txt",sep="\t")
volcanoplot = volcano(DE.res)
ggsave(volcanoplot,filename = "Volcanoplot.pdf")
entrezid = bitr(rownames(DE.res),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db") 
DE.res$SYMBOL = rownames(DE.res)
data = merge(entrezid,DE.res,by="SYMBOL",all = F)
data = data[which(data$change != "not"),]

#GO分析
UPenriGo<-enrichGO(gene = data$ENTREZID[data$change=="up"],OrgDb = "org.Hs.eg.db",ont = "ALL" ,pAdjustMethod = "fdr",pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2,readable = T,keyType = "ENTREZID")
DOWNenriGo<-enrichGO(gene = data$ENTREZID[data$change=="down"],OrgDb = "org.Hs.eg.db",ont = "ALL" ,pAdjustMethod = "fdr",pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2,readable = T,keyType = "ENTREZID")
UPgo.res = data.frame(UPenriGo)
UPgo.res$category = rep("up",nrow(UPgo.res))
DOWNgo.res = data.frame(DOWNenriGo)
DOWNgo.res$category = rep("down",nrow(DOWNgo.res))
DEGgo = rbind(UPgo.res,DOWNgo.res)
write.csv(DEGgo,file ="encichGOofDEG.csv" ,quote = FALSE)
if(nrow(UPgo.res)>10){
  UPgo.res = UPgo.res[1:10,]
}
if(nrow(DOWNgo.res)>10){
  DOWNgo.res = DOWNgo.res[1:10,]
}
DEGgo = rbind(UPgo.res,DOWNgo.res)
double_enrich_plot(DEGgo,"GO")

#KEGG分析
UPenriKegg<- enrichKEGG(gene = data$ENTREZID[data$change=="up"],keyType = "kegg",organism  = "hsa",pvalueCutoff  = 0.05,pAdjustMethod  = "BH",qvalueCutoff  = 0.2)
UPkk = setReadable(UPenriKegg,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
DOWNenriKegg<- enrichKEGG(gene = data$ENTREZID[data$change=="down"],keyType = "kegg",organism  = "hsa",pvalueCutoff  = 0.05,pAdjustMethod  = "BH",qvalueCutoff  = 0.2)
DOWNkk = setReadable(DOWNenriKegg,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
UPkegg = data.frame(UPkk)
UPkegg$category = rep("up",nrow(UPkegg))
DOWNkegg = data.frame(DOWNkk)
DOWNkegg$category = rep("down",nrow(DOWNkegg))
DEGkegg = rbind(UPkegg,DOWNkegg)
write.csv(DEGkegg,file ="encichKEGGofDEG.csv" ,quote = FALSE)
if(nrow(UPkegg)>10){
  UPkegg = UPkegg[1:10,]
}
if(nrow(DOWNkegg)>10){
  DOWNkegg = DOWNkegg[1:10,]
}
DEGkegg = rbind(UPkegg,DOWNkegg)
double_enrich_plot(DEGkegg,"KEGG")

