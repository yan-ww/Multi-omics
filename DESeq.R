setwd("E:/BIO/EP/analysis/")
library(DESeq2)
library(ggplot2)

count<-read.csv("count-FvN.csv")
rownames(count)=count[,1]
count<- count[,-1]
group<-read.csv("group-all.CSV")
condition<-group$GROUP
people<-group$PEOPLE
coldata<-data.frame(row.names=colnames(count),condition,people)  
dds <- DESeqDataSetFromMatrix(DEgenecount, coldata, design = ~people+condition)
dds <- DESeq(dds) 
res = results(dds, contrast=c("condition", "F", "N"))  
res = res[order(res$pvalue),]
summary(res)  
res$change<-as.factor(ifelse((res$pvalue) < 0.05 & (res$log2FoldChange) >0.6,"UP",
                             ifelse(res$pvalue < 0.05 & res$log2FoldChange<(-0.6),"DOWN","NOT")))
write.csv(res,"result.csv")
table(res$change)
volcano<-res[c("pvalue","log2FoldChange","change")]
volcano$gene<-rownames(res)
volcano<-as.data.frame(volcano)
ggplot(data=volcano, aes(x=log2FoldChange, y=-log10 (pvalue),color=change)) + 
  geom_point(alpha=0.5, size=1.75) +
  theme_set(theme_bw(base_size=20))+ 
  scale_color_manual(values=c('blue','black','red'))+
  xlab("log[2] fold change(F/N)") +
  ##也可以写成(x=expression(log[2](F/N)), y=expression( -log[10](P.Adjusted.Value))) 
  ylab("-log[10] pvalue") + 
  theme(plot.title = element_text(size=15,hjust = 0.5),legend.position="right", 
        legend.title = element_blank()) + 
  geom_hline(yintercept = 1.30,linetype="dashed")+
  geom_vline(xintercept = c(-0.6,0.6),linetype="dashed")+
  geom_text(data = volcano[volcano$pvalue<0.05&abs(volcano$log2FoldChange)>1,],
            aes(label = gene),
            size = 3,color = "black",show.legend = FALSE )
ggsave("volcano-FvN.png")
dev.off()
rm(list=ls())
