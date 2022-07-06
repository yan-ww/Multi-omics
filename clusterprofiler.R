##clusterprofiler
setwd("E:/BIO/EP/data/new/DEbybrainpart/new10.28/")
library(clusterProfiler)
library(org.Hs.eg.db)

gene<-volcano$gene[volcano$change!="NOT"]
library(clusterProfiler)
library(org.Hs.eg.db)
entrezid<-bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
enrich.BP<-enrichGO(gene = entrezid$ENTREZID,OrgDb = "org.Hs.eg.db",
                    ont = "BP" ,pAdjustMethod = "fdr",pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2,readable = T,keyType = "ENTREZID")
enrich.MF<-enrichGO(gene = entrezid$ENTREZID,OrgDb = "org.Hs.eg.db",
                    ont = "MF" ,pAdjustMethod = "fdr",pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2,readable = T,keyType = "ENTREZID")
ego <- enrichKEGG(gene = entd$ENTREZID,keyType = "kegg",organism  = 'hsa',pvalueCutoff  = 0.05,pAdjustMethod  = "BH",qvalueCutoff  = 0.05)
##KEGG

dotplot(enrich.BP)
write.csv(enrich.BP,"encichBP.csv")




