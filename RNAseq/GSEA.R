# get args----
args = commandArgs(T)
options(stringsAsFactors = FALSE)
setwd(args[1])
DEGpath=args[2]


# library----
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(enrichplot))
options(download.file.method="wget")


# GO
gseadata = read.table(DEGpath,sep="\t",header=T,row.names=1)
gseadata = gseadata[order(gseadata$log2FoldChange,decreasing = T),]
gene_list = gseadata$log2FoldChange
names(gene_list) = rownames(gseadata)
gene_list = sort(gene_list,decreasing = T)
gseago.res = gseGO(gene_list,ont = "BP",OrgDb = org.Hs.eg.db ,pAdjustMethod = "BH",pvalueCutoff = 0.05,keyType = "SYMBOL")
png(filename = "GSEAGO.png",width = 800, height = 1000, units = "px", pointsize = 12,bg="white")
dotplot(gseago.res,showCategory=10,split=".sign",font.size=20) + facet_grid(.~.sign)
dev.off()
pdf(file = "GSEAGO.pdf",width = 8, height = 10)
dotplot(gseago.res,showCategory=10,split=".sign",font.size=20) + facet_grid(.~.sign)
dev.off()
gseago.res = data.frame(gseago.res)
write.csv(gseago.res,"GSEAGO.CSV",quote =FALSE)


#KEGG
gseadata = read.table(DEGpath,sep="\t",header=T,row.names=1)
gene <- bitr(rownames(gseadata), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
data <- merge(gene,gseadata,by.x = "SYMBOL",by.y = "row.names")
data = data[order(data$log2FoldChange,decreasing = T),]
gene_list <- data$log2FoldChange
names(gene_list) <- data$ENTREZID
gene_list <- sort(gene_list,decreasing = T)
res <- gseKEGG(gene_list, organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")
png(filename = "GSEAKEGG.png",width = 800, height = 1000, units = "px", pointsize = 12,bg="white")
dotplot(res,showCategory=10,split=".sign",font.size=20) + facet_grid(.~.sign)
dev.off()
pdf(file = "GSEAKEGG.pdf",width = 8, height = 10)
dotplot(res,showCategory=10,split=".sign",font.size=20) + facet_grid(.~.sign)
dev.off()
kk = setReadable(res,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
write.csv(kk,"GSEKEGG.CSV",quote =FALSE)
png(filename ="NeuroGSEA.png",width=16,height=12,unit="in",res=300, pointsize = 12,bg="white")
gseaplot2(res, title = res$Description[1], geneSetID = 1)
dev.off()



