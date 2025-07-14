library(GenomicFeatures)
txdb = makeTxDbFromGFF("database/hg38gtf/Homo_sapiens.GRCh38.114.gtf",format = "gtf")
exons_gene = exonsBy(txdb,by = "gene")
exons_gene_lens = lapply(exons_gene,function(x){sum(width(reduce(x)))})
exons_gene_lens[1:10]
gene_length = sapply(exons_gene_lens,function(x){x})
id_length = as.data.frame(gene_length)

# 在ensembl biomart下载geneid和name
name = read.table("database/hg38gtf/geneid_name.txt",sep = "\t",header = T)
name1 = name[match(rownames(id_length),name$Gene.stable.ID),]
id_length$gene_name = name1$Gene.name
id_length$gene_id = rownames(id_length)
write.table(id_length,"database/hg38gtf/genename_length.txt",row.names = F,sep = "\t",quote = F)
