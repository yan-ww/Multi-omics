# 从ensemble官网下载hg19或hg38的gtf文件
cat  Homo_sapiens.GRCh38.114.gtf | grep -v '^#' | awk -F "\t" '$3=="gene"{print $9}'|grep 'protein_coding'  > protein_coding_gene.txt
awk -F ";" '{print $1"\t"$3}' protein_coding_gene.txt |grep 'gene_name' | sed 's/"//g;s/gene_id //;s/gene_name //;s/\.\S*//' > tmp.txt
awk 'BEGIN{OFS="\t"}{gsub(/[:.:].+/,"",$1);print $0}' tmp.txt > genesymbol.txt
awk -F "\t" 'BEGIN{OFS="\t"} {sub(/^ /, "", $2)} 1' genesymbol.txt > geneid_name_fromgtf_mRNA.txt
sed -i '/ENSG00000280987/d;/ENSG00000285188/d;/ENSG00000284762/d;/ENSG00000285458/d;/ENSG00000285437/d;/ENSG00000254093/d' geneid_name_fromgtf_mRNA.txt
awk -F "\t" 'NR==FNR{gene[$1]=$2; next} FNR==1{print "GeneSymbol\t"$0; next} NR>FNR{if ($1 in gene){print gene[$1]"\t"$0}}' geneid_name_fromgtf_mRNA.txt gene_COUNT.txt > count.txt
rm protein_coding_gene.txt tmp.txt genesymbol.txt
# 一些重复基因
# 小鼠
# ‘ 4933427D14Rik’, ‘ Pakap’
#删除 ENSMUSG00000107877 ENSMUSG00000089945 ENSMUSG00000090053

# HUMAN
# "MATR3"   "PDE4C"   "PDE8B"   "C4orf36" "POLR2J3" "PINX1"
# 对应删除ENSG00000280987  ENSG00000285188  ENSG00000284762  ENSG00000285458  ENSG00000285437  ENSG00000254093