cd ~/mouse
zcat gencode.vM33.basic.annotation.gtf | grep -v '^#' | awk -F "\t" '$3=="gene"{print $9}'  > protein_coding_gene.txt
awk -F ";" '{print $1"\t"$3}' protein_coding_gene.txt | sed 's/"//g;s/gene_id //;s/gene_name //;s/\.\S*//' > tmp.txt
awk 'BEGIN{OFS="\t"}{gsub(/[:.:].+/,"",$1);print $0}' tmp.txt > genesymbol.txt
awk -F "\t" 'NR==FNR{gene[$1]=$2; next} FNR==1{print "GeneSymbol\t"$0; next} NR>FNR{if ($1 in gene){print gene[$1]"\t"$0}}' genesymbol.txt gene_COUNT.txt > count.txt
#grep 'protein_coding' 
# ‘ 4933427D14Rik’, ‘ Pakap’
#删除 ENSMUSG00000107877 ENSMUSG00000089945 ENSMUSG00000090053