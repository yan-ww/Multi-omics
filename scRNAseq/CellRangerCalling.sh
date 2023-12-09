############prepare data################

#download raw RNA-seq dataset from SRA through prefetch
prefetch SRR10004175
#convert .sra file to fastq format
fastq-dump --split-files --gzip SRR8795637.sra -O ./
#check info of reads in fastq.gz files
zless SRR8795637_1.fastq.gz
#rename all fastq.gz files to meet the need of cell ranger
mv SRR8795637_1.fastq.gz SRR8795637_S1_L001_I1_001.fastq.gz
mv SRR8795637_2.fastq.gz SRR8795637_S1_L001_R2_001.fastq.gz
mv SRR8795637_3.fastq.gz SRR8795637_S1_L001_R3_001.fastq.gz

############cell ranger################

#download mm10 genome & its gtf file
wget http://ftp.ensembl.org/pub/release-105/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz
wget http://ftp.ensembl.org/pub/release-105/gtf/mus_musculus/Mus_musculus.GRCm39.105.gtf.gz
#check extra lines in gtf
cat Mus_musculus.GRCm39.105.gtf | grep -v "#" | awk -v FS='gene_biotype' 'NF>1{print $2}' | awk -F ";" '{print $1}' | sort | uniq -c  
#filter them 
cellranger mkgtf \
Mus_musculus.GRCm39.105.gtf \
Mus_musculus.GRCm39.105.filtered.gtf \
--attribute=gene_biotype:protein_coding \
--attribute=gene_biotype:lincRNA \
--attribute=gene_biotype:antisense \
--attribute=gene_biotype:IG_LV_gene \
--attribute=gene_biotype:IG_V_gene \
--attribute=gene_biotype:IG_V_pseudogene \
--attribute=gene_biotype:IG_D_gene \
--attribute=gene_biotype:IG_J_gene \
--attribute=gene_biotype:IG_J_pseudogene \
--attribute=gene_biotype:IG_C_gene \
--attribute=gene_biotype:IG_C_pseudogene \
--attribute=gene_biotype:TR_V_gene \
--attribute=gene_biotype:TR_V_pseudogene \
--attribute=gene_biotype:TR_D_gene \
--attribute=gene_biotype:TR_J_gene \
--attribute=gene_biotype:TR_J_pseudogene \
--attribute=gene_biotype:TR_C_gene \
#construct reference file 
cellranger mkref \
--genome=Homosapiens.GRCh39.105 \
--fasta=Mus_musculus.GRCm39.dna_sm.primary_assembly.fa \
--genes=Mus_musculus.GRCm39.105.filtered.gtf
#gene expression matrix calling
cellranger count \
--id=run_count_mysample 
--transcriptome=/omics/Homosapiens.GRCh39.105 \
--fastqs=/omics/SRR8795637/ 
#check output from cell ranger
firefox /run_count_mysample/outs/web_summary.html












