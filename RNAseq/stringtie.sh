#! /bin/bash

#软件和索引文件路径
stringtie=~{input.stringtie}
homo_gtf_path=~{input.homo_gtf}
script_dir=~{input.script_dir}

#工作路径及变量
workdir=~{input.out_dir}

sample_name=~{input.sample_name}
sample_bam=~{sample_bam}
logfile=${workdir}/log/stringtie.log

# 主要步骤：stringtie计数--> mRNA和lncRNA分离--->提取count值
$stringtie \
    -G ${homo_gtf_path} \
    -p 4 \
    -B \
    -e \
    -l ${sample_name} \
    -o ${sample_name}_transcripts.gtf \
    -A ${sample_name}_geneAbundance.tpm \
    ${sample_bam} \
    1>bam2allTranscript_do.stdout 2>bam2allTranscript_do.stderr \
    && touch bam2allTranscript_do.SUCCESS

# 判断是否有链特异性。链特异性的样本需要使用--rf参数，需要在bamqc的时候计算
# if [ ! -f bam2allTranscript_do.SUCCESS ] ; then
#     strandlib=$(awk -v a=${sample_name} '$1==a{print $2}' ${all_strand})
#     if [[ "$strandlib" == "true" ]] ; then
#         /usr/local/bin/stringtie \
#         -G ${homo_gtf_path} \
#         -p 4 \
#         -B \
#         -e \
#         --rf \
#         -l ${sample_name} \
#         -o ${sample_name}_transcripts.gtf \
#         -A ${sample_name}_geneAbundance.tpm \
#         ${sample_bam} \
#         1>bam2allTranscript_do.stdout 2>bam2allTranscript_do.stderr \
#         && touch bam2allTranscript_do.SUCCESS
#     else
#         /usr/local/bin/stringtie \
#         -G ${homo_gtf_path} \
#         -p 4 \
#         -B \
#         -e \
#         -l ${sample_name} \
#         -o ${sample_name}_transcripts.gtf \
#         -A ${sample_name}_geneAbundance.tpm \
#         ${sample_bam} \
#         1>bam2allTranscript_do.stdout 2>bam2allTranscript_do.stderr \
#         && touch bam2allTranscript_do.SUCCESS
#     fi

# 分离mRNA和lncRNA
python3 ${script_dir}/split_coding_nocoding.py \
			-t ${sample_name}_geneAbundance.tpm \
			-g ${sample_name}_transcripts.gtf \
			-b ${Gene2Biotype_path} \
			-s ${sample_name}

####step 2
		awk '$NF>1' ${sample_name}_lncRNA.tpm > lncRNA_gene_expression.tsv
		awk '$NF>1' ${sample_name}_mRNA.tpm > mRNA_gene_expression.tsv

		perl  ${script_dir}/transcript_sample_stat.pl -g ${Gene2Biotype_path} -n ${sample_name}_transcripts.gtf -s ${sample_name} -o transcript_classfication.xls

		echo -e "${sample_name}\t${out_dir}/${sample_name}/${sample_name}_mRNA.tpm" > cfg_mRNA_tpm.list
		echo -e "${sample_name}\t${out_dir}/${sample_name}/${sample_name}_lncRNA.tpm" > cfg_lncRNA_tpm.list

		echo -e "${sample_name} bam2TPM_do_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> bam2TPM_do.time

echo -e ${sample_name}"\t"${gtf_file_mRNA} > cfg_mRNA_gtf.list

read_len="100"
if [ -f ${out_dir}/../${sample_name}/02_fqqc/${sample_name}/${sample_name}.json ] ; then
    read_len=$(cat ${out_dir}/../${sample_name}/02_fqqc/${sample_name}/${sample_name}.json|grep read1_mean_length|cut -d ':' -f 2|cut -d ',' -f 1|head -1)
fi

####step_1:run tpm2readcount_do
if [ ! -f tpm2readcount_do.SUCCESS ] ; then
    ${script_dir}/prepDE.py \
    -i cfg_mRNA_gtf.list \
    -l $read_len \
    1>tpm2readcount_do.stdout 2>tpm2readcount_do.stderr \
    && touch tpm2readcount_do.SUCCESS

    cat tpm2readcount_do.stdout tpm2readcount_do.stderr
fi
if [ ! -f tpm2readcount_do.SUCCESS ] ; then > ${out_dir}/${sample_name}/mRNA.FAIL ; exit -1 ; fi

echo -e "${sample_name}\t${out_dir}/${sample_name}/mRNA/transcript_count_matrix.csv" > cfg_mRNA_transcript_count.list
echo -e "${sample_name}\t${out_dir}/${sample_name}/mRNA/gene_count_matrix.csv" > cfg_mRNA_genecount.list

# 提取转录本为单位的FPKM和TPM
if [ ! -f gtf2tpm_fpkm.SUCCESS ] ; then
    /usr/local/bin/miniconda/miniconda3/bin/python ${script_dir}/transcript_FPKM_TPM.py ${gtf_file_mRNA} ${sample_name} && touch gtf2tpm_fpkm.SUCCESS
fi
if [ ! -f gtf2tpm_fpkm.SUCCESS ] ; then exit -1 ; fi


# 2 进行lncRNA分析 -------------------------------------------------------------------------------------
if [ -f ${out_dir}/${sample_name}/lncRNA.FAIL ] ; then rm ${out_dir}/${sample_name}/lncRNA.FAIL ; fi
mkdir -p ${out_dir}/${sample_name}/lncRNA/ && cd ${out_dir}/${sample_name}/lncRNA/

echo -e ${sample_name}"\t"${gtf_file_lncRNA} > cfg_lncRNA_gtf.list
if [ ! -f tpm2readcount_do.SUCCESS ] ; then
    ${script_dir}/prepDE.py \
    -i cfg_lncRNA_gtf.list \
    -l $read_len \
    1>tpm2readcount_do.stdout 2>tpm2readcount_do.stderr \
    && touch tpm2readcount_do.SUCCESS

    cat tpm2readcount_do.stdout tpm2readcount_do.stderr
fi
if [ ! -f tpm2readcount_do.SUCCESS ] ; then > ${out_dir}/${sample_name}/lncRNA.FAIL ; exit -1 ; fi

echo -e "${sample_name}\t${out_dir}/${sample_name}/lncRNA/transcript_count_matrix.csv" > cfg_lncRNA_transcript_count.list
echo -e "${sample_name}\t${out_dir}/${sample_name}/lncRNA/gene_count_matrix.csv" > cfg_lncRNA_genecount.list

# 提取转录本为单位的FPKM和TPM
if [ ! -f gtf2tpm_fpkm.SUCCESS ] ; then
    python ${script_dir}/transcript_FPKM_TPM.py ${gtf_file_lncRNA} ${sample_name} && touch gtf2tpm_fpkm.SUCCESS
fi
if [ ! -f gtf2tpm_fpkm.SUCCESS ] ; then exit -1 ; fi

