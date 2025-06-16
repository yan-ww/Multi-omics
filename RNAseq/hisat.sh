#! /bin/bash

#软件路径
hisat=~{input.hisat2}
samtools=~{input.samtools}

#工作路径及变量
workdir=~{input.out_dir}
fq_config_file=~{input.fq_config_file}
index=~{input.hisat2_index}

if [ ! -d ${workdir}/log ]; then mkdir -p ${workdir}/log; fi
logfile=${workdir}/log/hisat2.log


cat fq_config_file|while samplename R1 R2;do
	logtime=$(date "+%Y-%m-%d %H:%M:%S")
	echo "------------------------------" >> $logfile
	echo "LOG TIME: $logtime" >> $logfile
	echo "FILE NAME: $samplename" >> $logfile
	echo "------------------------------" >> $logfile
	if [ ! -d ${workdir}/${samplename} ]; then mkdir -p ${workdir}/${samplename}; fi
	outdir=${workdir}/${samplename}/Align
	if [ ! -d ${outdir} ]; then mkdir -p ${outdir}; fi
	if [ ! -d ${workdir}/log ]; then mkdir -p ${workdir}/log; fi
	$hisat -p 6 -x {index} \
		-1 ${R1} \
		-2 ${R2} \
		-S ${outdir}/${samplename}.sam \
		2 >> $logfile
	$samtools view -bS ${outdir}/${samplename}.sam > ${outdir}/${samplename}.bam
	$samtools sort -@ 1 ${outdir}/${samplename}.bam -o ${outdir}/${samplename}.sorted.bam
	$samtools index -@ 1 ${outdir}/${samplename}.sorted.bam		
	files=$(ls -1v ${outdir}/*.sorted.bam 2>/dev/null) # 列出所有排序后的BAM文件
	printf "%s\t%s\n" "$samplename" "$(echo "$files" | tr '\n' '\t' | sed 's/\t$//')" >>${workdir}/bam.lst # 将样本名和文件路径写入bam.lst，使用制表符分隔
	echo "------------end---------------" >>${logfile}
done



