#!/bin/bash

#软件路径
trim=~{input.trim_glore}

#工作路径及变量
workdir=~{input.out_dir}
logfile=${workdir}/log/trim_galore.log
fq_config_file=~{input.fq_config_file}

rm -rf ${workdir}/cleanfq.lst
cat fq_config_file|while samplename file;do
	logtime=$(date "+%Y-%m-%d %H:%M:%S")
	echo "------------------------------" >> $logfile
	echo "FILENAME: ${samplename}" >> $logfile
	echo "LOG TIME: $logtime" >> $logfile
	echo "------------------------------" >> $logfile
	cat $file|while read a b f1 f2;do
		if [ ! -d ${workdir}/${samplename} ]; then mkdir -p ${workdir}/${samplename}; fi
		outdir=${workdir}/${samplename}/clean
		if [ ! -d ${outdir} ]; then mkdir -p ${outdir}; fi
		if [ ! -d ${workdir}/log ]; then mkdir -p ${workdir}/log; fi
		echo "Processing: $f1 and $f2" >> ${logfile}
		${trim} -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --paired -o ${outdir} $f1 $f2
		files=$(ls -1v ${outdir}/* 2>/dev/null) # 列出所有处理后的文件
		printf "%s\t%s\n" "$samplename" "$(echo "$files" | tr '\n' '\t' | sed 's/\t$//')" >>${workdir}/cleanfq.lst # 将样本名和文件路径写入cleanfq.lst，使用制表符分隔
		echo "Trimmed: $f1 and $f2" >> ${logfile}
		echo "------------------------------" >> $logfile	
	done
done




	




