#! /bin/bash

#软件路径
fastqc=~{input.fastqc}
multiqc=~{input.multiqc}

#工作路径及变量
workdir=~{input.out_dir}
fq_config_file=~{input.fq_config_file}

if [ ! -d ${workdir}/log ]; then mkdir -p ${workdir}/log; fi
logfile=${workdir}/log/fastqclog2st.log

cat fq_config_file|while samplename R1 R2;do
	logtime=$(date "+%Y-%m-%d %H:%M:%S")
	echo "------------------------------" >> $logfile
	echo "LOG TIME: $logtime" >> $logfile
	echo "FILE NAME: $samplename" >> $logfile
	echo "------------------------------" >> $logfile
	if [ ! -d ${workdir}/${samplename} ]; then mkdir -p ${workdir}/${samplename}; fi
	outdir=${workdir}/${samplename}/fastqc2nd
	if [ ! -d ${outdir} ]; then mkdir -p ${outdir}; fi
	$fastqc $R1 -o ${outdir} 2>>${logfile}
	$fastqc $R2 -o ${outdir} 2>>${logfile}
	echo "FastQC for $samplename completed.\n" >> ${logfile}

done 

$multiqc ${workdir}/*/fastqc2nd/ -o ${workdir} --fullnames clean_fastqc_report


