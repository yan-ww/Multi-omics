#软件路径
fastqc=~{input.fastqc}
multiqc=~{input.multiqc}

#工作路径及变量
workdir=~{input.out_dir}
fq_config_file=~{input.fq_config_file}

if [ ! -d ${workdir}/log ]; then mkdir -p ${workdir}/log; fi
logfile=${workdir}/log/fastqclog1st.log

cat fq_config_file|while samplename file;do
	sample=$samplename
	logtime=$(date "+%Y-%m-%d %H:%M:%S")
	echo "------------------------------" >> $logfile
	echo "LOG TIME: $logtime" >> $logfile
	echo "FILE NAME: $sample" >> $logfile
	echo "------------------------------" >> $logfile
	cat $file|while read a b c d;do
		if [ ! -d ${workdir}/${sample} ]; then mkdir -p ${workdir}/${sample}; fi
		outdir=${workdir}/${sample}/fastqc
		if [ ! -d ${outdir} ]; then mkdir -p ${outdir}; fi
		$fastqc $file -o ${outdir} 2>>${logfile}
		echo "FastQC for $sample completed.\n" >> ${logfile}
	done
done

$multiqc ${workdir}/*/fastqc/ -o ${workdir} --fullnames raw_fastqc_report


