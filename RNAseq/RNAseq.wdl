workflow RNAseq{
    ### pathways of analysis software
    String fastqc
    String multiqc
    String trim_glore

    ###run modules
    Boolean run_allStep_boolean=false
	Boolean run_bam2TPM_boolean=true
    Boolean run_pca=if (run_allStep_boolean || run_pca_boolean || run_deggene_boolean || run_deg2pathway_boolean) then true else false

    ### pathways of input and output
    fq_config_file=fq_config_file
    ## fq_config_file has four columns，including samplename,group (case or control), fq1 path, fq2 path
    String out_dir
    String script_dir

    call Task_Allocation{
        input:
            fq_config_file=fq_config_file
    }

    call Fastqc{
        input:
            fq_config_file=Task_Allocation.cfg_task_Allocation.list,
            fastqc=fastqc,
            multiqc=multiqc,
            out_dir=out_dir,
            script_dir=script_dir      
    }

    call Trim_glore{
        input:
            fq_config_file=Task_Allocation.cfg_task_Allocation.list,
            trim=trim_glore,
            out_dir=out_dir,
            script_dir=script_dir
    }

    call Fastqc2nd{
        input:
            fq_config_file=Trim_glore.clean_fq_list,
            fastqc=fastqc,
            multiqc=multiqc,
            out_dir=out_dir,
            script_dir=script_dir
    }

    call Hisat2{
        input:
           fq_config_file=Trim_glore.clean_fq_list,
           hisat2=hisat2,
           samtools=samtools,
           hisat2_index=hisat2_index,
           out_dir=out_dir,
           script_dir=script_dir
    }

    call Stringtie{
        input:
            stringtie=stringtie,
            bam_list=Hisat2.bam_list,
            out_dir=out_dir,
            script_dir=script_dir
    }
}

task Task_Allocation{
    String fq_config_file

    command <<<
        set -e
        cat ${fq_config_file} | cut -f 1 | sort -u | while read sampleName ; do
			mkdir -p ${out_dir}/$sampleName/00.cfg_dir
			cat ${fq_config_file} | sort -u | awk -v n=$sampleName '$1==n' > ${out_dir}/$sampleName/00.cfg_dir/cfg_task.$sampleName
			echo -e "$sampleName\t${out_dir}/$sampleName/00.cfg_dir/cfg_task.$sampleName"
		done |sort -u > ${out_dir}/cfg_task_Allocation.list
	>>>
    output{
        String allocated_list=${out_dir}/cfg_task_Allocation.list
    }
}

task Fastqc{
    String fq_config_file
    String fastqc
    String multiqc
    String out_dir
    String script_dir

    command <<<
        set -e
        sh ${script_dir}/fastqc.sh
    >>>
    output{
        raw_fq_report=${out_dir}/raw_fastqc_report.html
    }
}

task Trim_glore{
    String fq_config_file
    String trim
    String out_dir
    String script_dir

    command <<<
        set -e
        sh ${script_dir}/trim_glore.sh
    >>>
    output{
        clean_fq_list=${out_dir}/cleanfq.lst
    }
}

task Fastqc2n{
    String fq_config_file
    String fastqc
    String multiqc
    String out_dir
    String script_dir

    command <<<
        set -e
        sh ${script_dir}/fastqc2nd.sh
    >>>
    output{
        clean_fq_report=${out_dir}/clean_fastqc_report.html
    }
}

task Hisat2{
    String fq_config_file
    String hisat2
    String samtools
    String hisat2_index
    String out_dir
    String script_dir

    command <<<
        set -e
        sh ${script_dir}/hisat2.sh
    >>>
    output{
        bam_list=${out_dir}/bam.lst
    }
}

task Stringtie{
    String stringtie
    String bam_list
    String out_dir
    String script_dir
    command <<<
        set -e
        sh ${script_dir}/stringtie.sh
    >>>
    output{
        
    }
}