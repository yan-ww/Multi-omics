workflow RNAseq{
    ### pathways of analysis software
    String fastqc
    String multiqc
    String trim_glore
    String docker

    ### index files
    String homo_gtf

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
            script_dir=script_dir,
            homo_gtf=homo_gtf
    }

    call Combine{
        input:
            cfg_mRNA_genecount_list=Stringtie.cfg_mRNA_genecount_list

    }
    
    call PCA{
        input:
            sample_info=sample_info,
            allRNA_genecount=Combine.allRNA_genecount,
            out_dir=out_dir,
            script_dir=script_dir
    }

    call DESeq2{
        input:
            sample_info=sample_info,
            mRNA_genecount=Combine.mRNA_genecount,
            lncRNA_genecount=Combine.lncRNA_genecount,
            allRNA_genecount=Combine.allRNA_genecount,
            out_dir=out_dir,
            script_dir=script_dir
    }

    call GSEA{
        input:
            DEGresult=DESeq2.DEGresult,
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

    runtime{
		docker:"gmseq/bioinfonda:v1.0.5"
		cpu:"1"
		memory:"1G"
		num_proc:"1"
		maxRetries:2
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
    runtime{
		docker:"gmseq/bioinfonda:v1.0.5"
		cpu:"1"
		memory:"1G"
		num_proc:"1"
		maxRetries:2
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
    runtime{
		docker:"gmseq/bioinfonda:v1.0.5"
		cpu:"1"
		memory:"1G"
		num_proc:"1"
		maxRetries:2
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
    runtime{
		docker:"gmseq/bioinfonda:v1.0.5"
		cpu:"1"
		memory:"1G"
		num_proc:"1"
		maxRetries:2
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
    runtime{
		docker:"gmseq/bioinfonda:v1.0.5"
		cpu:"1"
		memory:"1G"
		num_proc:"1"
		maxRetries:2
    }
}

task Stringtie{
    String stringtie
    String bam_list
    String out_dir
    String script_dir
    String homo_gtf
    command <<<
        set -e
        sh ${script_dir}/stringtie.sh
    >>>
    output{
        String gtf="${out_dir}/${sample_name}/08_bam2allTrascript/${sample_name}_transcripts.gtf"
    }
    runtime{
		docker:"gmseq/bioinfonda:v1.0.5"
		cpu:"1"
		memory:"1G"
		num_proc:"1"
		maxRetries:2
    }
}

task Combine{
    String cfg_mRNA_genecount_list
    command <<<
        set -e
        sh ${script_dir}/Combine.sh
    >>>
    output{
        String mRNA_genecount=
        String lncRNA_genecount=
        String allRNA_genecount=
    }
    runtime{
		docker:"gmseq/bioinfonda:v1.0.5"
		cpu:"1"
		memory:"1G"
		num_proc:"1"
		maxRetries:2
    }
}

task PCA{
    String sample_info
    String allRNA_genecount
    String out_dir
    String script_dir
    command <<<
        set -e
        Rscript ${script_dir}/PCA.R ${out_dir} ${allRNA_genecount}
    >>>
    output{

    }
    runtime{
		docker:"gmseq/bioinfonda:v1.0.5"
		cpu:"1"
		memory:"1G"
		num_proc:"1"
		maxRetries:2
    }
}

task DESeq2{
    String sample_info
    String mRNA_genecount
    String lncRNA_genecount
    String allRNA_genecount
    String out_dir
    String script_dir
    command <<<
        set -e
        Rscript ${script_dir}/DifferentialExpressionAnalysis.R ${out_dir} ${mRNA_genecount}
    >>>
    output{
        String DEGresult=
    }
    runtime{
		docker:"gmseq/bioinfonda:v1.0.5"
		cpu:"1"
		memory:"1G"
		num_proc:"1"
		maxRetries:2
    }
}

task GSEA{
    String DEGresult
    String out_dir
    String script_dir
    command <<<
        set -e
        Rscript ${script_dir}/GSEA.R ${out_dir} ${DEGresult}
    >>>
    output{
        String gseago=
    }
    runtime{
		docker:"gmseq/bioinfonda:v1.0.5"
		cpu:"1"
		memory:"1G"
		num_proc:"1"
		maxRetries:2
    }
}