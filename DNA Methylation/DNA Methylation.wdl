workflow Methylation{
    ### pathways of analysis software
    String fastp
    String Bismark

    ### index files

    ### run modules

    ### pathways of input and output
    fq_config_file=fq_config_file
    ## fq_config_file has four columns，including samplename,group (case or control), fq1 path, fq2 path
    String out_dir
    String script_dir

    call Task_Allocation{
        input:
            fq_config_file=fq_config_file
    }
}