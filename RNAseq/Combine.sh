# 1 mRNA层面的基因ID和基因symbol表达定量
mkdir -p ${out_dir}/02_geneCountAndTpm_cluster/mRNA/ && cd ${out_dir}/02_geneCountAndTpm_cluster/mRNA/
####step_1:run readcount_geneID_do
if [ ! -f readcount_geneID_do.SUCCESS ] ; then
    perl ${mnt_script_dir}/combine.readcount.pl \
    ${cfg_mRNA_genecount_list} \
    ./geneID_mRNA_readcount.xls \
    1>readcount_geneID_do.stdout 2>readcount_geneID_do.stderr \
    && touch readcount_geneID_do.SUCCESS

    cat readcount_geneID_do.stdout readcount_geneID_do.stderr
fi
if [ ! -f readcount_geneID_do.SUCCESS ] ; then >${out_dir}/02_geneCountAndTpm_cluster.FAIL ; exit -1 ; fi

####step2: readcount_geneSymbol_do
if [ ! -f readcount_geneSymbol_do.SUCCESS ] ; then
    perl ${mnt_script_dir}/change.name.v2.pl \
    ${gtf_path} \
    ./geneID_mRNA_readcount.xls \
    ./geneSymbol_mRNA_readcount.xls \
    1>readcount_geneSymbol_do.stdout 2>readcount_geneSymbol_do.stderr \
    && touch readcount_geneSymbol_do.SUCCESS

    cat readcount_geneSymbol_do.stdout readcount_geneSymbol_do.stderr
fi
if [ ! -f readcount_geneSymbol_do.SUCCESS ] ; then >${out_dir}/02_geneCountAndTpm_cluster.FAIL ; exit -1 ; fi

# 2 lncRNA层面的基因ID和基因symbol表达定量
mkdir -p ${out_dir}/02_geneCountAndTpm_cluster/lncRNA/ && cd ${out_dir}/02_geneCountAndTpm_cluster/lncRNA/
####step_1:run readcount_geneID_do
if [ ! -f readcount_geneID_do.SUCCESS ] ; then
    perl ${mnt_script_dir}/combine.readcount.pl \
    ${cfg_lncRNA_genecount_list} \
    ./geneID_lncRNA_readcount.xls \
    1>readcount_geneID_do.stdout 2>readcount_geneID_do.stderr \
    && touch readcount_geneID_do.SUCCESS

    cat readcount_geneID_do.stdout readcount_geneID_do.stderr
fi
if [ ! -f readcount_geneID_do.SUCCESS ] ; then >${out_dir}/02_geneCountAndTpm_cluster.FAIL ; exit -1 ; fi

####step2: readcount_geneSymbol_do
if [ ! -f readcount_geneSymbol_do.SUCCESS ] ; then
    perl ${mnt_script_dir}/change.name.v2.pl \
    ${gtf_path} \
    ./geneID_lncRNA_readcount.xls \
    ./geneSymbol_lncRNA_readcount.xls \
    1>readcount_geneSymbol_do.stdout 2>readcount_geneSymbol_do.stderr \
    && touch readcount_geneSymbol_do.SUCCESS

    cat readcount_geneSymbol_do.stdout readcount_geneSymbol_do.stderr
fi