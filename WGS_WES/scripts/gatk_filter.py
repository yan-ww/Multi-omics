#!/usr/bin/env python3
"""
GATK变异过滤脚本
Author: Your Name
Date: 2024
"""

import os
import subprocess
import gzip

def filter_snps(input_vcf, output_snps, output_indels, params):
    """
    使用GATK VariantFiltration过滤SNPs和INDELs
    """
    # SNPs过滤条件
    snp_filters = [
        f"QD < {params['qd']}",
        f"FS > {params['fs']}",
        f"MQ < {params['mq']}",
        f"MQRankSum < {params['mqr']}",
        f"SOR > {params['sor']}"
    ]
    
    # INDELs过滤条件
    indel_filters = [
        f"QD < {params['qd']}",
        f"FS > {params['fs']}",
        f"SOR > {params['sor']}"
    ]
    
    # 过滤SNPs
    cmd_snp = f"""
    gatk VariantFiltration \
        -V {input_vcf} \
        -O {output_snps} \
        --filter-name "SNP_filter" \
        --filter-expression "{' && '.join(snp_filters)}" \
        --genotype-filter-name "GQ_filter" \
        --genotype-filter-expression "GQ < 20"
    """
    
    # 过滤INDELs
    cmd_indel = f"""
    gatk VariantFiltration \
        -V {input_vcf} \
        -O {output_indels} \
        --filter-name "INDEL_filter" \
        --filter-expression "{' && '.join(indel_filters)}" \
        --genotype-filter-name "GQ_filter" \
        --genotype-filter-expression "GQ < 20"
    """
    
    try:
        subprocess.run(cmd_snp, shell=True, check=True, executable='/bin/bash')
        subprocess.run(cmd_indel, shell=True, check=True, executable='/bin/bash')
        print("Variant filtering completed successfully")
    except subprocess.CalledProcessError as e:
        print(f"Error in variant filtering: {e}")
        raise

def main():
    # 获取Snakemake提供的输入输出
    input_vcf = snakemake.input.raw_vcf
    output_snps = snakemake.output.filtered_snps
    output_indels = snakemake.output.filtered_indels
    params = snakemake.params
    
    filter_snps(input_vcf, output_snps, output_indels, params)

if __name__ == "__main__":
    main()