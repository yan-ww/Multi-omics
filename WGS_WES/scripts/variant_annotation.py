#!/usr/bin/env python3
"""
变异注释脚本 - 使用SnpEff和自定义过滤
Author: Your Name
Date: 2024
"""

import os
import subprocess
import pandas as pd
import gzip
import csv

def annotate_with_snpeff(input_vcf, output_vcf, db_name):
    """
    使用SnpEff进行变异注释
    """
    cmd = f"""
    snpeff ann {db_name} \
        {input_vcf} \
        -v \
        -noStats \
        -c snpEff.config \
        > {output_vcf}
    """
    
    try:
        subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
        # 索引注释后的VCF
        subprocess.run(f"tabix -p vcf {output_vcf}", shell=True, check=True)
        print("SnpEff annotation completed")
    except subprocess.CalledProcessError as e:
        print(f"Error in SnpEff annotation: {e}")
        raise

def parse_annotation_summary(vcf_file, output_csv):
    """
    解析注释后的VCF文件，提取关键信息生成摘要
    """
    variants = []
    
    # 打开VCF文件（支持gz压缩）
    if vcf_file.endswith('.gz'):
        f = gzip.open(vcf_file, 'rt')
    else:
        f = open(vcf_file, 'r')
    
    for line in f:
        if line.startswith('#'):
            continue
        
        fields = line.strip().split('\t')
        chrom = fields[0]
        pos = fields[1]
        ref = fields[3]
        alt = fields[4]
        qual = fields[5]
        info = fields[7]
        
        # 提取SnpEff注释
        ann_fields = {}
        if 'ANN=' in info:
            ann_part = [x for x in info.split(';') if x.startswith('ANN=')][0]
            ann = ann_part.split('=')[1].split('|')
            
            ann_fields = {
                'allele': ann[0],
                'annotation': ann[1],
                'impact': ann[2],
                'gene_name': ann[3],
                'gene_id': ann[4],
                'transcript': ann[6],
                'exon': ann[9],
                'protein_change': ann[11]
            }
        
        # 提取过滤状态
        filter_status = fields[6]
        
        variant_info = {
            'chromosome': chrom,
            'position': pos,
            'ref': ref,
            'alt': alt,
            'quality': qual,
            'filter': filter_status,
            **ann_fields
        }
        
        variants.append(variant_info)
    
    f.close()
    
    # 转换为DataFrame并保存
    df = pd.DataFrame(variants)
    
    # 按impact分类统计
    if 'impact' in df.columns:
        summary = df.groupby('impact').size().reset_index(name='count')
        print("\nVariant impact summary:")
        print(summary)
    
    # 保存到CSV
    df.to_csv(output_csv, index=False)
    print(f"Variant summary saved to {output_csv}")
    
    # 生成统计报告
    generate_statistics_report(df, output_csv.replace('.csv', '_stats.txt'))

def generate_statistics_report(df, output_file):
    """
    生成变异统计报告
    """
    with open(output_file, 'w') as f:
        f.write("=== Variant Statistics Report ===\n\n")
        
        f.write(f"Total variants: {len(df)}\n\n")
        
        if 'impact' in df.columns:
            f.write("Variants by impact:\n")
            impact_counts = df['impact'].value_counts()
            for impact, count in impact_counts.items():
                f.write(f"  {impact}: {count}\n")
        
        if 'annotation' in df.columns:
            f.write("\nTop variant types:\n")
            anno_counts = df['annotation'].value_counts().head(10)
            for anno, count in anno_counts.items():
                f.write(f"  {anno}: {count}\n")
        
        if 'gene_name' in df.columns:
            f.write("\nTop affected genes:\n")
            gene_counts = df['gene_name'].value_counts().head(20)
            for gene, count in gene_counts.items():
                if gene and gene != '.':
                    f.write(f"  {gene}: {count}\n")
        
        f.write("\nQuality distribution:\n")
        f.write(f"  Mean QUAL: {df['quality'].astype(float).mean():.2f}\n")
        f.write(f"  Median QUAL: {df['quality'].astype(float).median():.2f}\n")

def main():
    # 从Snakemake获取参数
    input_vcf = snakemake.input.vcf
    output_vcf = snakemake.output.annotated_vcf
    output_csv = snakemake.output.summary_csv
    db_name = snakemake.params.snpeff_db
    
    # 执行注释
    annotate_with_snpeff(input_vcf, output_vcf, db_name)
    
    # 解析结果并生成摘要
    parse_annotation_summary(output_vcf, output_csv)
    
    print("Variant annotation pipeline completed successfully!")

if __name__ == "__main__":
    main()