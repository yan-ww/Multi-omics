#!/usr/bin/env python3
# scripts/process_counts.py

import pandas as pd
import numpy as np
import os

def load_featurecounts_output(filepath):
    """加载featureCounts输出并清理"""
    # featureCounts输出有注释行以#开头
    df = pd.read_csv(filepath, sep='\t', comment='#')
    df = df.rename(columns={'Geneid': 'gene_id'})
    
    # 提取基因名称（如果需要）
    df['gene_name'] = df['gene_id'].str.split('.').str[0]
    
    # 获取counts列（以.bam结尾的列）
    counts_cols = [col for col in df.columns if col.endswith('.bam')]
    
    # 简化列名（去掉路径和.bam扩展名）
    sample_names = [os.path.basename(col).replace('.sorted.bam', '') for col in counts_cols]
    df_counts = df[['gene_id', 'gene_name'] + counts_cols].copy()
    df_counts.columns = ['gene_id', 'gene_name'] + sample_names
    
    return df_counts

def calculate_tpm(df_counts, gene_lengths=None):
    """计算TPM (Transcripts Per Million)"""
    if gene_lengths is None:
        # 如果没有基因长度文件，使用featureCounts的长度列
        # 需要从原始featureCounts输出中读取长度
        print("Warning: Using dummy gene lengths. For accurate TPM, provide gene lengths.")
        gene_lengths = np.ones(df_counts.shape[0])
    
    counts = df_counts.iloc[:, 2:].values
    # 每个基因的reads per kilobase (RPK)
    rpk = counts / (gene_lengths.reshape(-1, 1) / 1000)
    # 每个样本的缩放因子
    scaling_factors = rpk.sum(axis=0) / 1e6
    # TPM
    tpm = rpk / scaling_factors
    
    tpm_df = pd.DataFrame(tpm, 
                          index=df_counts['gene_id'], 
                          columns=df_counts.columns[2:])
    tpm_df.insert(0, 'gene_name', df_counts['gene_name'].values)
    tpm_df.insert(0, 'gene_id', df_counts['gene_id'].values)
    
    return tpm_df

def main():
    # 读取raw counts
    counts_raw = pd.read_csv(snakemake.input[0], sep='\t', comment='#')
    
    # 处理counts矩阵
    counts_matrix = load_featurecounts_output(snakemake.input[0])
    
    # 保存counts矩阵
    counts_matrix.to_csv(snakemake.output.counts_csv, index=False)
    
    # 计算TPM（简化版本，实际使用需要基因长度）
    # 这里使用featureCounts提供的长度（如果可用）
    if 'Length' in counts_raw.columns:
        gene_lengths = counts_raw['Length'].values
    else:
        gene_lengths = None
    
    tpm_matrix = calculate_tpm(counts_matrix, gene_lengths)
    tpm_matrix.to_csv(snakemake.output.tpm_csv, index=False)
    
    print(f"Counts matrix saved to {snakemake.output.counts_csv}")
    print(f"TPM matrix saved to {snakemake.output.tpm_csv}")

if __name__ == "__main__":
    main()
