#!/usr/bin/env python3
"""加载数据和初始QC"""

import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
from anndata import AnnData

# 设置
sc.settings.figdir = "results/figures"
sc.set_figure_params(dpi=80, figsize=(8, 6))

def main(adata_input, adata_output, qc_plot, sample, min_genes, max_mito):
    # 读取数据
    print(f"Loading data for sample: {sample}")
    adata = sc.read_h5ad(adata_input)
    
    # 基本统计
    print(f"原始数据: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # 计算线粒体基因比例
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, 
                               log1p=False, inplace=True)
    
    # 初始QC可视化
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # 小提琴图
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                jitter=0.4, multi_panel=True, show=False, ax=axes[0])
    axes[0].set_title(f'QC Metrics - {sample}')
    
    # 散点图
    sc.pl.scatter(adata, 'total_counts', 'n_genes_by_counts', 
                 color='pct_counts_mt', show=False, ax=axes[1])
    axes[1].set_title('Genes vs Counts')
    
    # 分布图
    axes[2].hist(adata.obs['n_genes_by_counts'], bins=50, alpha=0.7)
    axes[2].axvline(min_genes, color='r', linestyle='--', label=f'min={min_genes}')
    axes[2].set_xlabel('Genes per cell')
    axes[2].set_ylabel('Frequency')
    axes[2].legend()
    axes[2].set_title('Gene distribution')
    
    plt.tight_layout()
    plt.savefig(qc_plot, dpi=150, bbox_inches='tight')
    plt.close()
    
    # 保存原始数据
    adata.write(adata_output)
    print(f"Saved raw data to {adata_output}")

if __name__ == "__main__":
    main(snakemake.input.h5, 
         snakemake.output.adata,
         snakemake.output.qc_plot,
         snakemake.params.sample,
         snakemake.params.min_genes,
         snakemake.params.max_mito)