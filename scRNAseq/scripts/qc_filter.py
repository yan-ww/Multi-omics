#!/usr/bin/env python3
"""质控过滤"""

import scanpy as sc
import pandas as pd
import numpy as np

def main(adata_input, adata_output, qc_report, sample, 
         min_genes, max_genes, max_mito, min_cells):
    # 读取数据
    adata = sc.read_h5ad(adata_input)
    
    print(f"Before filtering: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # 记录过滤前的统计
    stats_before = {
        'total_cells': adata.n_obs,
        'total_genes': adata.n_vars,
        'mean_genes_per_cell': adata.obs['n_genes_by_counts'].mean(),
        'median_genes_per_cell': adata.obs['n_genes_by_counts'].median(),
        'mean_mito_percent': adata.obs['pct_counts_mt'].mean()
    }
    
    # 细胞过滤
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, max_genes=max_genes)
    adata = adata[adata.obs['pct_counts_mt'] < max_mito, :]
    
    # 基因过滤
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    print(f"After filtering: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # 记录过滤后的统计
    stats_after = {
        'total_cells': adata.n_obs,
        'total_genes': adata.n_vars,
        'mean_genes_per_cell': adata.obs['n_genes_by_counts'].mean(),
        'median_genes_per_cell': adata.obs['n_genes_by_counts'].median(),
        'mean_mito_percent': adata.obs['pct_counts_mt'].mean()
    }
    
    # 保存QC报告
    qc_stats = pd.DataFrame({
        'Metric': list(stats_before.keys()),
        'Before': list(stats_before.values()),
        'After': list(stats_after.values())
    })
    qc_stats.to_csv(qc_report, index=False)
    print(f"Saved QC report to {qc_report}")
    
    # 保存过滤后的数据
    adata.write(adata_output)
    print(f"Saved filtered data to {adata_output}")

if __name__ == "__main__":
    main(snakemake.input.adata,
         snakemake.output.adata,
         snakemake.output.qc_report,
         snakemake.params.sample,
         snakemake.params.min_genes,
         snakemake.params.max_genes,
         snakemake.params.max_mito,
         snakemake.params.min_cells)