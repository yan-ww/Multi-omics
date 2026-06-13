#!/usr/bin/env python3
"""标准化和高变基因选择"""

import scanpy as sc
import matplotlib.pyplot as plt

def main(adata_input, adata_output, hvg_plot, n_top_genes):
    # 读取数据
    adata = sc.read_h5ad(adata_input)
    
    # 标准化
    print("Performing normalization...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # 识别高变基因
    print(f"Identifying top {n_top_genes} highly variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key=None)
    
    # 可视化高变基因
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # 高变基因分布
    sc.pl.highly_variable_genes(adata, show=False, ax=axes[0])
    axes[0].set_title('Highly Variable Genes')
    
    # 标准化前后的对比
    axes[1].hist(adata.X.flatten(), bins=50, alpha=0.7)
    axes[1].set_xlabel('Log-normalized expression')
    axes[1].set_ylabel('Frequency')
    axes[1].set_title('Expression Distribution After Normalization')
    
    plt.tight_layout()
    plt.savefig(hvg_plot, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"Number of HVGs: {adata.var['highly_variable'].sum()}")
    
    # 保存数据
    adata.write(adata_output)
    print(f"Saved normalized data to {adata_output}")

if __name__ == "__main__":
    main(snakemake.input.adata,
         snakemake.output.adata,
         snakemake.output.hvg_plot,
         snakemake.params.n_top_genes)