#!/usr/bin/env python3
"""最终可视化"""

import scanpy as sc
import matplotlib.pyplot as plt

def main(adata_input, umap_output, marker_output):
    # 读取数据
    adata = sc.read_h5ad(adata_input)
    
    # 最终UMAP图（更高分辨率）
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 聚类
    sc.pl.umap(adata, color=['leiden'], ax=axes[0, 0],
              title='UMAP - Cell Clusters', show=False,
              legend_loc='on data', legend_fontsize=10,
              frameon=False)
    
    # 基因数量
    sc.pl.umap(adata, color=['n_genes_by_counts'], ax=axes[0, 1],
              title='UMAP - Genes per Cell', show=False,
              color_map='viridis', frameon=False)
    
    # 线粒体比例
    sc.pl.umap(adata, color=['pct_counts_mt'], ax=axes[1, 0],
              title='UMAP - Mitochondrial %', show=False,
              color_map='plasma', frameon=False)
    
    # 总表达量
    sc.pl.umap(adata, color=['total_counts'], ax=axes[1, 1],
              title='UMAP - Total Expression', show=False,
              color_map='inferno', frameon=False)
    
    plt.suptitle(f'Final UMAP Visualization\n{adata.n_obs} cells, {len(adata.obs["leiden"].unique())} clusters', 
                 fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(umap_output, dpi=200, bbox_inches='tight')
    plt.close()
    
    # 选择一些已知的PBMC标记基因进行可视化
    known_markers = ['CD3D', 'CD14', 'MS4A1', 'CD8A', 'FCGR3A', 
                    'NKG7', 'CD79A', 'CD4', 'CD68', 'CD1C']
    
    # 过滤出实际存在的基因
    existing_markers = [g for g in known_markers if g in adata.var_names]
    
    if existing_markers:
        # 标记基因UMAP图
        n_plots = len(existing_markers)
        n_cols = min(3, n_plots)
        n_rows = (n_plots + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3*n_rows))
        if n_plots == 1:
            axes = [axes]
        else:
            axes = axes.flatten()
        
        for i, gene in enumerate(existing_markers):
            sc.pl.umap(adata, color=[gene], ax=axes[i],
                      title=f'{gene}', show=False,
                      color_map='Reds', frameon=False,
                      size=20)
            axes[i].set_xlabel('')
            axes[i].set_ylabel('')
        
        # 隐藏多余的子图
        for i in range(len(existing_markers), len(axes)):
            axes[i].axis('off')
        
        plt.suptitle('Marker Gene Expression on UMAP', fontsize=14, y=1.02)
        plt.tight_layout()
        plt.savefig(marker_output, dpi=200, bbox_inches='tight')
        plt.close()
    
    print("Visualization complete!")

if __name__ == "__main__":
    main(snakemake.input.adata,
         snakemake.output.umap_by_cluster,
         snakemake.output.marker_genes_plot)