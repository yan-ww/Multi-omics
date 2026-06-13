#!/usr/bin/env python3
"""PCA、UMAP和聚类"""

import scanpy as sc
import matplotlib.pyplot as plt

def main(adata_input, adata_output, pca_plot, umap_plot, 
         n_pcs, n_neighbors, resolution):
    # 读取数据
    adata = sc.read_h5ad(adata_input)
    
    # 只使用高变基因进行PCA
    print(f"Running PCA with {n_pcs} components...")
    sc.tl.pca(adata, n_comps=n_pcs, use_highly_variable=True)
    
    # PCA方差解释率图
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # 方差解释率
    axes[0].plot(range(1, len(adata.uns['pca']['variance_ratio']) + 1), 
                 adata.uns['pca']['variance_ratio'], 'bo-')
    axes[0].set_xlabel('PC Components')
    axes[0].set_ylabel('Variance Explained Ratio')
    axes[0].set_title('PCA Variance Explained')
    axes[0].grid(True, alpha=0.3)
    
    # 累积方差解释率
    cumsum = np.cumsum(adata.uns['pca']['variance_ratio'])
    axes[1].plot(range(1, len(cumsum) + 1), cumsum, 'ro-')
    axes[1].set_xlabel('PC Components')
    axes[1].set_ylabel('Cumulative Variance Explained')
    axes[1].set_title('PCA Cumulative Variance')
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(pca_plot, dpi=150, bbox_inches='tight')
    plt.close()
    
    # 计算邻域图
    print(f"Computing neighborhood graph (n_neighbors={n_neighbors})...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    
    # UMAP
    print("Running UMAP...")
    sc.tl.umap(adata)
    
    # Leiden聚类
    print(f"Running Leiden clustering (resolution={resolution})...")
    sc.tl.leiden(adata, resolution=resolution)
    
    print(f"Number of clusters found: {len(adata.obs['leiden'].unique())}")
    
    # UMAP可视化
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # UMAP按聚类着色
    sc.pl.umap(adata, color=['leiden'], ax=axes[0], 
              title='UMAP - Leiden Clusters', show=False,
              legend_loc='on data', legend_fontsize=8)
    
    # UMAP按基因数量着色
    sc.pl.umap(adata, color=['n_genes_by_counts'], ax=axes[1],
              title='UMAP - Genes per Cell', show=False,
              color_map='viridis')
    
    plt.tight_layout()
    plt.savefig(umap_plot, dpi=150, bbox_inches='tight')
    plt.close()
    
    # 保存聚类信息
    cluster_counts = adata.obs['leiden'].value_counts().sort_index()
    print("Cluster sizes:")
    for cluster, count in cluster_counts.items():
        print(f"  Cluster {cluster}: {count} cells")
    
    # 保存数据
    adata.write(adata_output)
    print(f"Saved clustered data to {adata_output}")

if __name__ == "__main__":
    import numpy as np
    main(snakemake.input.adata,
         snakemake.output.adata,
         snakemake.output.pca_plot,
         snakemake.output.umap_plot,
         snakemake.params.n_pcs,
         snakemake.params.n_neighbors,
         snakemake.params.resolution)