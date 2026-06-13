#!/usr/bin/env python3
"""标记基因分析"""

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

def main(adata_input, markers_output, heatmap_output, dotplot_output, 
         n_markers, min_logfc):
    # 读取数据
    adata = sc.read_h5ad(adata_input)
    
    print(f"Finding marker genes for {len(adata.obs['leiden'].unique())} clusters...")
    
    # 计算每个cluster的标记基因
    sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon',
                           key_added='rank_genes')
    
    # 提取标记基因结果
    marker_results = []
    for cluster in adata.obs['leiden'].unique():
        cluster_markers = sc.get.rank_genes_groups_df(adata, group=cluster)
        cluster_markers = cluster_markers.head(n_markers)
        cluster_markers['cluster'] = cluster
        marker_results.append(cluster_markers)
    
    all_markers = pd.concat(marker_results, ignore_index=True)
    all_markers.to_csv(markers_output, index=False)
    print(f"Saved marker genes to {markers_output}")
    
    # 获取top标记基因列表
    top_markers = []
    for cluster in adata.obs['leiden'].unique():
        cluster_top = all_markers[all_markers['cluster'] == cluster].head(3)
        top_markers.extend(cluster_top['names'].tolist())
    top_markers = list(dict.fromkeys(top_markers))  # 去重
    
    # 热图
    print("Generating marker gene heatmap...")
    fig, ax = plt.subplots(figsize=(12, 8))
    sc.pl.matrixplot(adata, top_markers, groupby='leiden', 
                    dendrogram=True, colorbar_title='Mean expression',
                    standard_scale='var', cmap='Reds',
                    show=False, ax=ax)
    plt.tight_layout()
    plt.savefig(heatmap_output, dpi=150, bbox_inches='tight')
    plt.close()
    
    # 点图
    print("Generating dotplot...")
    fig, ax = plt.subplots(figsize=(12, 6))
    sc.pl.dotplot(adata, top_markers[:20], groupby='leiden',
                 dendrogram=True, colorbar_title='Expression',
                 cmap='viridis', show=False, ax=ax)
    plt.tight_layout()
    plt.savefig(dotplot_output, dpi=150, bbox_inches='tight')
    plt.close()
    
    # 打印每个cluster的top标记基因
    print("\nTop marker genes per cluster:")
    for cluster in adata.obs['leiden'].unique():
        cluster_markers = all_markers[all_markers['cluster'] == cluster]
        top_genes = cluster_markers.head(5)['names'].tolist()
        print(f"  Cluster {cluster}: {', '.join(top_genes)}")

if __name__ == "__main__":
    main(snakemake.input.adata,
         snakemake.output.markers,
         snakemake.output.heatmap,
         snakemake.output.dotplot,
         snakemake.params.n_markers,
         snakemake.params.min_logfc)