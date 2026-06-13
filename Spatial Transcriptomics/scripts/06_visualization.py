#!/usr/bin/env python3
# 综合可视化报告

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys

# 读取参数
input_rds = sys.argv[1]
input_normalized = sys.argv[2]
input_clustered = sys.argv[3]
input_svg = sys.argv[4]
input_deconv = sys.argv[5]
output_pdf = sys.argv[6]
sample_name = sys.argv[7]

# 创建多页PDF
from matplotlib.backends.backend_pdf import PdfPages

with PdfPages(output_pdf) as pdf:
    
    # 页1: 标题
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.text(0.5, 0.5, f'Spatial Transcriptomics Analysis Report\nSample: {sample_name}',
            ha='center', va='center', fontsize=20)
    ax.axis('off')
    pdf.savefig(fig)
    plt.close()
    
    # 页2: QC统计
    fig, ax = plt.subplots(figsize=(10, 6))
    qc_data = pd.DataFrame({
        'Metric': ['Spots after QC', 'Genes per spot', 'UMIs per spot'],
        'Value': [1500, 3500, 12000]  # 示例数值
    })
    ax.bar(qc_data['Metric'], qc_data['Value'])
    ax.set_title('Quality Control Summary')
    ax.set_ylabel('Count')
    plt.xticks(rotation=45)
    pdf.savefig(fig)
    plt.close()
    
    # 页3: SVG总结
    fig, ax = plt.subplots(figsize=(10, 6))
    svg_data = pd.read_csv(input_svg)
    top_svg = svg_data.head(20)
    ax.barh(range(len(top_svg)), top_svg['morans_I'])
    ax.set_yticks(range(len(top_svg)))
    ax.set_yticklabels(top_svg['gene'])
    ax.set_xlabel("Moran's I")
    ax.set_title('Top 20 Spatially Variable Genes')
    pdf.savefig(fig)
    plt.close()
    
    # 页4: 细胞类型组成
    fig, ax = plt.subplots(figsize=(10, 6))
    deconv_data = pd.read_csv(input_deconv, index_col=0)
    celltype_means = deconv_data.mean().sort_values(ascending=True)
    ax.barh(celltype_means.index, celltype_means.values)
    ax.set_xlabel('Average Proportion')
    ax.set_title('Cell Type Composition')
    pdf.savefig(fig)
    plt.close()
    
    # 页5: 方法总结
    fig, ax = plt.subplots(figsize=(10, 6))
    methods = ['Normalization: SCT', 'Clustering resolution: 0.5', 
               'SVG detection: SPARK', 'Deconvolution: RCTD']
    for i, method in enumerate(methods):
        ax.text(0.1, 0.9 - i*0.1, method, fontsize=12, transform=ax.transAxes)
    ax.axis('off')
    ax.set_title('Analysis Methods Used')
    pdf.savefig(fig)
    plt.close()

print(f"Report generated: {output_pdf}")