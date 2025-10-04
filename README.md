# Multi-omics Analysis Pipeline

This repository provides a collection of workflows for **multi-omics data analysis**, integrating transcriptomics, genomics, DNA methylation, and other molecular layers. It is designed to illustrate systematic approaches for dissecting biological complexity and can serve as a reference for computational biology research.

**This repository is under continuous development.**
---

## Overview

The repository includes step-by-step pipelines for different omics data types, with a focus on:

* **Transcriptomics (RNA-seq)**
* **Genomics (variant analysis)**
* **Epigenomics (DNA methylation)**

Each pipeline is modular and reproducible, implemented mainly in **R, Python, and Linux-based workflows**. The methods are motivated by research questions in systems biology, disease mechanism exploration, and biomarker discovery.

---

## Analysis Workflows

### 1. Transcriptomics

* Preprocessing and quality control of raw RNA-seq data
* Alignment and quantification of gene expression
* Differential expression analysis
* Functional enrichment (GO, KEGG, GSEA)
* Visualization: heatmaps, PCA, volcano plots

### 2. Genomics

* Quality control and trimming of DNA sequencing reads
* Alignment to the reference genome
* Variant calling and filtering (SNPs, indels)
* Annotation of variants with functional impact
* Summary statistics and visualization

### 3. DNA Methylation

* Preprocessing and normalization of methylation arrays or bisulfite sequencing
* Differential methylation analysis (CpG sites, regions)
* Integration with gene expression changes
* Visualization: methylation heatmaps, boxplots, Manhattan plots

---

## Requirements

The workflows make use of the following tools and environments:

* **Programming:** Python, R, Linux shell
* **Workflow management:** WDL
* **Packages:** DESeq2, edgeR, clusterProfiler, GATK, Bismark, scikit-learn, TensorFlow, PyTorch

Please refer to each subfolder for specific installation and usage instructions.

---

## Repository Structure

```
Multi-omics/
│── transcriptomics/    # RNA-seq analysis scripts
│── genomics/           # Variant calling and genomic analysis
│── methylation/        # DNA methylation workflows
│── utils/              # Helper scripts and functions
```

---

## Citation

If you find this repository useful, please cite it or contact me for collaboration.

---

## Contact

**Yanping Weng**

Email: yanping.weng@outlook.com

LinkedIn: https://www.linkedin.com/in/wengyanping/




