RNA-seq Differential Expression Analysis with DESeq2

This repository contains an R script to perform differential gene expression analysis using DESeq2 on RNA-seq count data, including handling of technical replicates, collapsing samples, and identifying up- and down-regulated genes.

Contents

rna_seq_analysis.R — Main R script for:

Reading raw count matrices

Creating metadata/sample tables

Collapsing technical replicates by donor/condition

Running DESeq2 differential expression analysis

Generating up- and down-regulated gene lists

combined_count_matrix.tsv — Example input count matrix (or path to your actual counts)

README.md — This file

Dependencies

The script requires the following R packages:

DESeq2

You can install missing packages via:

install.packages(c("ggplot2", "RColorBrewer", "pheatmap"))
BiocManager::install(c("DESeq2", "org.Hs.eg.db"))

Usage

Set working directories

Modify the paths in the script to point to your count matrix and featureCounts output:

setwd('/path/to/working/directory/')
directory <- '/path/to/feature_counts/'


Prepare metadata

Your metadata table (combined_count_matrix.tsv) should include:

sample	file	group	donor	Run

group: experimental condition

donor: biological replicate

Run: technical replicate/run ID

Run the script

The script will:

Collapse technical replicates (collapseReplicates)

Run DESeq2 differential expression analysis

Identify up- and down-regulated genes (based on adjusted p-value and log2 fold change)

Save results to separate TSV files:

upregulated_genes_DSC_vs_MZVH.tsv
downregulated_genes_DSC_vs_MZVH.tsv

Output

upregulated_genes_*.tsv — Genes with log2FoldChange > 0 and padj ≤ threshold

downregulated_genes_*.tsv — Genes with log2FoldChange < 0 and padj ≤ threshold

Optional outputs if added:

Variance-stabilized counts (assay(vsd))

QC plots (PCA, heatmaps)

Notes

Ensure that condition factor levels are correctly set; the reference level is automatically set to "DSC_control" in this workflow.

Low-count filtering is recommended to remove genes with negligible counts before DESeq2 analysis.

Adjust the significance threshold (alpha) based on your control conditions.

