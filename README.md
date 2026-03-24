# 🧬 Bacterial RNA-seq DESeq2 Pipeline

A reproducible RNA-seq analysis pipeline for **bacterial transcriptomics**, integrating:

- Feature quantification (`featureCounts`)
- Differential expression analysis (`DESeq2`)
- Visualization (PCA + Volcano plot)
- GFF-based gene annotation integration

---

## 📌 Overview

This pipeline performs end-to-end RNA-seq analysis starting from BAM files to differential gene expression results.

### Workflow


---

## ⚙️ Requirements



## 📂 Input Files

### 1. Metadata file (`meta_data.csv`)

Required columns:

| Column name | Description |
|------------|------------|
| sample_name | Project/sample name |
| Sample_ID | Sample ID (used in DESeq2) |
| Group | Experimental group |
| Bam_file_path | Path to BAM files |
| Output | Output directory |
| minimum_fragment | featureCounts parameter |
| maximum_fragment | featureCounts parameter |
| count_features | e.g. CDS |
| gff_annotation | GFF attribute (e.g. ID) |
| GFF_file_path | GFF file |
| padj_cutoff | e.g. 0.05 |
| log2fc_cutoff | e.g. 1 |

Optional:

- `Compare*` columns (e.g. `Compare1`, `Compare2`) for defining comparisons.

---

### 2. BAM files

- Sorted BAM files
- Must match `Sample_ID`

---



## 🚀 Usage

### Step 1: Run pipeline

```bash
bash Run_Deseq2_pipeline.sh
```