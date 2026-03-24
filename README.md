# 🧬 Bacterial RNA-seq DESeq2 Pipeline

A downstream RNA-seq analysis pipeline for **bacterial transcriptomics**, integrating:

- Feature quantification (`featureCounts`)
- Differential expression analysis (`DESeq2`)
- Visualization (PCA + Volcano plot)

---

## 📌 Overview
This pipeline performs end-to-end RNA-seq analysis starting from BAM files to differential gene expression results. First, run the [Bulk_RNA_for_Long_reads](https://github.com/Zehui312/Bulk_RNA_for_Long_reads) pipeline to generate BAM files, which will be located in the `5_total_stat` directory.

## 1.💡Workflow
<img src="/img/workflow.png" width="500">


## 2. ⚙️ Create Environment

The running environment is the same as [Bulk_RNA_for_Long_reads](https://github.com/Zehui312/Bulk_RNA_for_Long_reads).


## 3. 📂Fill meta_data.csv

### 1. Metadata file (`meta_data.csv`)

Required columns:

| Column name | Description |
|------------|------------|
| sample_name | Project name |
| Output | Output directory |
| Bam_file_path | Path to BAM files |
| Sample_ID | Sample ID (used in DESeq2) |
| Group | Experimental group |
| GFF_file_path | GFF file |
| count_features | e.g. CDS,ncRNA,tmRNA,regulatory_region,oriT,oriC |
| gff_annotation | GFF attribute (e.g. ID) |
| minimum_fragment | featureCounts parameter -d |
| maximum_fragment | featureCounts parameter -D |
| padj_cutoff | e.g. 0.05 |
| log2fc_cutoff | e.g. 1.8 |

Optional:

- `Compare*` columns (e.g. `Compare1`, `Compare2`) for defining comparisons. You can add Compare3...

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

## 📊 Output
├── 1_feature_count
│   ├── JB251030_gene.count
│   ├── JB251030_gene.count.summary
│   └── rename_JB251030_gene.count
└── 2_deseq2
    ├── 0_gff_annotation.csv
    ├── 1_PCA_JB251030.pdf
    ├── P_vs_F1
    │   ├── 1_DESeq2_results_P_vs_F1.csv
    │   ├── 2_DESeq2_diff_genes_P_vs_F1.csv
    │   ├── 3_DESeq2_normalized_counts_P_vs_F1.csv
    │   └── P_vs_F1.pdf
    ├── P_vs_F2
    │   ├── 1_DESeq2_results_P_vs_F2.csv
    │   ├── 2_DESeq2_diff_genes_P_vs_F2.csv
    │   ├── 3_DESeq2_normalized_counts_P_vs_F2.csv
    │   └── P_vs_F2.pdf
