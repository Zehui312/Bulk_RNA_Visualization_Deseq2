# 🧬 Bacterial RNA-seq DESeq2 Pipeline

This pipeline performs end-to-end RNA-seq analysis from BAM files to differential gene expression results. First, run the [Bulk_RNA_for_Long_reads](https://github.com/Zehui312/Bulk_RNA_for_Long_reads) pipeline to generate BAM files in the `5_total_stat` directory.

**Features:**
- Feature quantification (`featureCounts`)
- Differential expression analysis (`DESeq2`)
- Visualization (PCA + Volcano plots)

---


## 1. 💡 Workflow

<img src="/img/workflow.png" width="500">

## 2. ⚙️ Create Environment

Use the same environment as [Bulk_RNA_for_Long_reads](https://github.com/Zehui312/Bulk_RNA_for_Long_reads).

## 3. 📂 Fill meta_data.csv

### Metadata File (`meta_data.csv`)

**Required columns:**

| Column | Description |
|--------|-------------|
| sample_name | Project name |
| Output | Output directory |
| Bam_file_path | Path to BAM files |
| Sample_ID | Sample ID for DESeq2 |
| Group | Experimental group |
| GFF_file_path | GFF file path |
| count_features | Features to count (e.g., CDS, ncRNA, tmRNA, regulatory_region, oriT, oriC) |
| gff_annotation | GFF attribute (e.g., ID) |
| minimum_fragment | featureCounts parameter `-d` |
| maximum_fragment | featureCounts parameter `-D` |
| padj_cutoff | Adjusted p-value threshold (e.g., 0.05) |
| log2fc_cutoff | Log2 fold-change threshold (e.g., 1.8) |

**Optional:**
- `Compare1`, `Compare2`, `Compare3`, etc. for defining comparisons

## 🚀 4. Usage

After filling **meta_data.csv**, run:

```bash
bash Run_Deseq2_pipeline.sh
```

## 📊 Output

The pipeline generates:

```
├── 1_feature_count/
│   ├── {sample}_gene.count
│   ├── {sample}_gene.count.summary
│   └── rename_{sample}_gene.count
└── 2_deseq2/
    ├── 0_gff_annotation.csv
    ├── 1_PCA_{sample}.pdf
    ├── {comparison1}/
    │   ├── 1_DESeq2_results_{comparison1}.csv
    │   ├── 2_DESeq2_diff_genes_{comparison1}.csv
    │   ├── 3_DESeq2_normalized_counts_{comparison1}.csv
    │   └── {comparison1}.pdf
    └── {comparison2}/
        ├── 1_DESeq2_results_{comparison2}.csv
        ├── 2_DESeq2_diff_genes_{comparison2}.csv
        ├── 3_DESeq2_normalized_counts_{comparison2}.csv
        └── {comparison2}.pdf
```

**Key outputs:**
- `1_feature_count/` – Raw feature counts
- `2_deseq2/` – DESeq2 results, PCA plots, and volcano plots for each comparison
