library(DESeq2)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(optparse)

option_list <- list(
  make_option(c("-c", "--count_table"), type = "character", default = NULL,
              help = "Path to count table file", metavar = "character"),
  make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = "Path to metadata CSV file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output directory path", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$count_table) || is.null(opt$metadata) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Missing required arguments: --count_table, --metadata, and --output are required", call. = FALSE)
}

count_table_file <- opt$count_table
metadata_file <- opt$metadata
output_path <- opt$output



count_table <- read.table(count_table_file, row.names = 1, header = TRUE, comment.char = "#")
all_metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

sample_name <- all_metadata$sample_name[1]
gff_file <- all_metadata$GFF_file_path[1]
padj_cutoff <- all_metadata$padj_cutoff[1]
log2fc_cutoff <- all_metadata$log2fc_cutoff[1]
print("#================Input parameters:===============================")
print(paste("Count table file:", count_table_file))
print(paste("Metadata file:", metadata_file))
print(paste("GFF file:", gff_file))
print(paste("Output path:", output_path))
print(paste("Sample name:", sample_name))
print(paste("Adjusted p-value cutoff:", padj_cutoff))
print(paste("Log2 fold change cutoff:", log2fc_cutoff))
#=================================================================
#+++++++++++++++++++++++Step 0-1 read gff files ++++++++++++++++++
#=================================================================

gff <- rtracklayer::import(gff_file)
gff <- as.data.frame(gff)
gff_frame <- subset(gff, select = c("ID","Name","gene","start","end","strand","seqnames"))
gff_frame <- mutate(gff_frame, gene = if_else(is.na(gene), Name, gene))
# rownames(gff_frame) <- gff_frame$ID
# gff_frame$ID <- NULL
# head(gff_frame)



setwd(output_path)
write.csv(gff_frame, "0_gff_annotation.csv", row.names = FALSE)
#=================================================================
#+++++++++++++++++++++++Step 0-2 user-defined functions ++++++++++
#=================================================================
# Function to create volcano plot
create_volcano_plot <- function(res_gff, group1, group2, cutoff_fc = 1.8) {
    df <- res_gff
    volcano_name <- paste0(group1, "_vs_", group2)
    
    # Add significance column
    df$significance <- "No"
    df$significance[df$padj < 0.05 & df$log2FoldChange > cutoff_fc] <- "Up"
    df$significance[df$padj < 0.05 & df$log2FoldChange < -cutoff_fc] <- "Down"
    df$label <- ifelse(df$padj < 0.01 & abs(df$log2FoldChange) > cutoff_fc, df$gene, NA)
    diff_gene_count <- nrow(subset(df, df$padj < 0.05 & abs(df$log2FoldChange) > cutoff_fc))
    print(paste0("Number of differentially expressed genes between ", group1, " and ", group2, ": ", diff_gene_count))
    highlight_points <- subset(df, df$padj < 0.01 & abs(df$log2FoldChange) > cutoff_fc)
    max_padj <- max(-log10(df$padj[is.finite(-log10(df$padj))]))
    
    p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
        geom_point(alpha = 0.8, size = 2.5) +
        geom_point(data = highlight_points, aes(x = log2FoldChange, y = -log10(padj)),
                             size = 4, alpha = 0.5, shape = 21, stroke = 2) +
        geom_vline(xintercept = c(-cutoff_fc, cutoff_fc), linetype = "dashed", color = "black") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
        geom_text_repel(aes(label = label),
                                        size = 4,
                                        max.overlaps = 10,
                                        box.padding = 0.3,
                                        point.padding = 0.2) +
        scale_color_manual(values = c("Up" = "#d20815", "Down" = "#0072B2", "No" = "gray80")) +
        theme_bw() +
        theme(
            plot.title = element_text(size = 16, face = "bold"),
            axis.title = element_text(size = 15, face = "bold"),
            axis.text = element_text(size = 16),
            legend.position = "right",
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 14)
        ) +
        labs(
            title = paste0(volcano_name, ": labeling Padj < 0.05 and |log2FC| > ", cutoff_fc),
            subtitle = paste0("Number of differentially expressed genes: ", diff_gene_count),
            x = "log2(Fold Change)",
            y = "-log10(Adjusted P-value)",
            color = "Regulation"
        ) +
        xlim(c(-10, 10)) +
        ylim(c(0, max_padj)) +
        scale_x_continuous(breaks = seq(-9, 9, by = 1))
    
    ggsave(paste0(volcano_name, ".pdf"), plot = p, width = 10, height = 6)
    return(p)
}

#=================================================================
#+++++++++++++++++++++++Step 1 Count matrix input ++++++++++++++++
#=================================================================

metadata <- subset(all_metadata, select = c("Sample_ID", "Group"))
row.names(metadata) <- metadata$Sample_ID


#=================================================================
#+++++++++++++++++++++++Step 2 Check the row and col name ++++++++
#=================================================================
id <- rownames(metadata)
cts <- count_table[, id]

coldata <- metadata[id, ]

check_row_col <- all(rownames(coldata) == colnames(cts))
print(paste("Row names in coldata match column names in cts:",as.character(check_row_col)))

#=================================================================
#+++++++++++++++++++++++Step 3 Read count table ++++++++++++++++++
#=================================================================
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Group)

#=================================================================
#+++++++++++++++++++++++Step 3 Filter low count genes ++++++++++++
#=================================================================        
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds <- DESeq(dds)


#=================================================================
#+++++++++++++++++++++++Step 4 PCA plot +++++++++++++++++++++++++++
#=================================================================
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

plotPCA(vsd, intgroup=c("Group"))

pcaData <- plotPCA(vsd, intgroup=c("Group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))




p <- ggplot(pcaData, aes(PC1, PC2, color=Group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_bw() +theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.line  = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(color = "black",linewidth = 1.2, fill = NA)
  ) +
  ggtitle(paste("The PCA plot for", sample_name)) +
  coord_fixed()
ggsave(filename = paste0("1_PCA_", sample_name, ".pdf"), plot = p, width = 6, height = 5)


#=================================================================
#+++++++++++++++++++++++Step 4 Differential expression +++++++++++
#=================================================================
compare_metadata <- subset(all_metadata,select = c("Sample_ID","Group", grep("^Compare", colnames(all_metadata), value = TRUE)))
compare_data <- colnames(compare_metadata)[-c(1,2)]

for (i in 1:length(compare_data)) {
#>>Step 4-1: Subset the metadata for the current comparison  
  print(paste("Comparing groups for:", compare_data[i]))
  df_compare <- subset(compare_metadata, select = c("Sample_ID", "Group", compare_data[i]))
  df_compare <- df_compare[df_compare[[compare_data[i]]] == "yes", ]
  group1 <- unique(df_compare$Group)[1]
  group2 <- unique(df_compare$Group)[2]
  comp <- paste0(group1, "_vs_", group2)
  dir.create(paste0(output_path,"/", comp), recursive = TRUE)
  setwd(paste0(output_path,"/", comp) )
# >>Step 4-2: Run DESeq2 for the current comparison 
  res <- results(dds, contrast = c("Group", group1, group2))
  res_df <- as.data.frame(res)
  res_df$ID <- rownames(res_df)
  res_gff <- left_join(res_df, gff_frame, by = "ID")
  write.csv(res_gff, file = paste0("1_DESeq2_results_", comp, ".csv"), row.names = FALSE)
# >>Step 4-3: Subset the results for significant differentially expressed genes and save to a new CSV file
  diff_gene_tab <- subset(res_gff, padj < padj_cutoff & abs(log2FoldChange) > log2fc_cutoff)
  diff_gene_tab$regulation <- ifelse(diff_gene_tab$log2FoldChange > 0, "Up", "Down")
  
  write.csv(diff_gene_tab, file = paste0("2_DESeq2_diff_genes_", comp, ".csv"), row.names = FALSE)

# >>Step 4-4:  Save the normalized counts for all genes to a new CSV file, including gene annotations from the GFF file
  norm_counts <- counts(dds, normalized = TRUE)
  norm_counts <- norm_counts[,df_compare$Sample_ID ]
  norm_counts <- as.data.frame(norm_counts)
  norm_counts$ID <- rownames(norm_counts)
  norm_counts <- left_join(norm_counts, gff_frame, by = "ID")
  write.csv(norm_counts, file = paste0("3_DESeq2_normalized_counts_", comp, ".csv"), row.names = FALSE)

# >>Step 4-5: Create volcano plot
  create_volcano_plot(res_gff, group1, group2, cutoff_fc = log2fc_cutoff)
}

