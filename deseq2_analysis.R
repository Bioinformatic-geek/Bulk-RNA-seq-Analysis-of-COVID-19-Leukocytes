# DESeq2 Differential Expression Analysis
# Dataset: GSE157103 — COVID-19 vs Healthy leukocytes
# ============================================================

# Load libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
# ============================================================
# STEP 1 — Load and clean the count matrix
# ============================================================

# Read the featureCounts output
# comment.char skips the first line which is just the command log
raw_counts <- read.table(
  "counts/raw_counts.txt",
  header = TRUE,
  sep = "\t",
  comment.char = "#"
)
# Keep only the gene ID and count columns
# Columns 1-6 are gene info (Chr, Start, End etc), columns 7-12 are counts
counts <- raw_counts[, c(1, 7:12)]
# Clean up column names — remove the long file path, keep just the SRR number
colnames(counts) <- c("gene_id",
                       "COVID_419", "COVID_420", "COVID_421",
                       "Healthy_527", "Healthy_528", "Healthy_529")

# Set gene IDs as row names (DESeq2 requires this)
rownames(counts) <- counts$gene_id
counts <- counts[, -1]  # remove the gene_id column now it's the row name
# Preview — should show 6 columns and many rows
head(counts)
dim(counts)  # shows rows x columns

# ============================================================
# STEP 2 — Create metadata table
# ============================================================

# This tells DESeq2 which samples are which condition
metadata <- data.frame(
  sample = c("COVID_419", "COVID_420", "COVID_421",
             "Healthy_527", "Healthy_528", "Healthy_529"),
  condition = c("COVID", "COVID", "COVID",
                "Healthy", "Healthy", "Healthy")
)
# Set sample names as row names
rownames(metadata) <- metadata$sample
metadata <- metadata[, -1, drop = FALSE]
# Make condition a factor — tells DESeq2 Healthy is the reference group
# (fold changes will be calculated as COVID vs Healthy)
metadata$condition <- factor(metadata$condition,
                              levels = c("Healthy", "COVID"))

# Verify samples match between counts and metadata
all(rownames(metadata) == colnames(counts))  # must print TRUE
# ============================================================
# STEP 3 — Create DESeq2 object and run analysis
# ============================================================

# Create the DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ condition   # compare by condition (COVID vs Healthy)
)
# Pre-filter — remove genes with very low counts across all samples
# Keeps genes with at least 10 counts in at least 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
cat("Genes remaining after filtering:", nrow(dds), "\n")

# Run DESeq2 — this one function does normalisation + statistics
dds <- DESeq(dds)

# ============================================================
# STEP 4 — Extract results
# ============================================================

# Get results table
res <- results(dds,
               contrast = c("condition", "COVID", "Healthy"),
               alpha = 0.05)   # significance threshold
# Summary of results
summary(res)
# Convert to dataframe for easier handling
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
# Add significance labels
res_df$significant <- "Not significant"
res_df$significant[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
res_df$significant[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"
# Count significant genes
table(res_df$significant)
# Top upregulated genes
cat("\nTOP 20 UPREGULATED GENES IN COVID:\n")
top_up <- res_df[res_df$significant == "Upregulated", ]
top_up <- top_up[order(top_up$log2FoldChange, decreasing = TRUE), ]
print(head(top_up[, c("gene_id", "log2FoldChange", "padj")], 20))
# Top downregulated genes
cat("\nTOP 20 DOWNREGULATED GENES IN COVID:\n")
top_down <- res_df[res_df$significant == "Downregulated", ]
top_down <- top_down[order(top_down$log2FoldChange, decreasing = FALSE), ]
print(head(top_down[, c("gene_id", "log2FoldChange", "padj")], 20))
# Save full results to file
write.csv(res_df,
          "counts/deseq2_results.csv",
          row.names = FALSE)
cat("\nFull results saved to deseq2_results.csv\n")
# ============================================================
# STEP 5 — Visualisations
# ============================================================

# --- Plot 1: Volcano Plot ---
# Shows every gene — x axis is fold change, y axis is significance
# Upregulated genes appear top right, downregulated top left

res_df_plot <- res_df[!is.na(res_df$padj), ]  # remove NA rows
volcano <- ggplot(res_df_plot,
                  aes(x = log2FoldChange,
                      y = -log10(padj),
                      color = significant)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not significant" = "grey"
  )) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot: COVID-19 vs Healthy",
    x = "log2 Fold Change (COVID / Healthy)",
    y = "-log10 adjusted p-value",
    color = "Expression"
  ) +
  theme_minimal()
ggsave("counts/volcano_plot.png",
       volcano, width = 8, height = 6, dpi = 150)
cat("Volcano plot saved\n")
# --- Plot 2: Heatmap of top 30 genes ---
# Shows expression pattern across all 6 samples
top_genes <- c(
  head(rownames(top_up), 15),
  head(rownames(top_down), 15)
)
# Get normalised counts for these genes
vsd <- vst(dds, blind = FALSE)  # variance stabilising transformation
mat <- assay(vsd)[top_genes, ]
# Scale each gene (so colour shows relative change, not absolute counts)
mat_scaled <- t(scale(t(mat)))
# Column labels
png("counts/heatmap.png",
    width = 800, height = 1000)
pheatmap(mat_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         annotation_col = metadata,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Top Differentially Expressed Genes\nCOVID-19 vs Healthy")
dev.off()
cat("Heatmap saved\n")
# --- Plot 3: PCA Plot ---
# Shows whether samples cluster by condition — healthy together, COVID together
pca_plot <- plotPCA(vsd, intgroup = "condition") +
  geom_label_repel(aes(label = name), size = 3) +
  labs(title = "PCA Plot: Sample Clustering by Condition") +
  theme_minimal()

ggsave("counts/pca_plot.png",
       pca_plot, width = 7, height = 5, dpi = 150)
cat("PCA plot saved\n")
cat("\n=== DESeq2 Analysis Complete ===\n")
cat("Output files in ~/rnaseq_project/counts/:\n")
cat("  deseq2_results.csv  — full statistical results\n")
cat("  volcano_plot.png    — volcano plot\n")
cat("  heatmap.png         — expression heatmap\n")
cat("  pca_plot.png        — PCA clustering plot\n")
