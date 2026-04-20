suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(openxlsx)
})

cat("=== DESeq2 Analysis ===\n")

# ── Load count matrix ──────────────────────────────────────────────────────────
counts_file <- "results/04_counts/counts_clean.txt"
cat("Loading counts from:", counts_file, "\n")

# Read header to get exact column names before R mangles them
all_lines   <- readLines(counts_file, n=2)
header_cols <- strsplit(all_lines[2], "\t")[[1]]
cell_ids    <- header_cols[7:length(header_cols)]  # skip Geneid + 5 annotation cols

# Read data with check.names=FALSE to preserve names
raw <- read.table(counts_file, header=TRUE, skip=1, row.names=1,
                  sep="\t", check.names=FALSE)
expr_cols  <- 6:ncol(raw)
count_mat  <- as.matrix(raw[, expr_cols])
storage.mode(count_mat) <- "integer"
colnames(count_mat) <- cell_ids

cat("Genes:", nrow(count_mat), "| Cells:", ncol(count_mat), "\n")
cat("First 3 cell IDs:", paste(colnames(count_mat)[1:3], collapse=", "), "\n")

# ── Load metadata ──────────────────────────────────────────────────────────────
meta_file <- "/Users/zachortega/Desktop/MS_SMARTseq_Pipeline/data/metadata.csv"
cat("Loading metadata from:", meta_file, "\n")
meta <- read.csv(meta_file, stringsAsFactors=FALSE)
rownames(meta) <- meta$cell_id
cat("First 3 metadata IDs:", paste(rownames(meta)[1:3], collapse=", "), "\n")

# Match cells
common_cells <- intersect(colnames(count_mat), rownames(meta))
cat("Cells with metadata:", length(common_cells), "\n")

if (length(common_cells) == 0) {
  cat("ERROR: No matching cell IDs!\n")
  cat("Count matrix cols (first 5):", paste(colnames(count_mat)[1:5], collapse=", "), "\n")
  cat("Metadata IDs (first 5):", paste(rownames(meta)[1:5], collapse=", "), "\n")
  stop("No matching cell IDs between count matrix and metadata.")
}

count_mat <- count_mat[, common_cells]
meta      <- meta[common_cells, ]
meta$condition <- relevel(factor(meta$condition), ref="CD45pos_MHCIIpos")
meta$strain    <- relevel(factor(meta$strain), ref="NOD")
meta$batch     <- factor(meta$batch)

# ── Filter low-count genes ─────────────────────────────────────────────────────
keep <- rowSums(count_mat) >= 10
cat("Genes after count filter (>=10 total):", sum(keep), "\n")
count_mat <- count_mat[keep, ]

# ── DESeq2 ─────────────────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = meta,
  design    = ~ batch + strain + condition
)
dds <- DESeq(dds, test="Wald", fitType="parametric")
cat("DESeq2 complete.\n")

# ── Results ────────────────────────────────────────────────────────────────────
out_dir <- "results/05_dge"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Condition: MHCIIlo vs CD45pos
res_lo <- results(dds, contrast=c("condition","CD45neg_MHCIIlo","CD45pos_MHCIIpos"), alpha=0.05)
res_hi <- results(dds, contrast=c("condition","CD45neg_MHCIIhi","CD45pos_MHCIIpos"), alpha=0.05)

df_lo <- as.data.frame(res_lo) %>% rownames_to_column("gene_id") %>% arrange(padj)
df_hi <- as.data.frame(res_hi) %>% rownames_to_column("gene_id") %>% arrange(padj)

sig_lo <- df_lo %>% filter(padj < 0.05, abs(log2FoldChange) >= 1)
sig_hi <- df_hi %>% filter(padj < 0.05, abs(log2FoldChange) >= 1)

cat("Significant DEGs (MHCIIlo vs CD45pos):", nrow(sig_lo), "\n")
cat("Significant DEGs (MHCIIhi vs CD45pos):", nrow(sig_hi), "\n")

# Save CSVs
write.csv(df_lo,  file.path(out_dir, "dge_MHCIIlo_vs_CD45pos.csv"),  row.names=FALSE)
write.csv(df_hi,  file.path(out_dir, "dge_MHCIIhi_vs_CD45pos.csv"),  row.names=FALSE)
write.csv(sig_lo, file.path(out_dir, "dge_sig_MHCIIlo_vs_CD45pos.csv"), row.names=FALSE)
write.csv(sig_hi, file.path(out_dir, "dge_sig_MHCIIhi_vs_CD45pos.csv"), row.names=FALSE)

# Excel workbook
wb <- createWorkbook()
addWorksheet(wb, "MHCIIlo_vs_CD45pos_all");   writeData(wb, 1, df_lo)
addWorksheet(wb, "MHCIIhi_vs_CD45pos_all");   writeData(wb, 2, df_hi)
addWorksheet(wb, "MHCIIlo_sig");              writeData(wb, 3, sig_lo)
addWorksheet(wb, "MHCIIhi_sig");              writeData(wb, 4, sig_hi)
saveWorkbook(wb, file.path(out_dir, "dge_results.xlsx"), overwrite=TRUE)

# ── VST for visualizations ─────────────────────────────────────────────────────
vst_data <- vst(dds, blind=FALSE)

# PCA
pca_data <- plotPCA(vst_data, intgroup=c("condition","strain"), returnData=TRUE)
pca_var  <- round(100 * attr(pca_data, "percentVar"))
p_pca <- ggplot(pca_data, aes(PC1, PC2, color=condition, shape=strain)) +
  geom_point(size=2.5, alpha=0.8) +
  xlab(paste0("PC1: ", pca_var[1], "% variance")) +
  ylab(paste0("PC2: ", pca_var[2], "% variance")) +
  theme_bw(base_size=12) +
  ggtitle("PCA — VST-normalized expression")
ggsave(file.path(out_dir, "pca_plot.pdf"), p_pca, width=8, height=6)
ggsave(file.path(out_dir, "pca_plot.png"), p_pca, width=8, height=6, dpi=200)

# Volcano — MHCIIlo vs CD45pos
make_volcano <- function(df, title, outname) {
  df <- df %>% mutate(sig=case_when(
    padj<0.05 & log2FoldChange>=1  ~ "Up",
    padj<0.05 & log2FoldChange<=-1 ~ "Down",
    TRUE ~ "NS"))
  top <- df %>% filter(sig!="NS") %>% slice_min(padj, n=20)
  p <- ggplot(df, aes(log2FoldChange, -log10(padj), color=sig)) +
    geom_point(size=0.8, alpha=0.6) +
    geom_text_repel(data=top, aes(label=gene_id), size=2.5, max.overlaps=15) +
    scale_color_manual(values=c(Up="firebrick",Down="steelblue",NS="grey70")) +
    geom_vline(xintercept=c(-1,1), linetype="dashed", color="grey50") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey50") +
    theme_bw(base_size=12) + ggtitle(title) + labs(color="Direction")
  ggsave(paste0(out_dir,"/",outname,".pdf"), p, width=7, height=5)
  ggsave(paste0(out_dir,"/",outname,".png"), p, width=7, height=5, dpi=200)
}
make_volcano(df_lo, "CD45neg MHCIIlo vs CD45pos", "volcano_MHCIIlo")
make_volcano(df_hi, "CD45neg MHCIIhi vs CD45pos", "volcano_MHCIIhi")

# Heatmap
all_sig <- union(head(sig_lo$gene_id,25), head(sig_hi$gene_id,25))
if (length(all_sig) >= 2) {
  mat_heat <- assay(vst_data)[all_sig, ]
  mat_heat <- t(scale(t(mat_heat)))
  ann_col  <- data.frame(condition=meta$condition, strain=meta$strain,
                          row.names=rownames(meta))
  pheatmap(mat_heat,
           annotation_col=ann_col,
           show_colnames=FALSE,
           fontsize_row=7,
           filename=file.path(out_dir, "heatmap_top_degs.pdf"),
           width=12, height=8)
}

cat("\n=== Analysis complete ===\n")
cat("Outputs in:", out_dir, "\n")
