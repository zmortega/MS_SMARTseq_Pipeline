suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

# ── Paths ──────────────────────────────────────────────────────────────────────
base_dir  <- "/Users/zachortega/Desktop/MS_SMARTseq_Pipeline"
# -- Self-copy into pipeline scripts/ for version control ----------------------
# Copies this script into the pipeline folder so it is tracked by git.
# To run on a different machine, update base_dir below.
local({
  this_script <- normalizePath(
    grep("--file=", commandArgs(trailingOnly=FALSE), value=TRUE) |>
      sub("--file=", "", x=_),
    mustWork=FALSE
  )
  if (length(this_script) == 1 && nchar(this_script) > 0) {
    dest_dir <- file.path(base_dir, "scripts")
    dir.create(dest_dir, showWarnings=FALSE, recursive=TRUE)
    dest <- file.path(dest_dir, "nod_plots.R")
    if (file.copy(this_script, dest, overwrite=TRUE)) {
      cat("Script copied to:", dest, "\n")
    } else {
      cat("Warning: could not copy script to scripts/\n")
    }
  }
})

counts_f  <- file.path(base_dir, "results/04_counts/counts_clean.txt")
meta_f    <- file.path(base_dir, "data/metadata.csv")
out_dir   <- file.path(base_dir, "results/05_dge/NOD_plots")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# ── Load data ──────────────────────────────────────────────────────────────────
counts <- read.table(counts_f, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
meta   <- read.csv(meta_f, stringsAsFactors=FALSE)
rownames(meta) <- meta$cell_id

# ── Ensembl → gene symbol map ──────────────────────────────────────────────────
sym_file <- "/tmp/ensembl_to_symbol.txt"
if (file.exists(sym_file)) {
  id2sym  <- read.table(sym_file, sep="\t", col.names=c("ensembl","symbol"),
                        stringsAsFactors=FALSE)
  sym_map <- setNames(id2sym$symbol, id2sym$ensembl)
  cat("Symbol map loaded:", length(sym_map), "entries\n")
} else {
  gtf_f <- file.path(base_dir, "reference/gencode.vM33.primary_assembly.annotation.gtf")
  cat("Building symbol map from GTF (one-time, may take ~30s)...\n")
  gtf_lines  <- readLines(gtf_f)
  gene_lines <- gtf_lines[grepl('\tgene\t', gtf_lines)]
  ensembl    <- sub('.*gene_id "([^"]+)".*', '\\1', gene_lines)
  symbol     <- sub('.*gene_name "([^"]+)".*', '\\1', gene_lines)
  sym_map    <- setNames(symbol, ensembl)
  write.table(data.frame(ensembl, symbol), sym_file, sep="\t",
              row.names=FALSE, col.names=FALSE, quote=FALSE)
  cat("Symbol map built:", length(sym_map), "genes\n")
}

# Helper: return gene symbol if available, otherwise keep Ensembl ID
to_sym <- function(ens) {
  sym <- sym_map[ens]
  ifelse(!is.na(sym) & sym != "", sym, ens)
}

# ── Filter to NOD only ─────────────────────────────────────────────────────────
nod_cells <- meta$cell_id[meta$strain == "NOD"]
shared    <- intersect(colnames(counts), nod_cells)
counts    <- counts[, shared]
meta      <- meta[shared, ]
cat("NOD cells:", ncol(counts), "\n")
cat("Conditions:", paste(table(meta$condition), names(table(meta$condition)), collapse=" | "), "\n")

# ── Filter low-count genes ─────────────────────────────────────────────────────
keep   <- rowSums(counts) >= 10
counts <- counts[keep, ]
cat("Genes passing filter:", nrow(counts), "\n")

# ── DESeq2 (NOD only, single batch so no batch term) ──────────────────────────
meta$condition <- relevel(factor(meta$condition), ref="CD45pos_MHCIIpos")

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData   = meta,
  design    = ~ condition
)
dds <- estimateSizeFactors(dds, type="poscounts")
dds <- DESeq(dds, test="Wald", fitType="parametric", quiet=TRUE)
cat("DESeq2 complete.\n")

# VST for violin plots
vsd  <- vst(dds, blind=FALSE)
expr <- assay(vsd)  # genes x cells

# ── Define contrasts ───────────────────────────────────────────────────────────
contrasts <- list(
  list(
    name   = "CD45pos_MHCIIpos_vs_CD45neg_MHCIIhi",
    label  = "NOD CD45+ MHCIIpos vs CD45- MHCIIhi",
    con    = c("condition", "CD45neg_MHCIIhi", "CD45pos_MHCIIpos"),
    groups = c("CD45pos_MHCIIpos", "CD45neg_MHCIIhi")
  ),
  list(
    name   = "CD45pos_MHCIIpos_vs_CD45neg_MHCIIlo",
    label  = "NOD CD45+ MHCIIpos vs CD45- MHCIIlo",
    con    = c("condition", "CD45neg_MHCIIlo", "CD45pos_MHCIIpos"),
    groups = c("CD45pos_MHCIIpos", "CD45neg_MHCIIlo")
  ),
  list(
    name   = "CD45neg_MHCIIhi_vs_CD45neg_MHCIIlo",
    label  = "NOD CD45- MHCIIhi vs CD45- MHCIIlo",
    con    = c("condition", "CD45neg_MHCIIhi", "CD45neg_MHCIIlo"),
    groups = c("CD45neg_MHCIIlo", "CD45neg_MHCIIhi")
  )
)

# ── Color palette ──────────────────────────────────────────────────────────────
cond_colors <- c(
  "CD45+ MHCIIpos" = "#2166AC",
  "CD45- MHCIIhi"  = "#D6604D",
  "CD45- MHCIIlo"  = "#4DAC26"
)

# ── Helper: clean condition labels for plots ───────────────────────────────────
clean_label <- function(x) {
  x <- gsub("CD45pos_MHCIIpos", "CD45+ MHCIIpos", x)
  x <- gsub("CD45neg_MHCIIhi",  "CD45- MHCIIhi",  x)
  x <- gsub("CD45neg_MHCIIlo",  "CD45- MHCIIlo",  x)
  x
}

# ── Loop over contrasts ────────────────────────────────────────────────────────
for (ct in contrasts) {

  cat("\n──", ct$label, "──\n")

  # ── DESeq2 results ───────────────────────────────────────────────────────────
  res <- as.data.frame(results(dds, contrast=ct$con, alpha=0.05))
  res$ensembl <- rownames(res)
  res$symbol  <- to_sym(res$ensembl)
  res <- res %>% filter(!is.na(padj), !is.na(log2FoldChange))

  sig <- res %>% filter(padj < 0.05, abs(log2FoldChange) >= 1)
  cat("Significant DEGs:", nrow(sig), "\n")

  # Top 10 up, top 10 down by padj — deduplicate symbols (keep best padj per symbol)
  top_up   <- res %>% filter(log2FoldChange > 0) %>%
    arrange(padj) %>% distinct(symbol, .keep_all=TRUE) %>% slice_head(n=10)
  top_down <- res %>% filter(log2FoldChange < 0) %>%
    arrange(padj) %>% distinct(symbol, .keep_all=TRUE) %>% slice_head(n=10)
  top20    <- bind_rows(top_up, top_down)
  cat("Top 20 genes:", paste(top20$symbol, collapse=", "), "\n")

  # ── Volcano plot ─────────────────────────────────────────────────────────────
  res <- res %>% mutate(
    color_group = case_when(
      padj < 0.05 & log2FoldChange >=  1 ~ "Up",
      padj < 0.05 & log2FoldChange <= -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

  label_df <- res %>% filter(ensembl %in% top20$ensembl)

  p_vol <- ggplot(res, aes(x=log2FoldChange, y=-log10(padj), color=color_group)) +
    geom_point(size=0.7, alpha=0.5) +
    geom_point(data=label_df, size=1.8, alpha=1) +
    geom_text_repel(
      data         = label_df,
      aes(label    = symbol),
      size         = 3,
      fontface     = "italic",
      box.padding  = 0.4,
      max.overlaps = 20,
      segment.size = 0.3,
      color        = "black"
    ) +
    geom_vline(xintercept=c(-1, 1), linetype="dashed", color="grey50", linewidth=0.4) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey50", linewidth=0.4) +
    scale_color_manual(
      values = c(Up="firebrick", Down="steelblue", NS="grey75"),
      labels = c(
        Up   = paste0("Up (n=",   sum(res$color_group=="Up"),   ")"),
        Down = paste0("Down (n=", sum(res$color_group=="Down"), ")"),
        NS   = paste0("NS (n=",   sum(res$color_group=="NS"),   ")")
      )
    ) +
    labs(
      title    = ct$label,
      subtitle = "NOD cells only | padj<0.05 & |log2FC|>1",
      x        = expression(log[2]~Fold~Change),
      y        = expression(-log[10]~(padj)),
      color    = NULL
    ) +
    theme_bw(base_size=12) +
    theme(
      plot.title       = element_text(face="bold", size=11),
      plot.subtitle    = element_text(size=9, color="grey40"),
      legend.position  = "top",
      panel.grid.minor = element_blank()
    )

  vol_path <- file.path(out_dir, paste0("volcano_NOD_", ct$name, ".pdf"))
  ggsave(vol_path, p_vol, width=7, height=6)
  cat("Saved:", vol_path, "\n")

  # ── Violin plots ─────────────────────────────────────────────────────────────
  grp_cells <- meta$cell_id[meta$condition %in% ct$groups]
  expr_sub  <- expr[top20$ensembl, grp_cells, drop=FALSE]
  rownames(expr_sub) <- top20$symbol   # use gene symbols as row names
  meta_sub  <- meta[grp_cells, ]

  # Long format
  vln_df <- as.data.frame(t(expr_sub)) %>%
    tibble::rownames_to_column("cell_id") %>%
    left_join(meta_sub[, c("cell_id","condition")], by="cell_id") %>%
    pivot_longer(-c(cell_id, condition), names_to="gene", values_to="VST") %>%
    mutate(
      condition_clean = clean_label(condition),
      gene = factor(gene, levels=top20$symbol)
    )

  # Build one violin per gene
  gene_plots <- lapply(top20$symbol, function(g) {
    df_g      <- vln_df %>% filter(gene == g)
    gene_fc   <- round(top20$log2FoldChange[top20$symbol == g][1], 2)
    gene_pj   <- signif(top20$padj[top20$symbol == g][1], 2)
    direction <- ifelse(gene_fc > 0, "\u25b2", "\u25bc")  # ▲ ▼

    ggplot(df_g, aes(x=condition_clean, y=VST, fill=condition_clean)) +
      geom_violin(trim=TRUE, scale="width", alpha=0.8, linewidth=0.3) +
      geom_jitter(width=0.15, size=0.5, alpha=0.4, color="grey20") +
      scale_fill_manual(values=cond_colors) +
      labs(
        title    = paste0(g, "  ", direction),
        subtitle = paste0("log2FC=", gene_fc, "  padj=", gene_pj),
        x        = NULL,
        y        = "VST expression"
      ) +
      theme_bw(base_size=9) +
      theme(
        plot.title       = element_text(face="bold.italic", size=9),
        plot.subtitle    = element_text(size=7, color="grey40"),
        legend.position  = "none",
        axis.text.x      = element_text(angle=30, hjust=1, size=7),
        panel.grid.minor = element_blank()
      )
  })

  # Arrange 4 rows x 5 cols (20 panels)
  vln_panel <- wrap_plots(gene_plots, nrow=4, ncol=5) +
    plot_annotation(
      title    = paste0("Top 20 DEGs — ", ct$label),
      subtitle = "VST-normalized expression | Top 10 up + Top 10 down by padj",
      theme    = theme(
        plot.title    = element_text(face="bold", size=13),
        plot.subtitle = element_text(size=9, color="grey40")
      )
    )

  vln_path <- file.path(out_dir, paste0("violins_NOD_", ct$name, ".pdf"))
  ggsave(vln_path, vln_panel, width=22, height=14)
  cat("Saved:", vln_path, "\n")
}

cat("\nAll 6 files written to:", out_dir, "\n")
