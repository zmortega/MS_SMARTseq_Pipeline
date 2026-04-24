suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

# ==============================================================================
# per_strain_plots.R
#
# Runs fully independent per-strain DESeq2 analyses (one per plate), each
# normalized using only that strain's CD45pos_MHCIIpos (A1-A12) wells for
# size factor estimation. Produces 7 output files per strain:
#   - 3 volcano PDFs
#   - 3 violin PDFs (20 panels each: top 10 up + top 10 down by padj)
#   - 1 expression matrix CSV (all genes, mean VST per condition group)
#
# DESTRUCTIVE: Wipes results/05_dge/ before writing new outputs.
# Upstream folders (01-04, qc_summary) are not touched.
#
# To run on a different machine: update base_dir below.
# ==============================================================================

# -- Paths ---------------------------------------------------------------------
base_dir <- "/Users/zachortega/Desktop/MS_SMARTseq_Pipeline"
counts_f <- file.path(base_dir, "results/04_counts/counts_clean.txt")
meta_f   <- file.path(base_dir, "data/metadata.csv")
dge_dir  <- file.path(base_dir, "results/05_dge")

# -- Self-copy into scripts/ for version control -------------------------------
local({
  this_script <- normalizePath(
    grep("--file=", commandArgs(trailingOnly=FALSE), value=TRUE) |>
      sub("--file=", "", x=_),
    mustWork=FALSE
  )
  if (length(this_script) == 1 && nchar(this_script) > 0) {
    dest_dir <- file.path(base_dir, "scripts")
    dir.create(dest_dir, showWarnings=FALSE, recursive=TRUE)
    dest <- file.path(dest_dir, "per_strain_plots.R")
    if (file.copy(this_script, dest, overwrite=TRUE)) {
      cat("Script copied to:", dest, "\n")
    }
  }
})

# -- Wipe old 05_dge/ and rebuild ----------------------------------------------
if (dir.exists(dge_dir)) {
  unlink(dge_dir, recursive=TRUE)
  cat("Removed old results/05_dge/\n")
}
dir.create(dge_dir, recursive=TRUE)
cat("Created fresh results/05_dge/\n\n")

# -- Load data -----------------------------------------------------------------
counts <- read.table(counts_f, header=TRUE, row.names=1, sep="\t",
                     check.names=FALSE)
meta   <- read.csv(meta_f, stringsAsFactors=FALSE)
rownames(meta) <- meta$cell_id
cat("Total cells in count matrix:", ncol(counts), "\n")
cat("Total cells in metadata:    ", nrow(meta),   "\n\n")

# -- Ensembl -> gene symbol map ------------------------------------------------
sym_file <- "/tmp/ensembl_to_symbol.txt"
if (file.exists(sym_file)) {
  id2sym  <- read.table(sym_file, sep="\t", col.names=c("ensembl","symbol"),
                        stringsAsFactors=FALSE)
  sym_map <- setNames(id2sym$symbol, id2sym$ensembl)
  cat("Symbol map loaded:", length(sym_map), "entries\n\n")
} else {
  gtf_f <- file.path(base_dir,
                     "reference/gencode.vM33.primary_assembly.annotation.gtf")
  cat("Building symbol map from GTF (one-time, ~30s)...\n")
  gtf_lines  <- readLines(gtf_f)
  gene_lines <- gtf_lines[grepl('\tgene\t', gtf_lines)]
  ensembl    <- sub('.*gene_id "([^"]+)".*',   '\\1', gene_lines)
  symbol     <- sub('.*gene_name "([^"]+)".*', '\\1', gene_lines)
  sym_map    <- setNames(symbol, ensembl)
  write.table(data.frame(ensembl, symbol), sym_file, sep="\t",
              row.names=FALSE, col.names=FALSE, quote=FALSE)
  cat("Symbol map built:", length(sym_map), "genes\n\n")
}

to_sym <- function(ens) {
  sym <- sym_map[ens]
  ifelse(!is.na(sym) & sym != "", sym, ens)
}

# -- Shared settings -----------------------------------------------------------
strains <- c("NOD", "B6G7", "B6MHCIIGFP", "NODPDL1")

contrasts <- list(
  list(
    name   = "CD45pos_vs_MHCIIhi",
    label  = "CD45+ MHCIIpos vs CD45- MHCIIhi",
    con    = c("condition", "CD45neg_MHCIIhi", "CD45pos_MHCIIpos"),
    groups = c("CD45pos_MHCIIpos", "CD45neg_MHCIIhi")
  ),
  list(
    name   = "CD45pos_vs_MHCIIlo",
    label  = "CD45+ MHCIIpos vs CD45- MHCIIlo",
    con    = c("condition", "CD45neg_MHCIIlo", "CD45pos_MHCIIpos"),
    groups = c("CD45pos_MHCIIpos", "CD45neg_MHCIIlo")
  ),
  list(
    name   = "MHCIIhi_vs_MHCIIlo",
    label  = "CD45- MHCIIhi vs CD45- MHCIIlo",
    con    = c("condition", "CD45neg_MHCIIhi", "CD45neg_MHCIIlo"),
    groups = c("CD45neg_MHCIIlo", "CD45neg_MHCIIhi")
  )
)

cond_colors <- c(
  "CD45+ MHCIIpos" = "#2166AC",
  "CD45- MHCIIhi"  = "#D6604D",
  "CD45- MHCIIlo"  = "#4DAC26"
)

clean_label <- function(x) {
  x <- gsub("CD45pos_MHCIIpos", "CD45+ MHCIIpos", x)
  x <- gsub("CD45neg_MHCIIhi",  "CD45- MHCIIhi",  x)
  x <- gsub("CD45neg_MHCIIlo",  "CD45- MHCIIlo",  x)
  x
}

# ==============================================================================
# Main loop: one independent analysis per strain
# ==============================================================================
for (strain in strains) {

  cat("==============================================================\n")
  cat("Strain:", strain, "\n")
  cat("==============================================================\n")

  out_dir <- file.path(dge_dir, paste0(strain, "_plots"))
  dir.create(out_dir, recursive=TRUE)

  # -- Subset to this strain ---------------------------------------------------
  strain_cells <- meta$cell_id[meta$strain == strain]
  shared       <- intersect(colnames(counts), strain_cells)
  s_counts     <- counts[, shared]
  s_meta       <- meta[shared, ]
  cat("Cells:", ncol(s_counts), "\n")
  cat("Conditions:",
      paste(names(table(s_meta$condition)), collapse=" | "), "\n")
  cat("Counts per condition:\n")
  print(table(s_meta$condition))

  # -- Filter low-count genes --------------------------------------------------
  keep     <- rowSums(s_counts) >= 10
  s_counts <- s_counts[keep, ]
  cat("Genes passing filter:", nrow(s_counts), "\n")

  # -- Build DESeq2 object -----------------------------------------------------
  s_meta$condition <- relevel(factor(s_meta$condition),
                               ref="CD45pos_MHCIIpos")

  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(s_counts)),
    colData   = s_meta,
    design    = ~ condition
  )

  # -- Size factors from CD45pos_MHCIIpos (A1-A12) wells only -----------------
  ref_cells <- s_meta$cell_id[s_meta$condition == "CD45pos_MHCIIpos"]
  cat("Reference cells for size factor estimation:", length(ref_cells), "\n")

  # Use reference-based size factor estimation: compute geometric mean from
  # CD45pos_MHCIIpos (A1-A12) wells only, then apply median-of-ratios to
  # all cells in this strain so every cell is on the same scale.
  ref_counts   <- as.matrix(s_counts[, ref_cells])
  # Only use genes with non-zero counts in ALL reference cells for geo mean
  ref_nonzero  <- rowSums(ref_counts == 0) == 0
  ref_sub      <- ref_counts[ref_nonzero, , drop=FALSE]
  ref_geo_mean <- exp(rowMeans(log(ref_sub)))

  sf_all <- apply(as.matrix(s_counts), 2, function(cell_col) {
    ratios <- cell_col[ref_nonzero] / ref_geo_mean
    median(ratios[is.finite(ratios) & ratios > 0], na.rm=TRUE)
  })
  sf_all[sf_all <= 0 | !is.finite(sf_all)] <- 1
  sizeFactors(dds) <- sf_all
  cat("Size factors (summary):\n")
  print(summary(sf_all))

  # -- Run DESeq2 --------------------------------------------------------------
  # sfType not passed so DESeq2 uses the sizeFactors we assigned above
  dds <- DESeq(dds, test="Wald", fitType="parametric", quiet=TRUE)
  cat("DESeq2 complete.\n")

  # -- VST for expression matrix and violins -----------------------------------
  vsd  <- varianceStabilizingTransformation(dds, blind=FALSE)
  expr <- assay(vsd)

  # -- Expression matrix: mean VST per condition, all genes -------------------
  cat("Writing expression matrix...\n")
  cond_levels <- c("CD45pos_MHCIIpos", "CD45neg_MHCIIhi", "CD45neg_MHCIIlo")
  expr_mat <- data.frame(
    ensembl_id  = rownames(expr),
    gene_symbol = to_sym(rownames(expr)),
    stringsAsFactors = FALSE
  )
  for (cond in cond_levels) {
    cond_cells <- s_meta$cell_id[s_meta$condition == cond]
    cond_cells <- intersect(cond_cells, colnames(expr))
    if (length(cond_cells) == 0) {
      expr_mat[[paste0("mean_VST_", cond)]] <- NA
    } else if (length(cond_cells) == 1) {
      expr_mat[[paste0("mean_VST_", cond)]] <- expr[, cond_cells]
    } else {
      expr_mat[[paste0("mean_VST_", cond)]] <- rowMeans(expr[, cond_cells])
    }
  }
  colnames(expr_mat) <- c(
    "ensembl_id", "gene_symbol",
    "mean_VST_CD45pos_MHCIIpos",
    "mean_VST_CD45neg_MHCIIhi",
    "mean_VST_CD45neg_MHCIIlo"
  )

  # -- Upregulation rankings: combined score = log2FC * -log10(padj) -----------
  # Pull results for MHCIIhi vs CD45pos and MHCIIlo vs CD45pos
  rank_res_hi <- as.data.frame(results(dds,
    contrast=c("condition", "CD45neg_MHCIIhi", "CD45pos_MHCIIpos"),
    alpha=0.05))
  rank_res_hi$ensembl_id <- rownames(rank_res_hi)
  rank_res_hi$score_hi <- ifelse(
    rank_res_hi$log2FoldChange > 0 & !is.na(rank_res_hi$padj) & rank_res_hi$padj > 0,
    rank_res_hi$log2FoldChange * -log10(rank_res_hi$padj),
    NA
  )

  rank_res_lo <- as.data.frame(results(dds,
    contrast=c("condition", "CD45neg_MHCIIlo", "CD45pos_MHCIIpos"),
    alpha=0.05))
  rank_res_lo$ensembl_id <- rownames(rank_res_lo)
  rank_res_lo$score_lo <- ifelse(
    rank_res_lo$log2FoldChange > 0 & !is.na(rank_res_lo$padj) & rank_res_lo$padj > 0,
    rank_res_lo$log2FoldChange * -log10(rank_res_lo$padj),
    NA
  )

  # Merge scores into expr_mat and assign ranks (rank 1 = highest score)
  # Genes with log2FC <= 0 or NA padj get NA rank (not upregulated)
  expr_mat <- expr_mat %>%
    left_join(rank_res_hi[, c("ensembl_id", "score_hi")], by="ensembl_id") %>%
    left_join(rank_res_lo[, c("ensembl_id", "score_lo")], by="ensembl_id") %>%
    mutate(
      rank_upregulated_MHCIIhi = ifelse(
        !is.na(score_hi),
        rank(-score_hi, ties.method="min", na.last="keep"),
        NA
      ),
      rank_upregulated_MHCIIlo = ifelse(
        !is.na(score_lo),
        rank(-score_lo, ties.method="min", na.last="keep"),
        NA
      )
    ) %>%
    select(-score_hi, -score_lo)

  expr_mat <- expr_mat[order(expr_mat$ensembl_id), ]
  expr_path <- file.path(out_dir, paste0("expression_matrix_", strain, ".csv"))
  write.csv(expr_mat, expr_path, row.names=FALSE)
  cat("Saved:", expr_path, "\n")

  # -- Per-contrast plots ------------------------------------------------------
  for (ct in contrasts) {

    cat("\n  --", ct$label, "--\n")

    res <- as.data.frame(results(dds, contrast=ct$con, alpha=0.05))
    res$ensembl <- rownames(res)
    res$symbol  <- to_sym(res$ensembl)
    res <- res %>% filter(!is.na(padj), !is.na(log2FoldChange))

    sig <- res %>% filter(padj < 0.05, abs(log2FoldChange) >= 1)
    cat("  Significant DEGs:", nrow(sig), "\n")

    # Top 10 up + top 10 down, deduplicated by symbol
    top_up   <- res %>% filter(log2FoldChange > 0) %>%
      arrange(padj) %>% distinct(symbol, .keep_all=TRUE) %>% slice_head(n=10)
    top_down <- res %>% filter(log2FoldChange < 0) %>%
      arrange(padj) %>% distinct(symbol, .keep_all=TRUE) %>% slice_head(n=10)
    top20    <- bind_rows(top_up, top_down)
    cat("  Top 20 genes:", paste(top20$symbol, collapse=", "), "\n")

    # -- Volcano ---------------------------------------------------------------
    res <- res %>% mutate(
      color_group = case_when(
        padj < 0.05 & log2FoldChange >=  1 ~ "Up",
        padj < 0.05 & log2FoldChange <= -1 ~ "Down",
        TRUE ~ "NS"
      )
    )

    label_df <- res %>% filter(ensembl %in% top20$ensembl)

    p_vol <- ggplot(res, aes(x=log2FoldChange, y=-log10(padj),
                              color=color_group)) +
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
      geom_vline(xintercept=c(-1,1), linetype="dashed",
                 color="grey50", linewidth=0.4) +
      geom_hline(yintercept=-log10(0.05), linetype="dashed",
                 color="grey50", linewidth=0.4) +
      scale_color_manual(
        values = c(Up="firebrick", Down="steelblue", NS="grey75"),
        labels = c(
          Up   = paste0("Up (n=",   sum(res$color_group=="Up"),   ")"),
          Down = paste0("Down (n=", sum(res$color_group=="Down"), ")"),
          NS   = paste0("NS (n=",   sum(res$color_group=="NS"),   ")")
        )
      ) +
      labs(
        title    = paste0(strain, " — ", ct$label),
        subtitle = paste0(strain, " cells only | padj<0.05 & |log2FC|>1"),
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

    vol_path <- file.path(out_dir,
                           paste0("volcano_", strain, "_", ct$name, ".pdf"))
    ggsave(vol_path, p_vol, width=7, height=6)
    cat("  Saved:", vol_path, "\n")

    # -- Violins ---------------------------------------------------------------
    grp_cells <- s_meta$cell_id[s_meta$condition %in% ct$groups]
    expr_sub  <- expr[top20$ensembl, grp_cells, drop=FALSE]
    rownames(expr_sub) <- top20$symbol
    meta_sub  <- s_meta[grp_cells, ]

    vln_df <- as.data.frame(t(expr_sub)) %>%
      tibble::rownames_to_column("cell_id") %>%
      left_join(meta_sub[, c("cell_id","condition")], by="cell_id") %>%
      pivot_longer(-c(cell_id, condition), names_to="gene", values_to="VST") %>%
      mutate(
        condition_clean = clean_label(condition),
        gene = factor(gene, levels=top20$symbol)
      )

    gene_plots <- lapply(top20$symbol, function(g) {
      df_g      <- vln_df %>% filter(gene == g)
      gene_fc   <- round(top20$log2FoldChange[top20$symbol == g][1], 2)
      gene_pj   <- signif(top20$padj[top20$symbol == g][1], 2)
      direction <- ifelse(gene_fc > 0, "\u25b2", "\u25bc")

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

    vln_panel <- wrap_plots(gene_plots, nrow=4, ncol=5) +
      plot_annotation(
        title    = paste0(strain, " — Top 20 DEGs — ", ct$label),
        subtitle = "VST-normalized | Top 10 up + Top 10 down by padj | Normalized to CD45+ MHCIIpos (A1-A12)",
        theme    = theme(
          plot.title    = element_text(face="bold", size=13),
          plot.subtitle = element_text(size=9, color="grey40")
        )
      )

    vln_path <- file.path(out_dir,
                           paste0("violins_", strain, "_", ct$name, ".pdf"))
    ggsave(vln_path, vln_panel, width=22, height=14)
    cat("  Saved:", vln_path, "\n")
  }

  cat("\n", strain, "complete. 7 files written to:", out_dir, "\n\n")
}

cat("==============================================================\n")
cat("All strains complete. Output structure:\n")
for (strain in strains) {
  cat(" results/05_dge/", strain, "_plots/  (7 files)\n", sep="")
}
cat("==============================================================\n")
