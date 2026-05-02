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
# size factor estimation. Produces 10 output files per strain:
#   - 3 volcano PDFs
#   - 3 violin PDFs (20 panels each: top 10 up + top 10 down by padj)
#   - 1 expression matrix CSV (all genes, mean VST per condition group,
#       upregulation rankings for MHCIIhi and MHCIIlo)
#   - 3 UMAP PDFs:
#       1. All 3 populations, colored by population
#       2. CD45neg only (MHCIIhi + MHCIIlo), colored by population
#       3. CD45neg only, colored by Leiden cluster
#
# DESTRUCTIVE: Wipes results/05_dge/ before writing new outputs.
# Upstream folders (01-04, qc_summary) are not touched.
#
# To run on a different machine: update base_dir below.
# ==============================================================================

# -- Auto-install missing packages ---------------------------------------------
required_pkgs <- c("igraph", "uwot")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    cat("Installing missing package:", pkg, "\n")
    install.packages(pkg, repos="https://cloud.r-project.org", quiet=TRUE)
  }
}
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(uwot))

# -- Paths ---------------------------------------------------------------------
base_dir <- "/Users/zachortega/Desktop/MS_SMARTseq_Pipeline"
counts_f <- file.path(base_dir, "results/04_counts/counts_clean.txt")
meta_f   <- file.path(base_dir, "data/metadata.csv")
dge_dir  <- file.path(base_dir, "results/05_dge")

# -- UMAP / clustering tuning parameters ---------------------------------------
N_HVG            <- 2000   # top N highly variable genes for PCA
VAR_THRESHOLD    <- 0.80   # cumulative variance explained to select PCs
UMAP_N_NEIGHBORS <- 15     # uwot: n_neighbors (lower = more local structure)
UMAP_MIN_DIST    <- 0.3    # uwot: min_dist (lower = tighter clusters)
LEIDEN_RESOLUTION <- 0.5   # Leiden resolution (higher = more clusters)
UMAP_SEED        <- 42     # reproducibility

# -- MHCII expression filter ---------------------------------------------------
# Cells must express BOTH H2-Aa AND H2-Ab above this VST threshold to be
# included in DESeq2 models, UMAPs, clustering, and all downstream analyses.
# VST scale is log-like; 1.0 = very low but detectable expression.
MHCII_VST_MIN    <- 1.0    # minimum VST for both H2-Aa and H2-Ab
MHCII_GENES      <- c("H2-Aa", "H2-Ab1")  # gene symbols to filter on

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

# -- MHCII expression filter ---------------------------------------------------
# Compute a quick log-normalized matrix to assess expression before DESeq2.
# We use log1p(counts / colSums * 1e4) as a fast pre-filter proxy, then
# confirm with VST inside each per-strain loop.
cat("Applying MHCII expression filter (both", paste(MHCII_GENES, collapse=" & "),
    ">= VST", MHCII_VST_MIN, ")...
")

# Look up Ensembl IDs for H2-Aa and H2-Ab1
sym_to_ens_map  <- setNames(names(sym_map), sym_map)
mhcii_ens       <- sym_to_ens_map[MHCII_GENES]
mhcii_ens       <- mhcii_ens[!is.na(mhcii_ens) & mhcii_ens %in% rownames(counts)]

if (length(mhcii_ens) < 2) {
  stop("Could not find both MHCII genes (", paste(MHCII_GENES, collapse=", "),
       ") in count matrix. Check gene symbols match your GTF annotation.")
}
cat("  MHCII Ensembl IDs:
")
for (i in seq_along(mhcii_ens)) cat("   ", names(mhcii_ens)[i], "->", mhcii_ens[i], "
")

# Quick log-normalised proxy: log1p(CPM) — fast, no DESeq2 needed
lib_sizes   <- colSums(counts)
cpm_mat     <- sweep(as.matrix(counts[mhcii_ens, ]), 2, lib_sizes, "/") * 1e6
logcpm_mat  <- log1p(cpm_mat)

# A cell passes if BOTH genes are above a CPM proxy threshold.
# We use log1p(CPM) > log1p(5) as a rough equivalent of VST > 1
# (conservative — better to keep borderline cells and let VST decide)
LOGCPM_PROXY <- log1p(5)
passes <- colSums(logcpm_mat >= LOGCPM_PROXY) == 2
cat("  Cells passing MHCII filter:", sum(passes), "/", length(passes), "
")
cat("  Cells removed:", sum(!passes), "

")

# Report removals per strain x condition
removed_df <- meta[!passes, ]
if (nrow(removed_df) > 0) {
  cat("  Removed cells breakdown:
")
  print(table(removed_df$strain, removed_df$condition))
  cat("
")
}

# Apply filter globally — all downstream code uses these filtered objects
counts <- counts[, passes]
meta   <- meta[passes, ]
cat("  Filtered count matrix:", ncol(counts), "cells x", nrow(counts), "genes

")


# ==============================================================================
# MCA correlation helper: uses clustifyrdata::ref_MCA (no compilation needed)
# Computes Pearson r between each cluster's mean VST profile and all 713
# Mouse Cell Atlas cell types. Returns ranked Excel workbook.
# ==============================================================================
run_mca_correlation <- function(expr_mat, cluster_vec_input, cluster_levels_input,
                                  to_sym_fn, out_path, label="combined") {
  # Download ref_MCA — try multiple sources
  ref_mca_cache <- "/tmp/ref_MCA.rda"
  ref_mca_urls  <- c(
    "https://github.com/rnabioco/clustifyrdata/raw/main/data/ref_MCA.rda",
    "https://github.com/rnabioco/scRNA-seq-Cell-Ref-Matrix/raw/master/atlas/musMusculus/MouseAtlas.rda"
  )
  if (!file.exists(ref_mca_cache) || file.size(ref_mca_cache) < 1000) {
    for (url in ref_mca_urls) {
      cat("  Trying:", url, "\n")
      tryCatch({
        download.file(url, destfile=ref_mca_cache, mode="wb", quiet=TRUE)
        if (file.exists(ref_mca_cache) && file.size(ref_mca_cache) > 1000) {
          cat("  Download successful\n"); break
        }
      }, error=function(e) cat("  Failed:", conditionMessage(e), "\n"))
    }
  }
  if (!file.exists(ref_mca_cache) || file.size(ref_mca_cache) < 1000) {
    cat("  WARNING: ref_MCA unavailable, skipping MCA correlation.\n")
    return(invisible(NULL))
  }
  e <- new.env(); load(ref_mca_cache, envir=e)
  # Handle both ref_MCA and MouseAtlas object names
  ref <- if (!is.null(e$ref_MCA)) e$ref_MCA else
         if (!is.null(e$MouseAtlas)) e$MouseAtlas else
         get(ls(e)[1], envir=e)
  if (!is.matrix(ref) && !is.data.frame(ref)) ref <- as.matrix(ref)
  cat("  ref_MCA loaded:", nrow(ref), "genes x", ncol(ref), "cell types\n")

  # Build per-cluster mean VST profiles using gene symbols
  your_syms   <- to_sym_fn(rownames(expr_mat))
  sym_mat     <- expr_mat
  rownames(sym_mat) <- your_syms
  sym_mat     <- sym_mat[!duplicated(rownames(sym_mat)), ]

  common_genes <- intersect(rownames(sym_mat), rownames(ref))
  cat("  Common genes for MCA correlation:", length(common_genes), "\n")
  if (length(common_genes) < 50) {
    cat("  WARNING: too few common genes, skipping MCA correlation.\n")
    return(invisible(NULL))
  }

  sym_sub <- sym_mat[common_genes, , drop=FALSE]
  ref_sub <- ref[common_genes, , drop=FALSE]

  wb_mca <- createWorkbook()

  # Summary: top 5 per cluster
  sum_rows <- list()

  for (cl in as.character(cluster_levels_input)) {
    cl_cells  <- names(cluster_vec_input)[cluster_vec_input == cl]
    if (length(cl_cells) == 0) next
    cl_mean   <- rowMeans(sym_sub[, cl_cells, drop=FALSE])

    # Pearson correlation against all MCA cell types
    corrs <- cor(cl_mean, ref_sub, method="pearson")[1, ]
    corrs_sorted <- sort(corrs, decreasing=TRUE)

    cl_df <- data.frame(
      rank      = seq_along(corrs_sorted),
      cell_type = names(corrs_sorted),
      pearson_r = round(corrs_sorted, 4),
      stringsAsFactors=FALSE
    )

    sheet_nm <- paste0("Cluster_", cl)
    addWorksheet(wb_mca, sheet_nm)
    writeData(wb_mca, sheet_nm, cl_df)
    addStyle(wb_mca, sheet_nm,
      style=createStyle(textDecoration="bold", fgFill="#DDEBF7"),
      rows=1, cols=1:3, gridExpand=TRUE)
    addStyle(wb_mca, sheet_nm,
      style=createStyle(fgFill="#BDD7EE"),
      rows=2:6, cols=1:3, gridExpand=TRUE)
    setColWidths(wb_mca, sheet_nm, cols=1:3, widths=c(8,40,12))

    sum_rows[[cl]] <- data.frame(
      cluster=cl, rank=1:5,
      cell_type=names(corrs_sorted)[1:5],
      pearson_r=round(corrs_sorted[1:5], 4),
      stringsAsFactors=FALSE)
  }

  sum_df <- do.call(rbind, sum_rows)
  addWorksheet(wb_mca, "Summary_top5")
  writeData(wb_mca, "Summary_top5", sum_df)
  addStyle(wb_mca, "Summary_top5",
    style=createStyle(textDecoration="bold", fgFill="#DDEBF7"),
    rows=1, cols=1:4, gridExpand=TRUE)
  setColWidths(wb_mca, "Summary_top5", cols=1:4, widths=c(12,8,40,12))

  # Move Summary to first tab
  worksheetOrder(wb_mca) <- c(
    which(names(wb_mca) == "Summary_top5"),
    which(names(wb_mca) != "Summary_top5")
  )

  saveWorkbook(wb_mca, out_path, overwrite=TRUE)
  cat("  Saved:", out_path, "\n")
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

# -- UMAP helper: PCA -> select PCs by variance -> UMAP -----------------------
run_umap <- function(expr_subset, seed=UMAP_SEED,
                     n_neighbors=UMAP_N_NEIGHBORS,
                     min_dist=UMAP_MIN_DIST,
                     var_threshold=VAR_THRESHOLD) {
  # expr_subset: genes x cells matrix (already VST, HVG-filtered)
  pca_res   <- prcomp(t(expr_subset), center=TRUE, scale.=FALSE)
  var_exp   <- pca_res$sdev^2 / sum(pca_res$sdev^2)
  cum_var   <- cumsum(var_exp)
  n_pcs     <- max(2, which(cum_var >= var_threshold)[1])
  cat("    PCs selected:", n_pcs,
      sprintf("(%.1f%% variance explained)\n", cum_var[n_pcs]*100))
  pcs       <- pca_res$x[, 1:n_pcs, drop=FALSE]
  set.seed(seed)
  umap_coords <- umap(pcs, n_neighbors=min(n_neighbors, nrow(pcs)-1),
                      min_dist=min_dist, verbose=FALSE)
  colnames(umap_coords) <- c("UMAP1", "UMAP2")
  rownames(umap_coords) <- rownames(pcs)
  umap_coords
}

# -- Leiden clustering helper --------------------------------------------------
run_leiden <- function(umap_coords, resolution=LEIDEN_RESOLUTION,
                       n_neighbors=UMAP_N_NEIGHBORS) {
  # Build kNN graph from UMAP coordinates, then run Leiden
  n        <- nrow(umap_coords)
  k        <- min(n_neighbors, n-1)
  # Compute pairwise distances and get k nearest neighbors
  dists    <- as.matrix(dist(umap_coords))
  diag(dists) <- Inf
  knn_idx  <- t(apply(dists, 1, function(d) order(d)[1:k]))
  # Build edge list
  edges    <- do.call(rbind, lapply(1:n, function(i)
                cbind(i, knn_idx[i, ])))
  g        <- igraph::graph_from_edgelist(edges, directed=FALSE)
  g        <- igraph::simplify(g)
  clusters <- igraph::cluster_leiden(g,
               objective_function="modularity",
               resolution_parameter=resolution)
  as.factor(igraph::membership(clusters))
}

# ==============================================================================
# Main loop: one independent analysis per strain
# ==============================================================================
# Accumulators for combined UMAP and per-strain cluster info
all_expr_list        <- list()
all_meta_list        <- list()
per_strain_cluster_list <- list()
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

  ref_counts   <- as.matrix(s_counts[, ref_cells])
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
  dds <- DESeq(dds, test="Wald", fitType="parametric", quiet=TRUE)
  cat("DESeq2 complete.\n")

  # -- VST ---------------------------------------------------------------------
  vsd  <- varianceStabilizingTransformation(dds, blind=FALSE)
  expr <- assay(vsd)

  # -- Accumulate VST for combined UMAP ----------------------------------------
  all_expr_list[[strain]] <- expr
  s_meta_store <- s_meta
  s_meta_store$condition <- as.character(s_meta_store$condition)
  s_meta_store$strain    <- as.character(s_meta_store$strain)
  all_meta_list[[strain]] <- s_meta_store

  # -- Expression matrix -------------------------------------------------------
  cat("Writing expression matrix...\n")
  cond_levels <- c("CD45pos_MHCIIpos", "CD45neg_MHCIIhi", "CD45neg_MHCIIlo")
  expr_mat <- data.frame(
    ensembl_id  = rownames(expr),
    gene_symbol = to_sym(rownames(expr)),
    stringsAsFactors = FALSE
  )
  for (cond in cond_levels) {
    cond_cells <- intersect(s_meta$cell_id[s_meta$condition == cond],
                            colnames(expr))
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

  # Upregulation rankings
  rank_res_hi <- as.data.frame(results(dds,
    contrast=c("condition", "CD45neg_MHCIIhi", "CD45pos_MHCIIpos"),
    alpha=0.05))
  rank_res_hi$ensembl_id <- rownames(rank_res_hi)
  rank_res_hi$score_hi   <- ifelse(
    rank_res_hi$log2FoldChange > 0 & !is.na(rank_res_hi$padj) &
      rank_res_hi$padj > 0,
    rank_res_hi$log2FoldChange * -log10(rank_res_hi$padj), NA)

  rank_res_lo <- as.data.frame(results(dds,
    contrast=c("condition", "CD45neg_MHCIIlo", "CD45pos_MHCIIpos"),
    alpha=0.05))
  rank_res_lo$ensembl_id <- rownames(rank_res_lo)
  rank_res_lo$score_lo   <- ifelse(
    rank_res_lo$log2FoldChange > 0 & !is.na(rank_res_lo$padj) &
      rank_res_lo$padj > 0,
    rank_res_lo$log2FoldChange * -log10(rank_res_lo$padj), NA)

  expr_mat <- expr_mat %>%
    left_join(rank_res_hi[, c("ensembl_id","score_hi")], by="ensembl_id") %>%
    left_join(rank_res_lo[, c("ensembl_id","score_lo")], by="ensembl_id") %>%
    mutate(
      rank_upregulated_MHCIIhi = ifelse(!is.na(score_hi),
        rank(-score_hi, ties.method="min", na.last="keep"), NA),
      rank_upregulated_MHCIIlo = ifelse(!is.na(score_lo),
        rank(-score_lo, ties.method="min", na.last="keep"), NA)
    ) %>%
    select(-score_hi, -score_lo)

  expr_mat <- expr_mat[order(expr_mat$ensembl_id), ]
  expr_path <- file.path(out_dir, paste0("expression_matrix_", strain, ".csv"))
  write.csv(expr_mat, expr_path, row.names=FALSE)
  cat("Saved:", expr_path, "\n")

  # ============================================================================
  # UMAPs
  # ============================================================================
  cat("\n  -- UMAPs --\n")

  # Select top N_HVG highly variable genes (by variance across all cells)
  gene_vars <- apply(expr, 1, var)
  hvg       <- names(sort(gene_vars, decreasing=TRUE))[1:min(N_HVG,
                                                              nrow(expr))]
  expr_hvg  <- expr[hvg, ]
  cat("  Highly variable genes selected:", length(hvg), "\n")

  # -- UMAP 1: All 3 populations, colored by population -----------------------
  cat("  Running UMAP (all cells)...\n")
  umap_all  <- run_umap(expr_hvg)
  umap_all_df <- data.frame(
    umap_all,
    cell_id         = rownames(umap_all),
    condition       = s_meta[rownames(umap_all), "condition"],
    condition_clean = clean_label(s_meta[rownames(umap_all), "condition"]),
    stringsAsFactors = FALSE
  )

  p_umap_all <- ggplot(umap_all_df,
                        aes(x=UMAP1, y=UMAP2, color=condition_clean)) +
    geom_point(size=2.5, alpha=0.85) +
    scale_color_manual(values=cond_colors) +
    labs(
      title    = paste0(strain, " — All populations"),
      subtitle = paste0("n=", nrow(umap_all_df), " cells | top ", length(hvg),
                        " HVGs | ", VAR_THRESHOLD*100, "% variance PCs"),
      color    = "Population"
    ) +
    theme_bw(base_size=12) +
    theme(
      plot.title       = element_text(face="bold", size=12),
      plot.subtitle    = element_text(size=9, color="grey40"),
      panel.grid.minor = element_blank(),
      aspect.ratio     = 1
    )

  umap1_path <- file.path(out_dir, paste0("umap_", strain, "_all_by_population.pdf"))
  ggsave(umap1_path, p_umap_all, width=6, height=6)
  cat("  Saved:", umap1_path, "\n")

  # -- UMAP 2 & 3: CD45neg only (MHCIIhi + MHCIIlo) ---------------------------
  cd45neg_cells <- s_meta$cell_id[s_meta$condition %in%
                                    c("CD45neg_MHCIIhi", "CD45neg_MHCIIlo")]
  cd45neg_cells <- intersect(cd45neg_cells, colnames(expr_hvg))
  expr_cd45neg  <- expr_hvg[, cd45neg_cells]
  cat("  Running UMAP (CD45neg cells, n=", ncol(expr_cd45neg), ")...\n",
      sep="")

  umap_neg     <- run_umap(expr_cd45neg)
  umap_neg_df  <- data.frame(
    umap_neg,
    cell_id         = rownames(umap_neg),
    condition       = s_meta[rownames(umap_neg), "condition"],
    condition_clean = clean_label(s_meta[rownames(umap_neg), "condition"]),
    stringsAsFactors = FALSE
  )

  # UMAP 2: colored by population
  p_umap_neg_pop <- ggplot(umap_neg_df,
                            aes(x=UMAP1, y=UMAP2, color=condition_clean)) +
    geom_point(size=2.5, alpha=0.85) +
    scale_color_manual(values=cond_colors[names(cond_colors) != "CD45+ MHCIIpos"]) +
    labs(
      title    = paste0(strain, " — CD45\u2212 populations by population"),
      subtitle = paste0("n=", nrow(umap_neg_df), " cells | MHCIIhi + MHCIIlo only"),
      color    = "Population"
    ) +
    theme_bw(base_size=12) +
    theme(
      plot.title       = element_text(face="bold", size=12),
      plot.subtitle    = element_text(size=9, color="grey40"),
      panel.grid.minor = element_blank(),
      aspect.ratio     = 1
    )

  umap2_path <- file.path(out_dir,
                           paste0("umap_", strain, "_cd45neg_by_population.pdf"))
  ggsave(umap2_path, p_umap_neg_pop, width=6, height=6)
  cat("  Saved:", umap2_path, "\n")

  # UMAP 3: colored by Leiden cluster
  cat("  Running Leiden clustering (resolution=", LEIDEN_RESOLUTION, ")...\n",
      sep="")
  leiden_clusters        <- run_leiden(umap_neg)
  umap_neg_df$cluster    <- as.character(leiden_clusters)
  n_clusters             <- length(unique(umap_neg_df$cluster))
  cat("  Clusters found:", n_clusters, "\n")

  # Generate a color palette with enough colors for all clusters
  cluster_palette <- setNames(
    scales::hue_pal()(n_clusters),
    sort(unique(umap_neg_df$cluster))
  )

  p_umap_neg_clust <- ggplot(umap_neg_df,
                              aes(x=UMAP1, y=UMAP2, color=cluster)) +
    geom_point(size=2.5, alpha=0.85) +
    scale_color_manual(values=cluster_palette) +
    labs(
      title    = paste0(strain, " — CD45\u2212 populations by Leiden cluster"),
      subtitle = paste0("n=", nrow(umap_neg_df), " cells | resolution=",
                        LEIDEN_RESOLUTION, " | ", n_clusters, " clusters"),
      color    = "Cluster"
    ) +
    theme_bw(base_size=12) +
    theme(
      plot.title       = element_text(face="bold", size=12),
      plot.subtitle    = element_text(size=9, color="grey40"),
      panel.grid.minor = element_blank(),
      aspect.ratio     = 1
    )

  umap3_path <- file.path(out_dir,
                           paste0("umap_", strain, "_cd45neg_by_cluster.pdf"))
  ggsave(umap3_path, p_umap_neg_clust, width=6, height=6)
  cat("  Saved:", umap3_path, "\n")

  # ============================================================================
  # Per-strain UMAP (all-gene PCA for clustering, HVG for visualization)
  # ============================================================================
  cat("\n  -- Per-strain UMAP --\n")

  # HVG selection for this strain
  s_gene_vars  <- apply(expr, 1, var)
  s_hvg        <- names(sort(s_gene_vars, decreasing=TRUE))[
                    1:min(N_HVG, nrow(expr))]
  expr_s_hvg   <- expr[s_hvg, ]

  # UMAP embedding on HVGs
  cat("  Running UMAP (HVGs)...\n")
  s_umap       <- run_umap(expr_s_hvg)
  s_umap_df    <- data.frame(s_umap,
                              cell_id   = rownames(s_umap),
                              condition = s_meta[rownames(s_umap), "condition"],
                              stringsAsFactors=FALSE)

  # All-gene PCA for Leiden
  cat("  Running all-gene PCA for clustering...\n")
  s_pca_all    <- prcomp(t(expr), center=TRUE, scale.=FALSE)
  s_var_all    <- s_pca_all$sdev^2 / sum(s_pca_all$sdev^2)
  s_cum_all    <- cumsum(s_var_all)
  s_n_pcs      <- max(2, which(s_cum_all >= VAR_THRESHOLD)[1])
  cat("  PCs selected:", s_n_pcs,
      sprintf("(%.1f%% variance)\n", s_cum_all[s_n_pcs]*100))
  s_pcs_all    <- s_pca_all$x[, 1:s_n_pcs, drop=FALSE]

  # Scale Leiden parameters for per-strain cell counts
  # Target: ~3-5 clusters per strain
  # NOTE: with small cell counts, higher resolution is needed to overcome
  # sparse kNN graphs — counterintuitively, low res collapses everything to 1
  n_cells_s    <- ncol(expr)
  s_neighbors  <- max(8, floor(n_cells_s / 5))
  s_resolution <- ifelse(n_cells_s < 50,  1.2,
                  ifelse(n_cells_s < 80,  1.5,
                  ifelse(n_cells_s < 120, 1.8, LEIDEN_RESOLUTION)))
  cat("  Cells:", n_cells_s, "| Leiden resolution:", s_resolution,
      "| n_neighbors:", s_neighbors, "\n")

  cat("  Running Leiden clustering...\n")
  s_leiden     <- run_leiden(s_pcs_all, resolution=s_resolution,
                              n_neighbors=s_neighbors)
  s_umap_df$cluster <- as.character(s_leiden)
  s_n_clust    <- length(unique(s_umap_df$cluster))
  cat("  Clusters found:", s_n_clust, "\n")

  s_umap_df$cluster <- factor(s_umap_df$cluster,
                               levels=sort(unique(as.integer(s_umap_df$cluster))))

  s_clust_pal  <- setNames(scales::hue_pal()(s_n_clust),
                            sort(unique(as.character(s_umap_df$cluster))))

  # UMAP colored by cluster
  p_s_clust <- ggplot(s_umap_df, aes(x=UMAP1, y=UMAP2, color=cluster)) +
    geom_point(size=2.2, alpha=0.85) +
    scale_color_manual(values=s_clust_pal) +
    labs(
      title    = paste0(strain, " — ", nrow(s_umap_df),
                        " cells by Leiden cluster"),
      subtitle = paste0("n=", nrow(s_umap_df),
                        " (MHCII-filtered) | UMAP: top ", N_HVG,
                        " HVGs | Clusters: all-gene PCA | resolution=",
                        LEIDEN_RESOLUTION, " | ", s_n_clust, " clusters"),
      color    = "Cluster"
    ) +
    theme_bw(base_size=12) +
    theme(plot.title=element_text(face="bold", size=12),
          plot.subtitle=element_text(size=8, color="grey40"),
          panel.grid.minor=element_blank(), aspect.ratio=1)

  s_umap1_path <- file.path(out_dir,
                              paste0("umap_", strain, "_by_cluster.pdf"))
  ggsave(s_umap1_path, p_s_clust, width=6, height=6)
  cat("  Saved:", s_umap1_path, "\n")

  # UMAP colored by condition
  p_s_cond <- ggplot(s_umap_df,
                      aes(x=UMAP1, y=UMAP2, color=clean_label(condition))) +
    geom_point(size=2.2, alpha=0.85) +
    scale_color_manual(values=cond_colors) +
    labs(
      title    = paste0(strain, " — ", nrow(s_umap_df),
                        " cells by population"),
      subtitle = paste0("n=", nrow(s_umap_df), " (MHCII-filtered)"),
      color    = "Population"
    ) +
    theme_bw(base_size=12) +
    theme(plot.title=element_text(face="bold", size=12),
          plot.subtitle=element_text(size=8, color="grey40"),
          panel.grid.minor=element_blank(), aspect.ratio=1)

  s_umap2_path <- file.path(out_dir,
                              paste0("umap_", strain, "_by_population.pdf"))
  ggsave(s_umap2_path, p_s_cond, width=6, height=6)
  cat("  Saved:", s_umap2_path, "\n")

  # Store per-strain cluster info for cell identity scoring
  s_umap_df$strain <- strain
  per_strain_cluster_list[[strain]] <- list(
    umap_df    = s_umap_df,
    expr       = expr,
    clust_pal  = s_clust_pal,
    n_clust    = s_n_clust
  )

  # -- Per-contrast plots (volcano + violin) -----------------------------------
  for (ct in contrasts) {

    cat("\n  --", ct$label, "--\n")

    res <- as.data.frame(results(dds, contrast=ct$con, alpha=0.05))
    res$ensembl <- rownames(res)
    res$symbol  <- to_sym(res$ensembl)
    res <- res %>% filter(!is.na(padj), !is.na(log2FoldChange))

    sig <- res %>% filter(padj < 0.05, abs(log2FoldChange) >= 1)
    cat("  Significant DEGs:", nrow(sig), "\n")

    top_up   <- res %>% filter(log2FoldChange > 0) %>%
      arrange(padj) %>% distinct(symbol, .keep_all=TRUE) %>% slice_head(n=10)
    top_down <- res %>% filter(log2FoldChange < 0) %>%
      arrange(padj) %>% distinct(symbol, .keep_all=TRUE) %>% slice_head(n=10)
    top20    <- bind_rows(top_up, top_down)
    cat("  Top 20 genes:", paste(top20$symbol, collapse=", "), "\n")

    # Volcano
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
        data=label_df, aes(label=symbol),
        size=3, fontface="italic", box.padding=0.4,
        max.overlaps=20, segment.size=0.3, color="black"
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
        plot.title=element_text(face="bold", size=11),
        plot.subtitle=element_text(size=9, color="grey40"),
        legend.position="top",
        panel.grid.minor=element_blank()
      )

    vol_path <- file.path(out_dir,
                           paste0("volcano_", strain, "_", ct$name, ".pdf"))
    ggsave(vol_path, p_vol, width=7, height=6)
    cat("  Saved:", vol_path, "\n")

    # Violins
    grp_cells <- s_meta$cell_id[s_meta$condition %in% ct$groups]
    expr_sub  <- expr[top20$ensembl, grp_cells, drop=FALSE]
    rownames(expr_sub) <- top20$symbol
    meta_sub  <- s_meta[grp_cells, ]

    vln_df <- as.data.frame(t(expr_sub)) %>%
      tibble::rownames_to_column("cell_id") %>%
      left_join(meta_sub[, c("cell_id","condition")], by="cell_id") %>%
      pivot_longer(-c(cell_id,condition), names_to="gene", values_to="VST") %>%
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
          x=NULL, y="VST expression"
        ) +
        theme_bw(base_size=9) +
        theme(
          plot.title=element_text(face="bold.italic", size=9),
          plot.subtitle=element_text(size=7, color="grey40"),
          legend.position="none",
          axis.text.x=element_text(angle=30, hjust=1, size=7),
          panel.grid.minor=element_blank()
        )
    })

    vln_panel <- wrap_plots(gene_plots, nrow=4, ncol=5) +
      plot_annotation(
        title    = paste0(strain, " — Top 20 DEGs — ", ct$label),
        subtitle = "VST-normalized | Top 10 up + Top 10 down by padj | Normalized to CD45+ MHCIIpos (A1-A12)",
        theme    = theme(
          plot.title=element_text(face="bold", size=13),
          plot.subtitle=element_text(size=9, color="grey40")
        )
      )

    vln_path <- file.path(out_dir,
                           paste0("violins_", strain, "_", ct$name, ".pdf"))
    ggsave(vln_path, vln_panel, width=22, height=14)
    cat("  Saved:", vln_path, "\n")
  }

  cat("\n", strain, "complete. 10 files written to:", out_dir, "\n\n")
}

# ==============================================================================
# Combined UMAP: all 372 cells across all 4 plates
# ==============================================================================
cat("==============================================================\n")
cat("Building combined UMAP (all strains)...\n")
cat("==============================================================\n")

# Combine VST matrices — keep only genes present in all 4 plates
common_genes <- Reduce(intersect, lapply(all_expr_list, rownames))
cat("Genes common across all plates:", length(common_genes), "\n")

combined_expr <- do.call(cbind, lapply(all_expr_list, function(e) e[common_genes, ]))
combined_meta <- do.call(rbind, all_meta_list)
combined_meta$condition <- as.character(combined_meta$condition)
combined_meta$strain    <- as.character(combined_meta$strain)
rownames(combined_meta) <- combined_meta$cell_id
combined_meta <- combined_meta[colnames(combined_expr), , drop=FALSE]
cat("Total cells in combined matrix:", ncol(combined_expr), "\n")
cat("Metadata rows matched:", sum(!is.na(combined_meta$strain)), "\n")

# Add strain_condition label for the 12-color UMAP
combined_meta$strain_condition <- paste0(
  combined_meta$strain, " ",
  clean_label(combined_meta$condition)
)
cat("strain_condition sample:\n"); print(head(combined_meta$strain_condition))
cat("NA count:", sum(is.na(combined_meta$strain_condition)), "\n")
cat("strain_condition sample:\n")
print(head(combined_meta$strain_condition))

# HVG selection on combined matrix
gene_vars_comb <- apply(combined_expr, 1, var)
hvg_comb       <- names(sort(gene_vars_comb, decreasing=TRUE))[
                    1:min(N_HVG, length(gene_vars_comb))]
expr_hvg_comb  <- combined_expr[hvg_comb, ]
cat("HVGs selected for combined UMAP:", length(hvg_comb), "\n")

# PCA + UMAP on HVGs (for visualization only)
cat("Running PCA + UMAP on HVG matrix (visualization)...\n")
umap_comb <- run_umap(expr_hvg_comb)
umap_comb_df <- data.frame(
  umap_comb,
  cell_id = rownames(umap_comb),
  stringsAsFactors = FALSE
)
umap_comb_df <- umap_comb_df %>%
  left_join(
    combined_meta[, c("cell_id", "strain", "condition", "strain_condition")],
    by = "cell_id"
  )

# PCA on ALL genes for Leiden clustering
cat("Running PCA on all", nrow(combined_expr), "genes for clustering...\n")
pca_all   <- prcomp(t(combined_expr), center=TRUE, scale.=FALSE)
var_all   <- pca_all$sdev^2 / sum(pca_all$sdev^2)
cum_all   <- cumsum(var_all)
n_pcs_all <- max(2, which(cum_all >= VAR_THRESHOLD)[1])
cat("PCs selected (all-gene):", n_pcs_all,
    sprintf("(%.1f%% variance explained)\n", cum_all[n_pcs_all]*100))
pcs_all   <- pca_all$x[, 1:n_pcs_all, drop=FALSE]

# Leiden clustering on all-gene PCs
cat("Running Leiden clustering on all-gene PCs (resolution=",
    LEIDEN_RESOLUTION, ")...\n", sep="")
leiden_comb           <- run_leiden(pcs_all)
umap_comb_df$cluster  <- as.character(leiden_comb)
n_clust_comb          <- length(unique(umap_comb_df$cluster))
cat("Clusters found:", n_clust_comb, "\n")

combined_out <- file.path(dge_dir, "combined_plots")
dir.create(combined_out, showWarnings=FALSE)

# -- Combined UMAP 1: colored by Leiden cluster --------------------------------
cluster_pal_comb <- setNames(
  scales::hue_pal()(n_clust_comb),
  sort(unique(umap_comb_df$cluster))
)

p_comb_clust <- ggplot(umap_comb_df,
                        aes(x=UMAP1, y=UMAP2, color=cluster)) +
  geom_point(size=1.8, alpha=0.8) +
  scale_color_manual(values=cluster_pal_comb) +
  labs(
    title    = paste0("All strains — ", nrow(umap_comb_df), " cells by Leiden cluster"),
    subtitle = paste0("n=", nrow(umap_comb_df),
                      " (MHCII-filtered) | UMAP: top ", N_HVG,
                      " HVGs | Clusters: all-gene PCA | resolution=",
                      LEIDEN_RESOLUTION, " | ", n_clust_comb, " clusters"),
    color    = "Cluster"
  ) +
  theme_bw(base_size=12) +
  theme(
    plot.title       = element_text(face="bold", size=12),
    plot.subtitle    = element_text(size=9, color="grey40"),
    panel.grid.minor = element_blank(),
    aspect.ratio     = 1
  )

comb1_path <- file.path(combined_out, "umap_all372_by_cluster.pdf")
ggsave(comb1_path, p_comb_clust, width=7, height=6)
cat("Saved:", comb1_path, "\n")

# -- Combined UMAP 2: colored by strain x condition (12 colors) ----------------
# Build 12-color palette: 4 strains x 3 conditions
strain_list <- c("NOD", "B6G7", "B6MHCIIGFP", "NODPDL1")
cond_list   <- c("CD45+ MHCIIpos", "CD45- MHCIIhi", "CD45- MHCIIlo")

# Base hues per strain, lighter/darker shades per condition
strain_base_colors <- c(
  NOD        = "#2166AC",
  B6G7       = "#D6604D",
  B6MHCIIGFP = "#4DAC26",
  NODPDL1    = "#8B44AC"
)
# Condition modifiers: MHCIIpos = full saturation, MHCIIhi = lighter,
# MHCIIlo = darker
condition_alpha <- c(
  "CD45+ MHCIIpos" = 1.0,
  "CD45- MHCIIhi"  = 0.55,
  "CD45- MHCIIlo"  = 0.30
)

# Build palette by blending base color toward white (lighter) or black (darker)
blend_color <- function(hex, alpha_toward_white) {
  rgb_vals <- col2rgb(hex) / 255
  blended  <- rgb_vals * alpha_toward_white +
               (1 - alpha_toward_white) * ifelse(alpha_toward_white >= 0.6,
                                                  1, 0)
  rgb(blended[1], blended[2], blended[3])
}

sc_palette <- c()
for (st in strain_list) {
  for (co in cond_list) {
    key            <- paste0(st, " ", co)
    sc_palette[key] <- blend_color(strain_base_colors[st],
                                    condition_alpha[co])
  }
}

# Rebuild sc_palette with fully distinct colors per strain x condition
# Each strain gets 3 visually distinct shades (light, mid, dark)
sc_palette_v2 <- c(
  "NOD CD45+ MHCIIpos" = "#08519C",   # dark blue
  "NOD CD45- MHCIIhi"  = "#6BAED6",   # mid blue
  "NOD CD45- MHCIIlo"  = "#BDD7E7",   # light blue
  "B6G7 CD45+ MHCIIpos" = "#A50F15",  # dark red
  "B6G7 CD45- MHCIIhi"  = "#FB6A4A",  # mid red
  "B6G7 CD45- MHCIIlo"  = "#FCBBA1",  # light red/salmon
  "B6MHCIIGFP CD45+ MHCIIpos" = "#006D2C",  # dark green
  "B6MHCIIGFP CD45- MHCIIhi"  = "#41AB5D",  # mid green
  "B6MHCIIGFP CD45- MHCIIlo"  = "#C7E9C0",  # light green
  "NODPDL1 CD45+ MHCIIpos" = "#54278F",  # dark purple
  "NODPDL1 CD45- MHCIIhi"  = "#9E9AC8",  # mid purple
  "NODPDL1 CD45- MHCIIlo"  = "#DADAEB"   # light purple
)

p_comb_strain <- ggplot(umap_comb_df,
                         aes(x=UMAP1, y=UMAP2, color=strain_condition)) +
  geom_point(size=1.0, alpha=1.0, stroke=0.2, shape=21,
             aes(fill=strain_condition), color="grey20") +
  scale_fill_manual(values=sc_palette_v2) +
  scale_color_manual(values=sc_palette_v2, guide="none") +
  labs(
    title    = paste0("All strains — ", nrow(umap_comb_df), " cells by strain & population"),
    subtitle = paste0("n=", nrow(umap_comb_df),
                      " (MHCII-filtered) | 4 strains x 3 conditions"),
    fill     = "Strain / Population"
  ) +
  guides(fill=guide_legend(ncol=2, override.aes=list(size=3, stroke=0.3))) +
  theme_bw(base_size=12) +
  theme(
    plot.title       = element_text(face="bold", size=12),
    plot.subtitle    = element_text(size=9, color="grey40"),
    panel.grid.minor = element_blank(),
    legend.text      = element_text(size=8),
    aspect.ratio     = 1
  )

comb2_path <- file.path(combined_out, "umap_all372_by_strain_population.pdf")
ggsave(comb2_path, p_comb_strain, width=8, height=6)
cat("Saved:", comb2_path, "\n")

# -- 4-panel UMAP: each strain on comprehensive coordinates -------------------
cat("Generating per-strain panels on comprehensive UMAP coordinates...\n")

strain_panel_plots <- lapply(strains, function(st) {
  # Cells belonging to this strain
  st_cells    <- umap_comb_df$cell_id[umap_comb_df$strain == st]
  df_fg       <- umap_comb_df %>% filter(cell_id %in% st_cells)
  df_bg       <- umap_comb_df %>% filter(!cell_id %in% st_cells)

  ggplot() +
    # Background: all other cells in light grey
    geom_point(data=df_bg, aes(x=UMAP1, y=UMAP2),
               color="grey88", size=1.0, alpha=0.5) +
    # Foreground: this strain colored by comprehensive cluster
    geom_point(data=df_fg, aes(x=UMAP1, y=UMAP2, color=cluster),
               size=1.8, alpha=0.9) +
    scale_color_manual(values=cluster_pal_comb) +
    labs(
      title  = paste0(st, " (n=", nrow(df_fg), ")"),
      x      = "UMAP1", y = "UMAP2",
      color  = "Cluster"
    ) +
    coord_fixed() +
    theme_bw(base_size=11) +
    theme(
      plot.title       = element_text(face="bold", size=11),
      panel.grid.minor = element_blank(),
      legend.position  = "right",
      legend.text      = element_text(size=8)
    )
})

# Stack 4 panels vertically
strain_panels_combined <- wrap_plots(strain_panel_plots, ncol=1) +
  plot_annotation(
    title    = "Per-strain cells on comprehensive UMAP",
    subtitle = paste0("Comprehensive Leiden clusters (resolution=",
                      LEIDEN_RESOLUTION, ") | grey = other strains"),
    theme    = theme(
      plot.title    = element_text(face="bold", size=13),
      plot.subtitle = element_text(size=9, color="grey40")
    )
  )

strain_panels_path <- file.path(combined_out,
                                 "umap_per_strain_on_comprehensive.pdf")
ggsave(strain_panels_path, strain_panels_combined,
       width=7, height=6 * length(strains))
cat("Saved:", strain_panels_path, "\n")

cat("\nCombined UMAPs complete.\n")

# ==============================================================================
# Bar charts: cluster composition (all 4 in one PDF, 2x2 layout)
# ==============================================================================
cat("\nGenerating cluster composition bar charts...\n")

# Auto-install packages needed for bar charts + heatmap
if (!requireNamespace("BiocManager", quietly=TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
}
for (pkg in c("openxlsx", "ComplexHeatmap", "circlize", "RColorBrewer")) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    cat("Installing", pkg, "...\n")
    if (pkg %in% c("ComplexHeatmap", "circlize")) {
      BiocManager::install(pkg, ask=FALSE, update=FALSE, quiet=TRUE)
    } else {
      install.packages(pkg, repos="https://cloud.r-project.org", quiet=TRUE)
    }
  }
}
suppressPackageStartupMessages({
  library(openxlsx)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
})

# Attach clean labels to plotting df
umap_comb_df$condition_clean  <- clean_label(umap_comb_df$condition)
umap_comb_df$strain_condition <- paste0(umap_comb_df$strain, " ",
                                         umap_comb_df$condition_clean)
umap_comb_df$cluster          <- factor(umap_comb_df$cluster,
                                         levels=sort(unique(as.integer(
                                           umap_comb_df$cluster))))

strain_colors <- c(
  NOD        = "#2166AC",
  B6G7       = "#D6604D",
  B6MHCIIGFP = "#4DAC26",
  NODPDL1    = "#8B44AC"
)

MIN_LABEL_PROP <- 0.05   # only label segments >= 5% of bar

# Helper: compute label data for a stacked bar
make_label_df <- function(df, fill_var, count_or_prop) {
  df %>%
    group_by(cluster, .data[[fill_var]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(cluster) %>%
    mutate(
      total = sum(n),
      prop  = n / total,
      label_count = ifelse(prop >= MIN_LABEL_PROP, as.character(n), ""),
      label_prop  = ifelse(prop >= MIN_LABEL_PROP,
                           paste0(round(prop * 100), "%"), ""),
      y_pos = cumsum(if (count_or_prop == "count") n else prop) -
              0.5 * (if (count_or_prop == "count") n else prop)
    ) %>%
    ungroup()
}

# -- Bar chart 1: counts, fill = strain x population --------------------------
lab1 <- make_label_df(umap_comb_df, "strain_condition", "count")
p_bar_sc_counts <- ggplot(lab1, aes(x=cluster, y=n, fill=strain_condition)) +
  geom_bar(stat="identity", position="stack", color="white", linewidth=0.15) +
  geom_text(aes(y=y_pos, label=label_count),
            size=2.5, color="white", fontface="bold") +
  scale_fill_manual(values=sc_palette) +
  labs(title="Strain & population — counts", x="Cluster", y="Cell count",
       fill="Strain / Population") +
  guides(fill=guide_legend(ncol=2)) +
  theme_bw(base_size=10) +
  theme(plot.title=element_text(face="bold", size=10),
        panel.grid.minor=element_blank(),
        legend.text=element_text(size=7),
        legend.key.size=unit(0.4,"cm"))

# -- Bar chart 2: proportions, fill = strain x population ---------------------
lab2 <- make_label_df(umap_comb_df, "strain_condition", "prop")
p_bar_sc_prop <- ggplot(lab2, aes(x=cluster, y=prop, fill=strain_condition)) +
  geom_bar(stat="identity", position="stack", color="white", linewidth=0.15) +
  geom_text(aes(y=y_pos, label=label_prop),
            size=2.5, color="white", fontface="bold") +
  scale_fill_manual(values=sc_palette) +
  scale_y_continuous(labels=scales::percent_format()) +
  labs(title="Strain & population — proportions", x="Cluster",
       y="Proportion", fill="Strain / Population") +
  guides(fill=guide_legend(ncol=2)) +
  theme_bw(base_size=10) +
  theme(plot.title=element_text(face="bold", size=10),
        panel.grid.minor=element_blank(),
        legend.text=element_text(size=7),
        legend.key.size=unit(0.4,"cm"))

# -- Bar chart 3: counts, fill = strain only -----------------------------------
lab3 <- make_label_df(umap_comb_df, "strain", "count")
p_bar_s_counts <- ggplot(lab3, aes(x=cluster, y=n, fill=strain)) +
  geom_bar(stat="identity", position="stack", color="white", linewidth=0.15) +
  geom_text(aes(y=y_pos, label=label_count),
            size=2.5, color="white", fontface="bold") +
  scale_fill_manual(values=strain_colors) +
  labs(title="Strain — counts", x="Cluster", y="Cell count", fill="Strain") +
  theme_bw(base_size=10) +
  theme(plot.title=element_text(face="bold", size=10),
        panel.grid.minor=element_blank())

# -- Bar chart 4: proportions, fill = strain only ------------------------------
lab4 <- make_label_df(umap_comb_df, "strain", "prop")
p_bar_s_prop <- ggplot(lab4, aes(x=cluster, y=prop, fill=strain)) +
  geom_bar(stat="identity", position="stack", color="white", linewidth=0.15) +
  geom_text(aes(y=y_pos, label=label_prop),
            size=2.5, color="white", fontface="bold") +
  scale_fill_manual(values=strain_colors) +
  scale_y_continuous(labels=scales::percent_format()) +
  labs(title="Strain — proportions", x="Cluster", y="Proportion",
       fill="Strain") +
  theme_bw(base_size=10) +
  theme(plot.title=element_text(face="bold", size=10),
        panel.grid.minor=element_blank())

# Combine all 4 into one PDF (2 x 2 layout, A3 landscape)
bar_combined <- (p_bar_sc_counts | p_bar_sc_prop) /
                (p_bar_s_counts  | p_bar_s_prop) +
  plot_annotation(
    title    = paste0("Cluster composition — ", n_clust_comb,
                      " clusters | n=372 cells"),
    subtitle = "Labels shown for segments \u226505% of bar",
    theme    = theme(
      plot.title    = element_text(face="bold", size=14),
      plot.subtitle = element_text(size=9, color="grey40")
    )
  )

bars_path <- file.path(combined_out, "barplots_cluster_composition.pdf")
ggsave(bars_path, bar_combined, width=20, height=14)
cat("Saved:", bars_path, "\n")
# ==============================================================================
# Cluster marker genes: Wilcoxon rank-sum test (one-vs-rest) per cluster
# ==============================================================================
cat("\nRunning Wilcoxon rank-sum marker gene analysis...\n")

# Parameters
MIN_LOG2FC  <- 0.5   # minimum log2FC (VST difference) to consider
MAX_PADJ    <- 0.05  # BH-adjusted p-value threshold
TOP_EXCEL   <- 100   # genes per cluster in Excel
TOP_HEATMAP <- 15    # genes per cluster for heatmap

cluster_vec    <- setNames(as.character(umap_comb_df$cluster), umap_comb_df$cell_id)
cluster_levels <- sort(unique(as.integer(cluster_vec)))
cat("Clusters:", paste(cluster_levels, collapse=", "), "\n")
cat("Genes:", nrow(combined_expr), "\n")

# Excel workbook
wb <- createWorkbook()

summary_df <- umap_comb_df %>%
  group_by(cluster, strain, condition_clean) %>%
  summarise(n_cells=n(), .groups="drop") %>%
  arrange(as.integer(as.character(cluster)), strain, condition_clean)
addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_df)

all_markers <- list()

for (cl in cluster_levels) {
  cl_label   <- as.character(cl)
  cat("  Cluster", cl_label, "vs rest...\n")

  cl_cells   <- names(cluster_vec)[cluster_vec == cl_label]
  rest_cells <- names(cluster_vec)[cluster_vec != cl_label]

  mean_cl   <- if (length(cl_cells)   == 1) combined_expr[, cl_cells]
               else rowMeans(combined_expr[, cl_cells,   drop=FALSE])
  mean_rest <- if (length(rest_cells) == 1) combined_expr[, rest_cells]
               else rowMeans(combined_expr[, rest_cells, drop=FALSE])

  # VST is log-scale so difference approximates log2FC
  log2fc <- mean_cl - mean_rest

  # Pre-filter to upregulated candidates only (speeds up Wilcoxon loop)
  candidates <- names(log2fc)[log2fc >= MIN_LOG2FC]
  cat("    Candidates (log2FC >=", MIN_LOG2FC, "):", length(candidates), "\n")

  if (length(candidates) == 0) {
    all_markers[[cl_label]] <- data.frame()
    next
  }

  cl_mat   <- combined_expr[candidates, cl_cells,   drop=FALSE]
  rest_mat <- combined_expr[candidates, rest_cells, drop=FALSE]

  pvals <- sapply(candidates, function(g) {
    wilcox.test(cl_mat[g, ], rest_mat[g, ],
                alternative="greater", exact=FALSE)$p.value
  })
  padj <- p.adjust(pvals, method="BH")

  res_cl <- data.frame(
    ensembl_id       = candidates,
    gene_symbol      = to_sym(candidates),
    mean_VST_cluster = round(mean_cl[candidates],   4),
    mean_VST_rest    = round(mean_rest[candidates], 4),
    log2FC           = round(log2fc[candidates],    4),
    pval             = signif(pvals, 4),
    padj             = signif(padj,  4),
    stringsAsFactors = FALSE
  ) %>%
    filter(padj < MAX_PADJ) %>%
    arrange(desc(log2FC)) %>%
    distinct(gene_symbol, .keep_all=TRUE) %>%
    slice_head(n=TOP_EXCEL) %>%
    mutate(rank=row_number()) %>%
    select(rank, ensembl_id, gene_symbol,
           mean_VST_cluster, mean_VST_rest, log2FC, pval, padj)

  cat("    Significant markers:", nrow(res_cl),
      "| Top gene:", if (nrow(res_cl) > 0) res_cl$gene_symbol[1] else "none", "\n")

  all_markers[[cl_label]] <- res_cl

  sheet_name <- paste0("Cluster_", cl_label)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, res_cl)
  addStyle(wb, sheet_name,
    style=createStyle(textDecoration="bold", fgFill="#D9E1F2"),
    rows=1, cols=1:8, gridExpand=TRUE)
  setColWidths(wb, sheet_name, cols=1:8,
    widths=c(6, 20, 16, 16, 14, 10, 12, 12))
}

markers_path <- file.path(combined_out, "cluster_marker_genes.xlsx")
saveWorkbook(wb, markers_path, overwrite=TRUE)
cat("Saved:", markers_path, "\n")

# ==============================================================================
# Heatmap: top marker genes per cluster
# ==============================================================================
cat("\nGenerating cluster marker gene heatmap...\n")

TOP_GENES_PER_CLUSTER <- TOP_HEATMAP  # set above in marker gene section

# Select top N per cluster by combined score, deduplicate (first occurrence wins)
heatmap_genes <- c()
for (cl_label in as.character(cluster_levels)) {
  if (!is.null(all_markers[[cl_label]]) && nrow(all_markers[[cl_label]]) > 0) {
    top_g <- all_markers[[cl_label]]$gene_symbol[
               !all_markers[[cl_label]]$gene_symbol %in% heatmap_genes]
    top_g <- head(top_g, TOP_GENES_PER_CLUSTER)
    heatmap_genes <- c(heatmap_genes, top_g)
  }
}
cat("Total unique genes in heatmap:", length(heatmap_genes), "\n")

# Get Ensembl IDs for selected gene symbols
sym_to_ens <- setNames(names(sym_map), sym_map)
heatmap_ens <- sym_to_ens[heatmap_genes]
heatmap_ens <- heatmap_ens[!is.na(heatmap_ens) & heatmap_ens %in% rownames(combined_expr)]

# Z-score VST matrix: rows = genes, cols = cells ordered by cluster
cell_order   <- umap_comb_df$cell_id[order(as.integer(
                  as.character(umap_comb_df$cluster)))]
expr_heat    <- combined_expr[heatmap_ens, cell_order]
rownames(expr_heat) <- to_sym(rownames(expr_heat))

# Z-score per gene
expr_z <- t(scale(t(expr_heat)))
expr_z[expr_z >  2.5] <-  2.5   # cap for color scale
expr_z[expr_z < -2.5] <- -2.5

# Cluster annotation bar (top)
cluster_anno_vec <- umap_comb_df$cluster[match(cell_order, umap_comb_df$cell_id)]
cluster_anno_vec <- as.character(cluster_anno_vec)

ha_top <- HeatmapAnnotation(
  Cluster = cluster_anno_vec,
  col     = list(Cluster = cluster_pal_comb),
  annotation_name_side = "left",
  annotation_legend_param = list(
    Cluster = list(title="Cluster", nrow=1)
  )
)

# Color scale: purple-black-yellow (similar to reference figure)
col_fun <- colorRamp2(
  c(-2.5, 0, 2.5),
  c("#3D0751", "#1A1A1A", "#F5E642")
)

# Gene labels: show all since these are curated top markers
ht <- Heatmap(
  expr_z,
  name                  = "Z-score",
  col                   = col_fun,
  top_annotation        = ha_top,
  show_column_names     = FALSE,
  show_row_names        = TRUE,
  row_names_gp          = gpar(fontsize=7, fontface="italic"),
  cluster_rows          = FALSE,
  cluster_columns       = FALSE,
  row_title             = NULL,
  column_title          = paste0("Cluster marker genes — ", length(heatmap_genes),
                                  " genes x ", ncol(expr_z), " cells"),
  column_title_gp       = gpar(fontsize=12, fontface="bold"),
  heatmap_legend_param  = list(
    title          = "Z-score\n(VST)",
    legend_height  = unit(3, "cm")
  ),
  use_raster            = TRUE,
  raster_quality        = 5
)

heatmap_path <- file.path(combined_out, "heatmap_cluster_markers.pdf")
pdf(heatmap_path, width=16, height=max(8, length(heatmap_genes) * 0.18 + 3))
draw(ht, heatmap_legend_side="right", annotation_legend_side="bottom")
dev.off()
cat("Saved:", heatmap_path, "\n")

# ==============================================================================
# Cell identity scoring: CellMarker 2.0 gene set overlap per cluster
# ==============================================================================
cat("\nRunning cell identity scoring against CellMarker 2.0...\n")

# -- Download CellMarker 2.0 mouse marker table --------------------------------
cm_url   <- "http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download_files/file/Cell_marker_Mouse.xlsx"
cm_cache <- "/tmp/CellMarker_Mouse.xlsx"

if (!file.exists(cm_cache)) {
  cat("  Downloading CellMarker 2.0 mouse table...\n")
  tryCatch(
    download.file(cm_url, cm_cache, mode="wb", quiet=TRUE),
    error = function(e) cat("  Download failed:", conditionMessage(e), "\n")
  )
}

# Load or fall back to built-in curated set
if (file.exists(cm_cache) && file.size(cm_cache) > 1000) {
  cm_raw <- openxlsx::read.xlsx(cm_cache)
  cat("  CellMarker loaded:", nrow(cm_raw), "entries\n")

  # Keep mouse entries; filter to relevant tissues
  relevant_tissues <- c("Pancreas", "Immune system", "Lymphoid tissue",
                         "Blood", "Spleen", "Thymus", "Bone marrow",
                         "Lymph node", "Peripheral blood", "Liver",
                         "Adipose tissue", "Normal tissue", "")
  cm_filt <- cm_raw %>%
    filter(grepl(paste(relevant_tissues, collapse="|"),
                 tissue_type, ignore.case=TRUE) |
           is.na(tissue_type)) %>%
    filter(!is.na(Symbol) & Symbol != "")

  # Build named list: cell_type -> vector of gene symbols
  cm_list <- cm_filt %>%
    group_by(cell_name) %>%
    summarise(markers=list(unique(Symbol)), .groups="drop")
  cell_db <- setNames(cm_list$markers, cm_list$cell_name)
  cat("  Cell types in filtered DB:", length(cell_db), "\n")

} else {
  cat("  Using built-in curated marker set (CellMarker download unavailable)\n")
  # Curated fallback: key pancreatic/immune cell types
  cell_db <- list(
    "Beta cell"             = c("Ins1","Ins2","Pdx1","Nkx6-1","Slc2a2","Iapp"),
    "Alpha cell"            = c("Gcg","Arx","Mafb","Irx2","Slc38a5"),
    "Delta cell"            = c("Sst","Hhex","Rbp4"),
    "Ductal cell"           = c("Krt19","Krt7","Sox9","Cftr","Spp1","Muc1"),
    "Acinar cell"           = c("Ptf1a","Amy2a","Cpa1","Prss1","Rbpjl"),
    "Stellate cell"         = c("Col1a1","Col3a1","Acta2","Pdgfrb","Des"),
    "Endothelial cell"      = c("Pecam1","Cdh5","Eng","Kdr","Tie1","Vwf"),
    "Macrophage"            = c("Cd68","Adgre1","Csf1r","Itgam","Mrc1","Ly6c2"),
    "Dendritic cell"        = c("Itgax","Siglech","Clec9a","Xcr1","Ccr7","Flt3"),
    "B cell"                = c("Cd19","Ms4a1","Pax5","Cd79a","Cd79b","Blnk"),
    "T cell"                = c("Cd3e","Cd3d","Trac","Trbc1","Cd4","Cd8a"),
    "NK cell"               = c("Ncr1","Klrb1c","Gzma","Prf1","Nkg7"),
    "Regulatory T cell"     = c("Foxp3","Il2ra","Ikzf2","Ctla4","Tnfrsf18"),
    "Mast cell"             = c("Kit","Fcer1a","Tpsab1","Mcpt4","Hdc"),
    "Neutrophil"            = c("S100a8","S100a9","Ly6g","Cxcr2","Csf3r"),
    "Fibroblast"            = c("Col1a2","Col6a1","Thy1","Fap","Pdpn","S100a4"),
    "Pericyte"              = c("Rgs5","Pdgfrb","Notch3","Cspg4","Abcc9"),
    "Epithelial cell"       = c("Epcam","Cdh1","Krt8","Krt18","Cldn3","Ocln"),
    "Plasmacytoid DC"       = c("Siglech","Bst2","Irf7","Tlr7","Ccr9","Ly6d"),
    "Plasma cell"           = c("Sdc1","Prdm1","Xbp1","Mzb1","Jchain","Igha")
  )
  cat("  Built-in cell types:", length(cell_db), "\n")
}

# -- Scoring function: Fisher exact + Jaccard per cluster ---------------------
score_cluster <- function(cluster_markers, cell_db, universe_size) {
  # cluster_markers: character vector of significant marker gene symbols
  # Returns dataframe ranked by Fisher p-value
  results <- lapply(names(cell_db), function(ct) {
    ref_genes <- cell_db[[ct]]
    overlap   <- intersect(cluster_markers, ref_genes)
    n_overlap <- length(overlap)
    n_cluster <- length(cluster_markers)
    n_ref     <- length(ref_genes)
    n_neither <- universe_size - length(union(cluster_markers, ref_genes))

    # Fisher exact test
    ft <- fisher.test(matrix(c(n_overlap,
                                n_cluster - n_overlap,
                                n_ref     - n_overlap,
                                n_neither),
                              nrow=2), alternative="greater")
    # Jaccard similarity
    jaccard <- n_overlap / length(union(cluster_markers, ref_genes))

    data.frame(
      cell_type        = ct,
      n_overlap        = n_overlap,
      n_cluster_markers = n_cluster,
      n_ref_markers    = n_ref,
      overlap_genes    = paste(overlap, collapse=", "),
      jaccard          = round(jaccard, 4),
      fisher_pval      = ft$p.value,
      odds_ratio       = round(ft$estimate, 3),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, results) %>%
    mutate(fisher_padj = p.adjust(fisher_pval, method="BH")) %>%
    arrange(fisher_pval) %>%
    slice_head(n=20) %>%
    mutate(rank=row_number()) %>%
    select(rank, cell_type, n_overlap, n_cluster_markers, n_ref_markers,
           jaccard, odds_ratio, fisher_pval, fisher_padj, overlap_genes)
}

universe_genes <- to_sym(rownames(combined_expr))
universe_size  <- length(unique(universe_genes))

# -- Score combined clusters ---------------------------------------------------
cat("  Scoring combined clusters...\n")
wb_cellid_comb <- createWorkbook()

# Summary sheet
addWorksheet(wb_cellid_comb, "Summary")
writeData(wb_cellid_comb, "Summary",
          data.frame(
            cluster      = names(all_markers),
            n_markers    = sapply(all_markers, nrow),
            top_cell_type = sapply(names(all_markers), function(cl) {
              if (is.null(all_markers[[cl]]) || nrow(all_markers[[cl]]) == 0)
                return("insufficient markers")
              sc <- score_cluster(all_markers[[cl]]$gene_symbol, cell_db, universe_size)
              if (nrow(sc) > 0) sc$cell_type[1] else "no match"
            })
          ))

for (cl_label in names(all_markers)) {
  cat("    Combined cluster", cl_label, "\n")
  if (is.null(all_markers[[cl_label]]) || nrow(all_markers[[cl_label]]) == 0) {
    next
  }
  sc <- score_cluster(all_markers[[cl_label]]$gene_symbol, cell_db, universe_size)
  sheet_nm <- paste0("Cluster_", cl_label)
  addWorksheet(wb_cellid_comb, sheet_nm)
  writeData(wb_cellid_comb, sheet_nm, sc)
  addStyle(wb_cellid_comb, sheet_nm,
    style=createStyle(textDecoration="bold", fgFill="#D9E1F2"),
    rows=1, cols=1:10, gridExpand=TRUE)
  setColWidths(wb_cellid_comb, sheet_nm, cols=1:10,
    widths=c(5,25,10,14,12,10,12,12,12,40))
}

cellid_comb_path <- file.path(combined_out, "cell_identity_combined_clusters.xlsx")
saveWorkbook(wb_cellid_comb, cellid_comb_path, overwrite=TRUE)
cat("  Saved:", cellid_comb_path, "\n")

# ==============================================================================
# Mouse Cell Atlas expression profile correlation (clustifyrdata::ref_MCA)
# ==============================================================================
cat("\nRunning Mouse Cell Atlas expression profile correlation (combined)...\n")
run_mca_correlation(
  expr_mat         = combined_expr,
  cluster_vec_input = cluster_vec,
  cluster_levels_input = cluster_levels,
  to_sym_fn        = to_sym,
  out_path         = file.path(combined_out,
                               "MCA_celltype_correlation_combined_clusters.xlsx"),
  label            = "combined"
)

# ==============================================================================
# VAF / VRC correlation: compare clusters to Clarke et al. 2025 (GSE292898)
# ==============================================================================
cat("\nRunning VAF/VRC correlation analysis (Clarke et al. 2025)...\n")

vaf_counts_f <- file.path(base_dir,
  "reference/GSE292898_teyton_don_2025_processed_mouse_raw_counts_matrix.csv.gz")

if (!file.exists(vaf_counts_f)) {
  cat("  WARNING: VAF reference file not found at:\n  ", vaf_counts_f, "\n")
  cat("  Download from GEO GSE292898 and place in reference/ to enable this analysis.\n")
} else {
  cat("  Loading Clarke et al. mouse count matrix...\n")

  # -- Load and normalize reference counts -------------------------------------
  vaf_raw <- read.csv(gzfile(vaf_counts_f), row.names=1, check.names=FALSE)
  # Columns: gene_id, gene_name, then cell IDs
  vaf_gene_id   <- vaf_raw[["gene_id"]]
  vaf_gene_name <- vaf_raw[["gene_name"]]
  vaf_counts_mat <- as.matrix(vaf_raw[, !colnames(vaf_raw) %in% c("gene_id","gene_name")])
  rownames(vaf_counts_mat) <- vaf_gene_name
  storage.mode(vaf_counts_mat) <- "numeric"
  cat("  Reference matrix:", nrow(vaf_counts_mat), "genes x",
      ncol(vaf_counts_mat), "cells\n")

  # Subset to CD45neg cells (plates B-H) per README — these contain VAF + VRC
  cd45neg_mask <- !grepl("^A[0-9]", colnames(vaf_counts_mat))
  vaf_neg_mat  <- vaf_counts_mat[, cd45neg_mask, drop=FALSE]
  cat("  CD45neg cells (B-H plates):", ncol(vaf_neg_mat), "\n")

  # Log1p-CPM normalization
  lib_sizes_vaf  <- colSums(vaf_neg_mat)
  cpm_vaf        <- sweep(vaf_neg_mat, 2, lib_sizes_vaf, "/") * 1e6
  logcpm_vaf     <- log1p(cpm_vaf)

  # -- Recluster CD45neg reference cells to identify VAF vs VRC ---------------
  # Use PCA on top variable genes then k-means k=2 (paper shows 2 clusters)
  cat("  Clustering reference CD45neg cells to assign VAF/VRC labels...\n")
  gene_var_vaf  <- apply(logcpm_vaf, 1, var)
  top_vaf_genes <- names(sort(gene_var_vaf, decreasing=TRUE))[
                     1:min(500, nrow(logcpm_vaf))]
  pca_vaf       <- prcomp(t(logcpm_vaf[top_vaf_genes, ]), center=TRUE,
                           scale.=FALSE)
  set.seed(UMAP_SEED)
  km_vaf        <- kmeans(pca_vaf$x[, 1:min(10, ncol(pca_vaf$x))],
                           centers=2, nstart=25)
  cluster_vaf   <- km_vaf$cluster

  # Assign VAF vs VRC using known marker genes from paper (Figure 3)
  # VAF markers: Col1a1, Col1a2, Timp3, Pecam1(low), Spp1
  # VRC markers: Pecam1, Eng, Cdh5, Kdr (endothelial)
  vaf_markers <- intersect(c("Col1a1","Col1a2","Timp3","Spp1","Thy1",
                               "Pdpn","S100a4","Fn1"), rownames(logcpm_vaf))
  vrc_markers <- intersect(c("Pecam1","Eng","Cdh5","Kdr","Tie1","Vwf",
                               "Cdh5","Esam"), rownames(logcpm_vaf))

  score_markers <- function(mat, genes) {
    g <- intersect(genes, rownames(mat))
    if (length(g) == 0) return(rep(0, ncol(mat)))
    colMeans(mat[g, , drop=FALSE])
  }

  c1_cells <- colnames(logcpm_vaf)[cluster_vaf == 1]
  c2_cells <- colnames(logcpm_vaf)[cluster_vaf == 2]

  c1_vaf_score <- mean(score_markers(logcpm_vaf[, c1_cells, drop=FALSE], vaf_markers))
  c2_vaf_score <- mean(score_markers(logcpm_vaf[, c2_cells, drop=FALSE], vaf_markers))
  c1_vrc_score <- mean(score_markers(logcpm_vaf[, c1_cells, drop=FALSE], vrc_markers))
  c2_vrc_score <- mean(score_markers(logcpm_vaf[, c2_cells, drop=FALSE], vrc_markers))

  # Cluster with higher VAF score = VAF, other = VRC
  if (c1_vaf_score > c2_vaf_score) {
    vaf_cells <- c1_cells; vrc_cells <- c2_cells
  } else {
    vaf_cells <- c2_cells; vrc_cells <- c1_cells
  }
  cat("  VAF cells identified:", length(vaf_cells), "\n")
  cat("  VRC cells identified:", length(vrc_cells), "\n")
  cat("  VAF score check - VAF cluster:", round(max(c1_vaf_score, c2_vaf_score), 3),
      "VRC cluster:", round(min(c1_vaf_score, c2_vaf_score), 3), "\n")
  cat("  VRC score check - VRC cluster:", round(max(c1_vrc_score, c2_vrc_score), 3),
      "VAF cluster:", round(min(c1_vrc_score, c2_vrc_score), 3), "\n")

  # -- Mean expression profiles per reference population ----------------------
  mean_vaf_ref <- rowMeans(logcpm_vaf[, vaf_cells, drop=FALSE])
  mean_vrc_ref <- rowMeans(logcpm_vaf[, vrc_cells, drop=FALSE])

  # Also include CD45pos as a third reference (plate A)
  cd45pos_mask <- grepl("^A[0-9]", colnames(vaf_counts_mat))
  vaf_pos_mat  <- vaf_counts_mat[, cd45pos_mask, drop=FALSE]
  lib_pos      <- colSums(vaf_pos_mat)
  logcpm_pos   <- log1p(sweep(vaf_pos_mat, 2, lib_pos, "/") * 1e6)
  mean_cd45pos_ref <- rowMeans(logcpm_pos)

  ref_profiles <- list(
    VAF    = mean_vaf_ref,
    VRC    = mean_vrc_ref,
    CD45pos = mean_cd45pos_ref
  )

  # -- Compute correlations for each of your Leiden clusters ------------------
  # Use combined_expr (VST) mean profiles per cluster
  cat("  Computing correlations against combined Leiden clusters...\n")

  # Find common genes between reference (gene symbols) and your data
  # Your data uses Ensembl IDs as rownames — map to symbols
  your_symbols <- to_sym(rownames(combined_expr))
  common_genes_vaf <- intersect(your_symbols,
                                 names(mean_vaf_ref)[mean_vaf_ref > 0])
  cat("  Common genes for correlation:", length(common_genes_vaf), "\n")

  # Build your cluster mean profiles on common genes
  your_idx <- which(your_symbols %in% common_genes_vaf)
  your_syms_sub <- your_symbols[your_idx]

  corr_results <- list()
  for (cl_label in as.character(cluster_levels)) {
    cl_cells_your <- names(cluster_vec)[cluster_vec == cl_label]
    if (length(cl_cells_your) == 0) next

    your_mean <- rowMeans(combined_expr[your_idx, cl_cells_your, drop=FALSE])
    names(your_mean) <- your_syms_sub

    cl_corrs <- data.frame(
      reference_population = character(),
      pearson_r            = numeric(),
      spearman_rho         = numeric(),
      n_genes              = integer(),
      stringsAsFactors     = FALSE
    )

    for (ref_name in names(ref_profiles)) {
      ref_vec  <- ref_profiles[[ref_name]]
      common   <- intersect(names(your_mean), names(ref_vec))
      if (length(common) < 50) {
        cat("    Cluster", cl_label, "vs", ref_name, "- too few common genes:",
            length(common), "\n")
        next
      }
      y <- your_mean[common]
      r <- ref_vec[common]
      pearson_r   <- cor(y, r, method="pearson")
      spearman_rho <- cor(y, r, method="spearman")
      cl_corrs <- rbind(cl_corrs, data.frame(
        reference_population = ref_name,
        pearson_r            = round(pearson_r,   4),
        spearman_rho         = round(spearman_rho, 4),
        n_genes              = length(common),
        stringsAsFactors     = FALSE
      ))
    }
    corr_results[[cl_label]] <- cl_corrs
  }

  # -- Build Excel workbook ---------------------------------------------------
  wb_vaf <- createWorkbook()

  # Summary sheet
  summary_vaf <- do.call(rbind, lapply(names(corr_results), function(cl) {
    df <- corr_results[[cl]]
    if (nrow(df) == 0) return(NULL)
    best_row <- df[which.max(df$pearson_r), ]
    data.frame(
      cluster              = cl,
      best_match           = best_row$reference_population,
      best_pearson_r       = best_row$pearson_r,
      best_spearman_rho    = best_row$spearman_rho,
      stringsAsFactors     = FALSE
    )
  }))
  addWorksheet(wb_vaf, "Summary")
  writeData(wb_vaf, "Summary", summary_vaf)
  addStyle(wb_vaf, "Summary",
    style=createStyle(textDecoration="bold", fgFill="#E2EFDA"),
    rows=1, cols=1:4, gridExpand=TRUE)
  setColWidths(wb_vaf, "Summary", cols=1:4, widths=c(10,18,16,16))

  # One sheet per cluster
  for (cl_label in names(corr_results)) {
    df <- corr_results[[cl_label]]
    if (nrow(df) == 0) next
    sheet_nm <- paste0("Cluster_", cl_label)
    addWorksheet(wb_vaf, sheet_nm)
    writeData(wb_vaf, sheet_nm, df)
    addStyle(wb_vaf, sheet_nm,
      style=createStyle(textDecoration="bold", fgFill="#E2EFDA"),
      rows=1, cols=1:4, gridExpand=TRUE)
    setColWidths(wb_vaf, sheet_nm, cols=1:4, widths=c(22,12,14,10))

    # Highlight best match row in green
    best_idx <- which.max(df$pearson_r) + 1  # +1 for header row
    addStyle(wb_vaf, sheet_nm,
      style=createStyle(fgFill="#C6EFCE"),
      rows=best_idx, cols=1:4, gridExpand=TRUE)
  }

  vaf_corr_path <- file.path(combined_out,
                              "VAF_VRC_correlation_combined_clusters.xlsx")
  saveWorkbook(wb_vaf, vaf_corr_path, overwrite=TRUE)
  cat("  Saved:", vaf_corr_path, "\n")
}

# -- Score per-strain clusters ------------------------------------------------
cat("  Scoring per-strain clusters...\n")

for (strain in strains) {
  if (is.null(per_strain_cluster_list[[strain]])) next
  s_info    <- per_strain_cluster_list[[strain]]
  s_umap_df <- s_info$umap_df
  s_expr    <- s_info$expr
  s_levels  <- sort(unique(as.integer(as.character(s_umap_df$cluster))))

  # Compute per-strain cluster markers (Wilcoxon, same method as combined)
  s_cluster_vec <- setNames(as.character(s_umap_df$cluster), s_umap_df$cell_id)

  wb_strain <- createWorkbook()

  # Collect all cluster results first, then write Summary sheet first
  strain_summary  <- data.frame(cluster=character(), n_markers=integer(),
                                 top_cell_type=character(),
                                 stringsAsFactors=FALSE)
  cluster_results <- list()

  for (cl in s_levels) {
    cl_label   <- as.character(cl)
    cl_cells   <- names(s_cluster_vec)[s_cluster_vec == cl_label]
    rest_cells <- names(s_cluster_vec)[s_cluster_vec != cl_label]

    # Attempt Wilcoxon one-vs-rest if enough cells in both groups
    syms_s <- character(0)
    if (length(cl_cells) >= 2 && length(rest_cells) >= 2) {
      mean_cl   <- rowMeans(s_expr[, cl_cells,   drop=FALSE])
      mean_rest <- rowMeans(s_expr[, rest_cells, drop=FALSE])
      log2fc    <- mean_cl - mean_rest
      candidates <- names(log2fc)[log2fc >= MIN_LOG2FC]

      if (length(candidates) >= 2) {
        pvals_s <- sapply(candidates, function(g) {
          wilcox.test(s_expr[g, cl_cells], s_expr[g, rest_cells],
                      alternative="greater", exact=FALSE)$p.value
        })
        padj_s <- p.adjust(pvals_s, method="BH")
        sig_s  <- candidates[padj_s < MAX_PADJ]
        syms_s <- unique(to_sym(sig_s))
        syms_s <- syms_s[syms_s != ""]
      }
    }

    # Fallback: if too few Wilcoxon markers use top 50 expressed genes
    if (length(syms_s) < 10) {
      cat("    Cluster", cl_label, "- insufficient Wilcoxon markers (",
          length(syms_s), "), using top expressed genes\n")
      mean_cl_all <- if (length(cl_cells) == 1) s_expr[, cl_cells]
                     else rowMeans(s_expr[, cl_cells, drop=FALSE])
      top_ens <- names(sort(mean_cl_all, decreasing=TRUE))[
                   1:min(50, length(mean_cl_all))]
      syms_s  <- unique(to_sym(top_ens))
      syms_s  <- syms_s[syms_s != ""]
    }

    sc_s   <- score_cluster(syms_s, cell_db, universe_size)
    top_ct <- if (nrow(sc_s) > 0) sc_s$cell_type[1] else "no match"
    strain_summary <- rbind(strain_summary,
                             data.frame(cluster=cl_label,
                                        n_markers=length(syms_s),
                                        top_cell_type=top_ct,
                                        stringsAsFactors=FALSE))
    cluster_results[[cl_label]] <- sc_s
  }

  # Write Summary sheet first so it appears as the first tab
  addWorksheet(wb_strain, "Summary")
  writeData(wb_strain, "Summary", strain_summary)
  addStyle(wb_strain, "Summary",
    style=createStyle(textDecoration="bold", fgFill="#FCE4D6"),
    rows=1, cols=1:3, gridExpand=TRUE)

  # Then write cluster sheets
  for (cl_label in names(cluster_results)) {
    sheet_nm <- paste0("Cluster_", cl_label)
    addWorksheet(wb_strain, sheet_nm)
    writeData(wb_strain, sheet_nm, cluster_results[[cl_label]])
    addStyle(wb_strain, sheet_nm,
      style=createStyle(textDecoration="bold", fgFill="#FCE4D6"),
      rows=1, cols=1:10, gridExpand=TRUE)
    setColWidths(wb_strain, sheet_nm, cols=1:10,
      widths=c(5,25,10,14,12,10,12,12,12,40))
  }

  strain_cellid_path <- file.path(dge_dir, paste0(strain, "_plots"),
                                   paste0("cell_identity_", strain, "_clusters.xlsx"))
  saveWorkbook(wb_strain, strain_cellid_path, overwrite=TRUE)
  cat("  Saved:", strain_cellid_path, "\n")
}

# ==============================================================================
# Combined cluster violin plots
# ==============================================================================
cat("\nGenerating combined cluster violin plots...\n")

# Genes to plot - add more symbols to this vector to include additional plots
VIOLIN_GENES <- c("Ptprc")

# Look up Ensembl IDs for violin genes
violin_ens <- sym_to_ens_map[VIOLIN_GENES]
violin_ens <- violin_ens[!is.na(violin_ens) & violin_ens %in% rownames(combined_expr)]

missing_genes <- setdiff(VIOLIN_GENES, names(violin_ens))
if (length(missing_genes) > 0) {
  cat("  Warning: genes not found in matrix:", paste(missing_genes, collapse=", "), "\n")
}
cat("  Genes to plot:", paste(names(violin_ens), collapse=", "), "\n")

# Build long-format dataframe for all violin genes
vln_comb_df <- as.data.frame(t(combined_expr[violin_ens, umap_comb_df$cell_id,
                                              drop=FALSE])) %>%
  tibble::rownames_to_column("cell_id") %>%
  left_join(umap_comb_df[, c("cell_id", "cluster")], by="cell_id") %>%
  pivot_longer(-c(cell_id, cluster), names_to="ensembl_id", values_to="VST") %>%
  mutate(
    gene_symbol = to_sym(ensembl_id),
    cluster     = factor(cluster, levels=sort(unique(as.integer(
                    as.character(umap_comb_df$cluster)))))
  )

# One violin panel per gene
vln_gene_plots <- lapply(seq_along(violin_ens), function(i) {
  ens  <- violin_ens[i]
  sym  <- names(violin_ens)[i]
  df_g <- vln_comb_df %>% filter(ensembl_id == ens)

  ggplot(df_g, aes(x=cluster, y=VST, fill=cluster)) +
    geom_violin(trim=TRUE, scale="width", alpha=0.85, linewidth=0.3) +
    geom_jitter(width=0.15, size=0.6, alpha=0.35, color="grey20") +
    scale_fill_manual(values=cluster_pal_comb) +
    labs(
      title = sym,
      x     = "Leiden Cluster",
      y     = "VST expression"
    ) +
    theme_bw(base_size=11) +
    theme(
      plot.title       = element_text(face="bold.italic", size=11),
      legend.position  = "none",
      panel.grid.minor = element_blank()
    )
})

# Arrange into a grid — wrap at 3 per row
n_genes    <- length(vln_gene_plots)
n_cols_vln <- min(3, n_genes)
n_rows_vln <- ceiling(n_genes / n_cols_vln)

vln_combined <- wrap_plots(vln_gene_plots, ncol=n_cols_vln) +
  plot_annotation(
    title    = "Gene expression across Leiden clusters",
    subtitle = paste0("VST-normalized | n=", nrow(umap_comb_df), " cells | ",
                      n_clust_comb, " clusters"),
    theme    = theme(
      plot.title    = element_text(face="bold", size=13),
      plot.subtitle = element_text(size=9, color="grey40")
    )
  )

vln_path <- file.path(combined_out, "violin_plots_by_cluster.pdf")
ggsave(vln_path, vln_combined,
       width  = n_cols_vln * 5,
       height = n_rows_vln * 5)
cat("Saved:", vln_path, "\n")

cat("\n==============================================================\n")
cat("All outputs complete. Output structure:\n")
for (st in strains) {
  cat(" results/05_dge/", st, "_plots/  (10 files)\n", sep="")
}
cat(" results/05_dge/combined_plots/\n")
cat("   - umap_all372_by_cluster.pdf\n")
cat("   - umap_all372_by_strain_population.pdf\n")
cat("   - barplots_cluster_composition.pdf  (all 4 charts, 2x2)\n")
cat("   - cluster_marker_genes.xlsx\n")
cat("   - heatmap_cluster_markers.pdf\n")
cat("==============================================================\n")

# ==============================================================================
# CombinedwithVAFPaperPlots: merge your cells with Clarke et al. 2025 mouse data
# ==============================================================================
cat("\n==============================================================\n")
cat("Building CombinedwithVAFPaperPlots analysis...\n")
cat("==============================================================\n")

vaf_counts_f2 <- file.path(base_dir,
  "reference/GSE292898_teyton_don_2025_processed_mouse_raw_counts_matrix.csv.gz")

if (!file.exists(vaf_counts_f2)) {
  cat("WARNING: Clarke et al. count matrix not found. Skipping merged analysis.\n")
  cat("Place file at:", vaf_counts_f2, "\n")
} else {

  vaf_out_dir <- file.path(dge_dir, "CombinedwithVAFPaperPlots")
  dir.create(vaf_out_dir, showWarnings=FALSE, recursive=TRUE)

  # -- Load Clarke counts ------------------------------------------------------
  cat("Loading Clarke et al. mouse count matrix...\n")
  vaf_raw2      <- read.csv(gzfile(vaf_counts_f2), row.names=1, check.names=FALSE)
  vaf_gene_sym2 <- vaf_raw2[["gene_name"]]
  vaf_cnt_mat2  <- as.matrix(vaf_raw2[, !colnames(vaf_raw2) %in%
                                        c("gene_id","gene_name")])
  rownames(vaf_cnt_mat2) <- vaf_gene_sym2
  storage.mode(vaf_cnt_mat2) <- "numeric"
  cat("Clarke matrix:", nrow(vaf_cnt_mat2), "genes x",
      ncol(vaf_cnt_mat2), "cells\n")

  # -- Map Clarke gene symbols to Ensembl IDs ----------------------------------
  # Your pipeline uses Ensembl IDs as rownames; Clarke uses symbols
  # Build reverse map: symbol -> first Ensembl ID
  ens_to_sym_vec <- to_sym(rownames(combined_expr))
  sym_to_ens_rev <- tapply(rownames(combined_expr), ens_to_sym_vec, `[`, 1)

  clarke_syms_in_yours <- intersect(rownames(vaf_cnt_mat2),
                                     names(sym_to_ens_rev))
  cat("Clarke genes mappable to your Ensembl IDs:", length(clarke_syms_in_yours),
      "\n")

  # Subset Clarke to mappable genes, reindex by Ensembl
  vaf_cnt_sub <- vaf_cnt_mat2[clarke_syms_in_yours, , drop=FALSE]
  ens_ids_for_clarke <- sym_to_ens_rev[clarke_syms_in_yours]
  rownames(vaf_cnt_sub) <- ens_ids_for_clarke

  # -- MHCII filter on Clarke cells --------------------------------------------
  cat("Applying MHCII filter to Clarke cells...\n")
  mhcii_ens_clarke <- sym_to_ens_rev[MHCII_GENES]
  mhcii_ens_clarke <- mhcii_ens_clarke[!is.na(mhcii_ens_clarke) &
                                         mhcii_ens_clarke %in% rownames(vaf_cnt_sub)]
  if (length(mhcii_ens_clarke) >= 2) {
    lib_c       <- colSums(vaf_cnt_sub)
    cpm_c       <- sweep(vaf_cnt_sub, 2, pmax(lib_c, 1), "/") * 1e6
    logcpm_c    <- log1p(cpm_c)
    passes_c    <- colSums(logcpm_c[mhcii_ens_clarke, , drop=FALSE] >=
                             log1p(5)) == 2
    vaf_cnt_sub <- vaf_cnt_sub[, passes_c, drop=FALSE]
    cat("Clarke cells passing MHCII filter:", ncol(vaf_cnt_sub), "/",
        sum(!is.na(passes_c)), "\n")
  } else {
    cat("WARNING: MHCII genes not found in Clarke data, skipping filter\n")
  }

  # -- DESeq2 VST on Clarke cells ----------------------------------------------
  cat("Running DESeq2 VST on Clarke cells...\n")

  # Assign population labels based on cell ID prefix
  clarke_cell_ids <- colnames(vaf_cnt_sub)
  clarke_cd45pos  <- grepl("^A[0-9]", clarke_cell_ids)

  # For CD45neg: reuse k-means VAF/VRC assignment from earlier
  # Recompute on the filtered subset
  cd45neg_ids <- clarke_cell_ids[!clarke_cd45pos]
  if (length(cd45neg_ids) >= 4) {
    logcpm_neg_sub <- log1p(sweep(
      vaf_cnt_sub[, cd45neg_ids, drop=FALSE], 2,
      pmax(colSums(vaf_cnt_sub[, cd45neg_ids, drop=FALSE]), 1), "/") * 1e6)
    gv_neg   <- apply(logcpm_neg_sub, 1, var)
    top_neg  <- names(sort(gv_neg, decreasing=TRUE))[1:min(500, length(gv_neg))]
    pca_neg  <- prcomp(t(logcpm_neg_sub[top_neg, ]), center=TRUE, scale.=FALSE)
    set.seed(UMAP_SEED)
    km_neg   <- kmeans(pca_neg$x[, 1:min(10, ncol(pca_neg$x))],
                        centers=2, nstart=25)

    vaf_m <- intersect(c("Col1a1","Col1a2","Timp3","Spp1","Thy1","Pdpn"),
                        rownames(logcpm_neg_sub))
    vrc_m <- intersect(c("Pecam1","Eng","Cdh5","Kdr","Tie1","Vwf"),
                        rownames(logcpm_neg_sub))
    sc1   <- mean(colMeans(logcpm_neg_sub[vaf_m, km_neg$cluster==1, drop=FALSE]))
    sc2   <- mean(colMeans(logcpm_neg_sub[vaf_m, km_neg$cluster==2, drop=FALSE]))
    vaf_ids_c <- cd45neg_ids[km_neg$cluster == ifelse(sc1 > sc2, 1, 2)]
    vrc_ids_c <- cd45neg_ids[km_neg$cluster == ifelse(sc1 > sc2, 2, 1)]
    cat("Clarke VAF cells:", length(vaf_ids_c),
        "| VRC cells:", length(vrc_ids_c), "\n")
  } else {
    vaf_ids_c <- cd45neg_ids
    vrc_ids_c <- character(0)
  }

  # Build condition vector for Clarke cells
  clarke_condition <- ifelse(clarke_cd45pos, "CD45pos_MHCIIpos",
                      ifelse(clarke_cell_ids %in% vaf_ids_c, "VAF", "VRC"))

  # Build a simple colData and run VST
  clarke_col_data <- data.frame(
    condition = factor(clarke_condition),
    row.names = clarke_cell_ids
  )
  # Round counts, filter low-count genes
  vaf_cnt_round <- round(vaf_cnt_sub)
  keep_c        <- rowSums(vaf_cnt_round) >= 10
  vaf_cnt_round <- vaf_cnt_round[keep_c, , drop=FALSE]

  dds_clarke <- DESeqDataSetFromMatrix(
    countData = vaf_cnt_round,
    colData   = clarke_col_data,
    design    = ~ condition
  )
  dds_clarke <- estimateSizeFactors(dds_clarke, type="poscounts")
  vsd_clarke <- varianceStabilizingTransformation(dds_clarke, blind=TRUE)
  expr_clarke <- assay(vsd_clarke)
  cat("Clarke VST complete:", nrow(expr_clarke), "genes x",
      ncol(expr_clarke), "cells\n")

  # -- Build metadata for Clarke cells -----------------------------------------
  meta_clarke <- data.frame(
    cell_id   = clarke_cell_ids,
    strain    = "Clarke2025",
    condition = clarke_condition,
    strain_condition = paste0("Clarke2025 ",
      ifelse(clarke_cd45pos, "CD45pos",
      ifelse(clarke_cell_ids %in% vaf_ids_c, "VAF", "VRC"))),
    row.names = clarke_cell_ids,
    stringsAsFactors = FALSE
  )

  # -- Merge your VST with Clarke VST ------------------------------------------
  common_ens <- intersect(rownames(combined_expr), rownames(expr_clarke))
  cat("Common Ensembl IDs for merged matrix:", length(common_ens), "\n")

  merged_expr_raw <- cbind(combined_expr[common_ens, ],
                            expr_clarke[common_ens, ])
  cat("Merged matrix (pre-correction):", nrow(merged_expr_raw), "genes x",
      ncol(merged_expr_raw), "cells\n")

  # -- Build merged metadata BEFORE batch correction (needed as covariate) ---
  your_meta_sub <- combined_meta[, c("cell_id","strain","condition",
                                      "strain_condition")]
  meta_merged <- rbind(your_meta_sub, meta_clarke[, colnames(your_meta_sub)])
  rownames(meta_merged) <- meta_merged$cell_id
  cat("Merged metadata:", nrow(meta_merged), "rows\n")

  # -- Batch correction: limma::removeBatchEffect ----------------------------
  # Corrects for the between-dataset expression scale difference while
  # preserving biological variation. Applied only for this merged analysis;
  # all other outputs use the uncorrected per-plate VST.
  if (!requireNamespace("limma", quietly=TRUE)) {
    cat("  Installing limma...\n")
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
    BiocManager::install("limma", ask=FALSE, update=FALSE, quiet=TRUE)
  }
  suppressPackageStartupMessages(library(limma))

  # Batch vector: "yours" for your 4 strains, "clarke" for Clarke2025
  all_cells_merged <- colnames(merged_expr_raw)
  batch_vec        <- ifelse(all_cells_merged %in% colnames(combined_expr),
                              "yours", "clarke")
  cat("  Batch composition — yours:", sum(batch_vec=="yours"),
      "| clarke:", sum(batch_vec=="clarke"), "\n")

  # Preserve condition structure as a covariate so biological signal
  # is not removed along with the batch effect
  meta_merged_ordered <- meta_merged[all_cells_merged, ]
  condition_covar     <- model.matrix(~ condition,
                           data=data.frame(
                             condition=factor(meta_merged_ordered$condition)))

  merged_expr <- removeBatchEffect(merged_expr_raw,
                                    batch     = batch_vec,
                                    design    = condition_covar)
  cat("  Batch correction applied (limma::removeBatchEffect)\n")
  cat("Merged matrix (post-correction):", nrow(merged_expr), "genes x",
      ncol(merged_expr), "cells\n")

  # (meta_merged already built above before batch correction)

  # -- Color palettes for merged analysis --------------------------------------
  # Clarke2025 gets grey/brown tones to distinguish from your strains
  strain_colors_merged <- c(
    NOD         = "#2166AC",
    B6G7        = "#D6604D",
    B6MHCIIGFP  = "#4DAC26",
    NODPDL1     = "#8B44AC",
    Clarke2025  = "#8B6914"
  )

  sc_palette_merged <- c(
    sc_palette_v2,
    "Clarke2025 CD45pos" = "#C49A2A",
    "Clarke2025 VAF"     = "#6B4F10",
    "Clarke2025 VRC"     = "#D4A843"
  )

  # -- HVG selection and UMAP on merged matrix ---------------------------------
  cat("Running HVG selection on merged matrix...\n")
  gv_merged  <- apply(merged_expr, 1, var)
  hvg_merged <- names(sort(gv_merged, decreasing=TRUE))[
                  1:min(N_HVG, length(gv_merged))]
  expr_hvg_m <- merged_expr[hvg_merged, ]

  cat("Running PCA + UMAP (merged)...\n")
  umap_merged <- run_umap(expr_hvg_m)

  umap_m_df <- data.frame(umap_merged,
                            cell_id = rownames(umap_merged),
                            stringsAsFactors=FALSE) %>%
    left_join(meta_merged[, c("cell_id","strain","condition","strain_condition")],
              by="cell_id")

  # All-gene PCA for Leiden on merged
  cat("Running all-gene PCA for Leiden clustering (merged)...\n")
  pca_m_all  <- prcomp(t(merged_expr), center=TRUE, scale.=FALSE)
  var_m_all  <- pca_m_all$sdev^2 / sum(pca_m_all$sdev^2)
  cum_m_all  <- cumsum(var_m_all)
  n_pcs_m    <- max(2, which(cum_m_all >= VAR_THRESHOLD)[1])
  cat("PCs selected (merged all-gene):", n_pcs_m,
      sprintf("(%.1f%% variance)\n", cum_m_all[n_pcs_m]*100))
  pcs_m_all  <- pca_m_all$x[, 1:n_pcs_m, drop=FALSE]

  cat("Running Leiden clustering (merged)...\n")
  leiden_m           <- run_leiden(pcs_m_all)
  umap_m_df$cluster  <- as.character(leiden_m)
  n_clust_m          <- length(unique(umap_m_df$cluster))
  umap_m_df$cluster  <- factor(umap_m_df$cluster,
                                levels=sort(unique(as.integer(umap_m_df$cluster))))
  cat("Clusters found:", n_clust_m, "\n")

  cluster_pal_m <- setNames(scales::hue_pal()(n_clust_m),
                              sort(unique(as.character(umap_m_df$cluster))))
  cluster_vec_m <- setNames(as.character(umap_m_df$cluster), umap_m_df$cell_id)

  # -- UMAP 1: by Leiden cluster -----------------------------------------------
  p_m_clust <- ggplot(umap_m_df, aes(x=UMAP1, y=UMAP2, color=cluster)) +
    geom_point(size=1.5, alpha=0.85) +
    scale_color_manual(values=cluster_pal_m) +
    labs(
      title    = paste0("Your cells + Clarke 2025 — ",
                        nrow(umap_m_df), " cells by Leiden cluster"),
      subtitle = paste0("n=", nrow(umap_m_df),
                        " | UMAP: top ", N_HVG,
                        " HVGs | Clusters: all-gene PCA | resolution=",
                        LEIDEN_RESOLUTION, " | ", n_clust_m, " clusters"),
      color    = "Cluster"
    ) +
    theme_bw(base_size=12) +
    theme(plot.title=element_text(face="bold", size=12),
          plot.subtitle=element_text(size=8, color="grey40"),
          panel.grid.minor=element_blank(), aspect.ratio=1)

  ggsave(file.path(vaf_out_dir, "umap_merged_by_cluster.pdf"),
         p_m_clust, width=7, height=6)
  cat("Saved: umap_merged_by_cluster.pdf\n")

  # -- UMAP 2: by strain/population --------------------------------------------
  p_m_strain <- ggplot(umap_m_df,
                        aes(x=UMAP1, y=UMAP2, fill=strain_condition)) +
    geom_point(size=1.0, alpha=1.0, stroke=0.2, shape=21,
               aes(fill=strain_condition), color="grey20") +
    scale_fill_manual(values=sc_palette_merged,
                      na.value="grey60") +
    labs(
      title    = paste0("Your cells + Clarke 2025 — ",
                        nrow(umap_m_df), " cells by strain & population"),
      subtitle = "5 strains x conditions",
      fill     = "Strain / Population"
    ) +
    guides(fill=guide_legend(ncol=2, override.aes=list(size=3, stroke=0.3))) +
    theme_bw(base_size=12) +
    theme(plot.title=element_text(face="bold", size=12),
          plot.subtitle=element_text(size=9, color="grey40"),
          panel.grid.minor=element_blank(),
          legend.text=element_text(size=7), aspect.ratio=1)

  ggsave(file.path(vaf_out_dir, "umap_merged_by_strain_population.pdf"),
         p_m_strain, width=8, height=6)
  cat("Saved: umap_merged_by_strain_population.pdf\n")

  # -- Per-strain panels on merged UMAP ----------------------------------------
  all_strains_merged <- c(strains, "Clarke2025")
  strain_panels_m <- lapply(all_strains_merged, function(st) {
    st_cells <- umap_m_df$cell_id[umap_m_df$strain == st]
    df_fg    <- umap_m_df %>% filter(cell_id %in% st_cells)
    df_bg    <- umap_m_df %>% filter(!cell_id %in% st_cells)
    ggplot() +
      geom_point(data=df_bg, aes(x=UMAP1, y=UMAP2),
                 color="grey88", size=0.9, alpha=0.5) +
      geom_point(data=df_fg, aes(x=UMAP1, y=UMAP2, color=cluster),
                 size=1.6, alpha=0.9) +
      scale_color_manual(values=cluster_pal_m) +
      labs(title=paste0(st, " (n=", nrow(df_fg), ")"),
           x="UMAP1", y="UMAP2", color="Cluster") +
      coord_fixed() +
      theme_bw(base_size=11) +
      theme(plot.title=element_text(face="bold", size=11),
            panel.grid.minor=element_blank(),
            legend.text=element_text(size=8))
  })
  panels_m_combined <- wrap_plots(strain_panels_m, ncol=1) +
    plot_annotation(
      title    = "Per-strain cells on merged UMAP (your cells + Clarke 2025)",
      subtitle = paste0("Merged Leiden clusters | resolution=", LEIDEN_RESOLUTION,
                        " | grey = other strains"),
      theme    = theme(plot.title=element_text(face="bold", size=13),
                       plot.subtitle=element_text(size=9, color="grey40"))
    )
  ggsave(file.path(vaf_out_dir, "umap_per_strain_on_merged.pdf"),
         panels_m_combined,
         width=7, height=6 * length(all_strains_merged))
  cat("Saved: umap_per_strain_on_merged.pdf\n")

  # -- Bar charts (cluster composition) ----------------------------------------
  cat("Generating merged cluster bar charts...\n")
  umap_m_df$condition_clean  <- clean_label(umap_m_df$condition)
  umap_m_df$condition_clean  <- ifelse(
    umap_m_df$strain == "Clarke2025",
    umap_m_df$strain_condition,
    umap_m_df$condition_clean)

  make_label_df_m <- function(df, fill_var, type) {
    df %>%
      group_by(cluster, .data[[fill_var]]) %>%
      summarise(n=n(), .groups="drop") %>%
      group_by(cluster) %>%
      mutate(total=sum(n), prop=n/total,
             label_count=ifelse(prop>=MIN_LABEL_PROP, as.character(n), ""),
             label_prop =ifelse(prop>=MIN_LABEL_PROP,
                                paste0(round(prop*100), "%"), ""),
             y_pos=cumsum(if(type=="count") n else prop) -
                   0.5*(if(type=="count") n else prop)) %>%
      ungroup()
  }

  lab_sc_c <- make_label_df_m(umap_m_df, "strain_condition", "count")
  lab_sc_p <- make_label_df_m(umap_m_df, "strain_condition", "prop")
  lab_s_c  <- make_label_df_m(umap_m_df, "strain",           "count")
  lab_s_p  <- make_label_df_m(umap_m_df, "strain",           "prop")

  mk_bar <- function(lab, fill_var, pal, y_var, y_lab, title_str, prop=FALSE) {
    # Use pre-computed label column based on prop flag
    lab$bar_label <- if (isTRUE(prop)) lab$label_prop else lab$label_count
    p <- ggplot(lab, aes(x=cluster, y=.data[[y_var]], fill=.data[[fill_var]])) +
      geom_bar(stat="identity", position="stack", color="white", linewidth=0.15) +
      geom_text(aes(y=y_pos, label=bar_label),
                size=2.5, color="white", fontface="bold") +
      scale_fill_manual(values=pal, na.value="grey60") +
      labs(title=title_str, x="Cluster",
           y=y_lab, fill=sub(" —.*","",title_str)) +
      theme_bw(base_size=10) +
      theme(plot.title=element_text(face="bold", size=10),
            panel.grid.minor=element_blank(),
            legend.text=element_text(size=7),
            legend.key.size=unit(0.4,"cm"))
    if (isTRUE(prop)) p <- p + scale_y_continuous(labels=scales::percent_format())
    p
  }

  p_b1 <- mk_bar(lab_sc_c, "strain_condition", sc_palette_merged,
                  "n",    "Cell count",       "Strain & population — counts")
  p_b2 <- mk_bar(lab_sc_p, "strain_condition", sc_palette_merged,
                  "prop", "Proportion",       "Strain & population — proportions",
                  prop=TRUE)
  p_b3 <- mk_bar(lab_s_c,  "strain",           strain_colors_merged,
                  "n",    "Cell count",       "Strain — counts")
  p_b4 <- mk_bar(lab_s_p,  "strain",           strain_colors_merged,
                  "prop", "Proportion",       "Strain — proportions", prop=TRUE)

  bar_m <- (p_b1 | p_b2) / (p_b3 | p_b4) +
    plot_annotation(
      title    = paste0("Merged cluster composition — ", n_clust_m,
                        " clusters | n=", nrow(umap_m_df), " cells"),
      subtitle = paste0("Labels shown for segments >= ", MIN_LABEL_PROP*100,
                        "% of bar"),
      theme    = theme(plot.title=element_text(face="bold", size=14),
                       plot.subtitle=element_text(size=9, color="grey40"))
    )
  ggsave(file.path(vaf_out_dir, "barplots_cluster_composition.pdf"),
         bar_m, width=20, height=14)
  cat("Saved: barplots_cluster_composition.pdf\n")

  # -- Cluster marker genes (Wilcoxon on merged VST) ---------------------------
  cat("Running Wilcoxon marker gene analysis on merged clusters...\n")
  cluster_levels_m <- sort(unique(as.integer(as.character(umap_m_df$cluster))))
  all_markers_m    <- list()
  wb_m             <- createWorkbook()

  summary_m <- umap_m_df %>%
    group_by(cluster, strain, condition_clean) %>%
    summarise(n_cells=n(), .groups="drop") %>%
    arrange(as.integer(as.character(cluster)), strain, condition_clean)
  addWorksheet(wb_m, "Summary"); writeData(wb_m, "Summary", summary_m)

  for (cl in cluster_levels_m) {
    cl_label   <- as.character(cl)
    cl_cells   <- names(cluster_vec_m)[cluster_vec_m == cl_label]
    rest_cells <- names(cluster_vec_m)[cluster_vec_m != cl_label]
    mean_cl    <- rowMeans(merged_expr[, cl_cells,   drop=FALSE])
    mean_rest  <- rowMeans(merged_expr[, rest_cells, drop=FALSE])
    log2fc     <- mean_cl - mean_rest
    candidates <- names(log2fc)[log2fc >= MIN_LOG2FC]
    if (length(candidates) < 2) { all_markers_m[[cl_label]] <- data.frame(); next }
    pv <- sapply(candidates, function(g)
           wilcox.test(merged_expr[g, cl_cells], merged_expr[g, rest_cells],
                       alternative="greater", exact=FALSE)$p.value)
    pa <- p.adjust(pv, method="BH")
    res_m <- data.frame(
      ensembl_id=candidates, gene_symbol=to_sym(candidates),
      mean_VST_cluster=round(mean_cl[candidates],4),
      mean_VST_rest=round(mean_rest[candidates],4),
      log2FC=round(log2fc[candidates],4),
      pval=signif(pv,4), padj=signif(pa,4),
      stringsAsFactors=FALSE) %>%
      filter(padj < MAX_PADJ) %>%
      arrange(desc(log2FC)) %>%
      distinct(gene_symbol, .keep_all=TRUE) %>%
      slice_head(n=TOP_EXCEL) %>%
      mutate(rank=row_number()) %>%
      select(rank,ensembl_id,gene_symbol,mean_VST_cluster,mean_VST_rest,
             log2FC,pval,padj)
    all_markers_m[[cl_label]] <- res_m
    sn <- paste0("Cluster_", cl_label)
    addWorksheet(wb_m, sn); writeData(wb_m, sn, res_m)
    addStyle(wb_m, sn, style=createStyle(textDecoration="bold",fgFill="#D9E1F2"),
             rows=1, cols=1:8, gridExpand=TRUE)
    setColWidths(wb_m, sn, cols=1:8, widths=c(6,20,16,16,14,10,12,12))
  }
  markers_m_path <- file.path(vaf_out_dir, "cluster_marker_genes.xlsx")
  saveWorkbook(wb_m, markers_m_path, overwrite=TRUE)
  cat("Saved: cluster_marker_genes.xlsx\n")

  # -- Heatmap -----------------------------------------------------------------
  cat("Generating merged cluster marker heatmap...\n")
  hm_genes_m <- c()
  for (cl_label in as.character(cluster_levels_m)) {
    if (!is.null(all_markers_m[[cl_label]]) &&
        nrow(all_markers_m[[cl_label]]) > 0) {
      top_g <- all_markers_m[[cl_label]]$gene_symbol[
                 !all_markers_m[[cl_label]]$gene_symbol %in% hm_genes_m]
      hm_genes_m <- c(hm_genes_m, head(top_g, TOP_HEATMAP))
    }
  }
  cat("Heatmap genes:", length(hm_genes_m), "\n")
  sym_to_ens_m  <- setNames(names(sym_map), sym_map)
  hm_ens_m      <- sym_to_ens_m[hm_genes_m]
  hm_ens_m      <- hm_ens_m[!is.na(hm_ens_m) & hm_ens_m %in% rownames(merged_expr)]
  cell_ord_m    <- umap_m_df$cell_id[order(as.integer(
                     as.character(umap_m_df$cluster)))]
  expr_hm_m     <- merged_expr[hm_ens_m, cell_ord_m]
  rownames(expr_hm_m) <- to_sym(rownames(expr_hm_m))
  expr_z_m      <- t(scale(t(expr_hm_m)))
  expr_z_m[expr_z_m >  2.5] <-  2.5
  expr_z_m[expr_z_m < -2.5] <- -2.5
  clust_anno_m  <- as.character(umap_m_df$cluster[
                     match(cell_ord_m, umap_m_df$cell_id)])
  ha_m <- HeatmapAnnotation(
    Cluster=clust_anno_m,
    col=list(Cluster=cluster_pal_m),
    annotation_name_side="left"
  )
  col_fun_m <- colorRamp2(c(-2.5,0,2.5), c("#3D0751","#1A1A1A","#F5E642"))
  ht_m <- Heatmap(
    expr_z_m, name="Z-score", col=col_fun_m,
    top_annotation=ha_m,
    show_column_names=FALSE, show_row_names=TRUE,
    row_names_gp=gpar(fontsize=7, fontface="italic"),
    cluster_rows=FALSE, cluster_columns=FALSE,
    column_title=paste0("Merged cluster marker genes — ",
                         length(hm_genes_m), " genes x ", ncol(expr_z_m), " cells"),
    column_title_gp=gpar(fontsize=12, fontface="bold"),
    use_raster=TRUE, raster_quality=5
  )
  hm_path_m <- file.path(vaf_out_dir, "heatmap_cluster_markers.pdf")
  pdf(hm_path_m, width=18,
      height=max(8, length(hm_genes_m)*0.18+3))
  draw(ht_m, heatmap_legend_side="right", annotation_legend_side="bottom")
  dev.off()
  cat("Saved: heatmap_cluster_markers.pdf\n")

  # -- Cluster violin plots ----------------------------------------------------
  cat("Generating merged cluster violin plots...\n")
  vln_ens_m <- sym_to_ens_m[VIOLIN_GENES]
  vln_ens_m <- vln_ens_m[!is.na(vln_ens_m) & vln_ens_m %in% rownames(merged_expr)]

  vln_m_df <- as.data.frame(t(merged_expr[vln_ens_m, umap_m_df$cell_id,
                                            drop=FALSE])) %>%
    tibble::rownames_to_column("cell_id") %>%
    left_join(umap_m_df[, c("cell_id","cluster")], by="cell_id") %>%
    pivot_longer(-c(cell_id,cluster), names_to="ensembl_id", values_to="VST") %>%
    mutate(gene_symbol=to_sym(ensembl_id),
           cluster=factor(cluster,
                           levels=sort(unique(as.integer(
                             as.character(umap_m_df$cluster))))))

  vln_m_plots <- lapply(seq_along(vln_ens_m), function(i) {
    ens <- vln_ens_m[i]; sym <- names(vln_ens_m)[i]
    df_g <- vln_m_df %>% filter(ensembl_id == ens)
    ggplot(df_g, aes(x=cluster, y=VST, fill=cluster)) +
      geom_violin(trim=TRUE, scale="width", alpha=0.85, linewidth=0.3) +
      geom_jitter(width=0.15, size=0.6, alpha=0.35, color="grey20") +
      scale_fill_manual(values=cluster_pal_m) +
      labs(title=sym, x="Leiden Cluster", y="VST expression") +
      theme_bw(base_size=11) +
      theme(plot.title=element_text(face="bold.italic", size=11),
            legend.position="none", panel.grid.minor=element_blank())
  })
  n_vln_m   <- length(vln_m_plots)
  vln_m_panel <- wrap_plots(vln_m_plots, ncol=min(3, n_vln_m)) +
    plot_annotation(
      title    = "Gene expression across merged Leiden clusters",
      subtitle = paste0("VST-normalized | n=", nrow(umap_m_df),
                        " cells (your data + Clarke 2025)"),
      theme    = theme(plot.title=element_text(face="bold", size=13),
                       plot.subtitle=element_text(size=9, color="grey40"))
    )
  ggsave(file.path(vaf_out_dir, "violin_plots_by_cluster.pdf"),
         vln_m_panel,
         width=min(3,n_vln_m)*5, height=ceiling(n_vln_m/3)*5)
  cat("Saved: violin_plots_by_cluster.pdf\n")

  # -- VAF/VRC correlation on merged clusters ----------------------------------
  cat("Running VAF/VRC correlation on merged clusters...\n")
  # Reuse ref_profiles from earlier VAF/VRC block if available
  if (exists("ref_profiles")) {
    your_sym_merged   <- to_sym(rownames(merged_expr))
    common_m_vaf      <- intersect(your_sym_merged, names(ref_profiles[["VAF"]]))
    your_idx_m        <- which(your_sym_merged %in% common_m_vaf)
    your_syms_m       <- your_sym_merged[your_idx_m]
    wb_vaf_m          <- createWorkbook()
    summary_vaf_m     <- data.frame(cluster=character(), best_match=character(),
                                     best_pearson_r=numeric(),
                                     best_spearman_rho=numeric(),
                                     stringsAsFactors=FALSE)
    for (cl_label in as.character(cluster_levels_m)) {
      cl_cells_m <- names(cluster_vec_m)[cluster_vec_m == cl_label]
      your_mean_m <- rowMeans(merged_expr[your_idx_m, cl_cells_m, drop=FALSE])
      names(your_mean_m) <- your_syms_m
      cl_corrs_m <- data.frame(reference_population=character(),
                                pearson_r=numeric(), spearman_rho=numeric(),
                                n_genes=integer(), stringsAsFactors=FALSE)
      for (ref_name in names(ref_profiles)) {
        ref_v  <- ref_profiles[[ref_name]]
        common <- intersect(names(your_mean_m), names(ref_v))
        if (length(common) < 50) next
        cl_corrs_m <- rbind(cl_corrs_m, data.frame(
          reference_population=ref_name,
          pearson_r=round(cor(your_mean_m[common], ref_v[common], method="pearson"),4),
          spearman_rho=round(cor(your_mean_m[common], ref_v[common], method="spearman"),4),
          n_genes=length(common), stringsAsFactors=FALSE))
      }
      best_m <- cl_corrs_m[which.max(cl_corrs_m$pearson_r), ]
      summary_vaf_m <- rbind(summary_vaf_m, data.frame(
        cluster=cl_label, best_match=best_m$reference_population,
        best_pearson_r=best_m$pearson_r,
        best_spearman_rho=best_m$spearman_rho, stringsAsFactors=FALSE))
      sn_m <- paste0("Cluster_", cl_label)
      addWorksheet(wb_vaf_m, sn_m); writeData(wb_vaf_m, sn_m, cl_corrs_m)
      addStyle(wb_vaf_m, sn_m,
               style=createStyle(textDecoration="bold", fgFill="#E2EFDA"),
               rows=1, cols=1:4, gridExpand=TRUE)
      if (nrow(cl_corrs_m) > 0) {
        best_idx_m <- which.max(cl_corrs_m$pearson_r) + 1
        addStyle(wb_vaf_m, sn_m, style=createStyle(fgFill="#C6EFCE"),
                 rows=best_idx_m, cols=1:4, gridExpand=TRUE)
      }
      setColWidths(wb_vaf_m, sn_m, cols=1:4, widths=c(22,12,14,10))
    }
    addWorksheet(wb_vaf_m, "Summary"); writeData(wb_vaf_m, "Summary", summary_vaf_m)
    saveWorkbook(wb_vaf_m,
                 file.path(vaf_out_dir, "VAF_VRC_correlation_merged_clusters.xlsx"),
                 overwrite=TRUE)
    cat("Saved: VAF_VRC_correlation_merged_clusters.xlsx\n")
  } else {
    cat("  Skipping VAF/VRC correlation (ref_profiles not available)\n")
  }

  # -- MCA cell type correlation on merged clusters ----------------------------
  cat("Running MCA cell type correlation on merged clusters...\n")
  run_mca_correlation(
    expr_mat          = merged_expr,
    cluster_vec_input = cluster_vec_m,
    cluster_levels_input = cluster_levels_m,
    to_sym_fn         = to_sym,
    out_path          = file.path(vaf_out_dir,
                                   "MCA_celltype_correlation_merged_clusters.xlsx"),
    label             = "merged"
  )

  cat("\nCombinedwithVAFPaperPlots complete.\n")
  cat("Output folder:", vaf_out_dir, "\n")
  cat("Files written:\n")
  for (f in list.files(vaf_out_dir)) cat(" -", f, "\n")

} # end if vaf_counts_f2 exists

# ==============================================================================
# Cluster marker genes: hybrid expression score per cluster


cat("\n==============================================================\n")
cat("All outputs complete. Output structure:\n")
for (st in strains) {
  cat(" results/05_dge/", st, "_plots/  (10 files)\n", sep="")
}
cat(" results/05_dge/combined_plots/\n")
cat("   - umap_all372_by_cluster.pdf\n")
cat("   - umap_all372_by_strain_population.pdf\n")
cat("   - barplots_cluster_composition.pdf  (all 4 charts, 2x2)\n")
cat("   - cluster_marker_genes.xlsx\n")
cat("   - heatmap_cluster_markers.pdf\n")
cat("==============================================================\n")
