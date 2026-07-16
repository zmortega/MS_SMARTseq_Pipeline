#!/usr/bin/env python3
"""
pipeline.py — SMARTseq HT scRNA-seq pipeline orchestrator
Plate-based full-length scRNA-seq | Mouse GRCm39 | Apple Silicon Mac

Usage:
    python pipeline.py --config config.yaml
    python pipeline.py --config config.yaml --steps qc,align,count,dge
    python pipeline.py --config config.yaml --resume          # skip completed steps
    python pipeline.py --config config.yaml --cells A01,A02   # subset of cells
"""

import argparse
import logging
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import yaml

# ── Logging setup ──────────────────────────────────────────────────────────────
def setup_logging(log_dir: Path) -> logging.Logger:
    log_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"pipeline_{ts}.log"

    fmt = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(
        level=logging.INFO,
        format=fmt,
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout),
        ],
    )
    log = logging.getLogger("smartseq_pipeline")
    log.info(f"Log file: {log_file}")
    return log


# ── Config loader ──────────────────────────────────────────────────────────────
def load_config(path: str) -> dict:
    with open(path) as f:
        cfg = yaml.safe_load(f)
    return cfg


# ── Cell discovery ─────────────────────────────────────────────────────────────
def discover_cells(fastq_dir: Path, r1_suffix: str, r2_suffix: str, log) -> list[str]:
    r1_files = sorted(fastq_dir.glob(f"*{r1_suffix}"))
    if not r1_files:
        log.error(f"No files matching *{r1_suffix} found in {fastq_dir}")
        sys.exit(1)

    cells = []
    for r1 in r1_files:
        cell_id = r1.name[: -len(r1_suffix)]
        r2 = fastq_dir / f"{cell_id}{r2_suffix}"
        if not r2.exists():
            log.warning(f"  Missing R2 for cell {cell_id}, skipping.")
            continue
        cells.append(cell_id)

    log.info(f"Discovered {len(cells)} cell(s) with matched R1/R2 pairs.")
    return cells


# ── Shell runner ───────────────────────────────────────────────────────────────
def run(cmd: str, log, check=True) -> subprocess.CompletedProcess:
    log.debug(f"CMD: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.stdout:
        log.debug(result.stdout.strip())
    if result.stderr:
        log.debug(result.stderr.strip())
    if check and result.returncode != 0:
        log.error(f"Command failed (exit {result.returncode}):\n{cmd}")
        log.error(result.stderr[-2000:])
        sys.exit(result.returncode)
    return result


def done_flag(path: Path) -> bool:
    return (path / ".done").exists()


def mark_done(path: Path):
    path.mkdir(parents=True, exist_ok=True)
    (path / ".done").touch()


# ── STEP 1: FastQC ─────────────────────────────────────────────────────────────
def step_fastqc(cells, cfg, log, resume):
    fastq_dir = Path(cfg["FASTQ_DIR"])
    qc_dir = Path(cfg["OUTPUT_DIR"]) / "01_fastqc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    log.info(f"── STEP 1: FastQC ({len(cells)} cells) ──")

    for cell in cells:
        cell_qc = qc_dir / cell
        if resume and done_flag(cell_qc):
            log.info(f"  [skip] {cell}")
            continue
        cell_qc.mkdir(parents=True, exist_ok=True)
        r1 = fastq_dir / f"{cell}{cfg['R1_SUFFIX']}"
        r2 = fastq_dir / f"{cell}{cfg['R2_SUFFIX']}"
        cmd = f"fastqc --threads 2 --outdir {cell_qc} {r1} {r2}"
        run(cmd, log)
        mark_done(cell_qc)
        log.info(f"  [done] {cell}")

    # MultiQC summary
    multiqc_out = Path(cfg["OUTPUT_DIR"]) / "01_fastqc_multiqc"
    cmd = f"multiqc {qc_dir} --outdir {multiqc_out} --filename fastqc_report -f"
    run(cmd, log, check=False)  # non-fatal
    log.info(f"  MultiQC report: {multiqc_out}/fastqc_report.html")


# ── STEP 2: Trimmomatic ────────────────────────────────────────────────────────
def step_trim(cells, cfg, log, resume):
    fastq_dir = Path(cfg["FASTQ_DIR"])
    trim_dir = Path(cfg["OUTPUT_DIR"]) / "02_trimmed"
    trim_dir.mkdir(parents=True, exist_ok=True)
    log.info(f"── STEP 2: Trimmomatic ({len(cells)} cells) ──")

    # Find adapter file bundled with trimmomatic conda install
    adapter_search = run(
        "find $CONDA_PREFIX -name 'TruSeq3-PE.fa' 2>/dev/null | head -1",
        log, check=False
    )
    adapter_file = adapter_search.stdout.strip()
    if not adapter_file:
        log.warning("Could not auto-find TruSeq3-PE.fa — using ILLUMINACLIP fallback")
        adapter_file = "TruSeq3-PE.fa"

    for cell in cells:
        cell_trim = trim_dir / cell
        if resume and done_flag(cell_trim):
            log.info(f"  [skip] {cell}")
            continue
        cell_trim.mkdir(parents=True, exist_ok=True)

        r1 = fastq_dir / f"{cell}{cfg['R1_SUFFIX']}"
        r2 = fastq_dir / f"{cell}{cfg['R2_SUFFIX']}"
        r1p = cell_trim / f"{cell}_R1_paired.fastq.gz"
        r1u = cell_trim / f"{cell}_R1_unpaired.fastq.gz"
        r2p = cell_trim / f"{cell}_R2_paired.fastq.gz"
        r2u = cell_trim / f"{cell}_R2_unpaired.fastq.gz"
        trim_log = cell_trim / f"{cell}_trim.log"

        cmd = (
            f"trimmomatic PE -threads 4 -phred33 "
            f"{r1} {r2} "
            f"{r1p} {r1u} {r2p} {r2u} "
            f"ILLUMINACLIP:{adapter_file}:2:30:10:2:keepBothReads "
            f"LEADING:{cfg['TRIM_LEADING']} "
            f"TRAILING:{cfg['TRIM_TRAILING']} "
            f"SLIDINGWINDOW:{cfg['TRIM_SLIDINGWINDOW']} "
            f"MINLEN:{cfg['TRIM_MINLEN']} "
            f"2> {trim_log}"
        )
        run(cmd, log)
        mark_done(cell_trim)
        log.info(f"  [done] {cell}")


# ── STEP 3: HISAT2 alignment ───────────────────────────────────────────────────
def step_align(cells, cfg, log, resume):
    trim_dir  = Path(cfg["OUTPUT_DIR"]) / "02_trimmed"
    align_dir = Path(cfg["OUTPUT_DIR"]) / "03_aligned"
    align_dir.mkdir(parents=True, exist_ok=True)
    log.info(f"── STEP 3: HISAT2 alignment ({len(cells)} cells) ──")

    for cell in cells:
        cell_align = align_dir / cell
        if resume and done_flag(cell_align):
            log.info(f"  [skip] {cell}")
            continue
        cell_align.mkdir(parents=True, exist_ok=True)

        r1p_gz = trim_dir / cell / f"{cell}_R1_paired.fastq.gz"
        r2p_gz = trim_dir / cell / f"{cell}_R2_paired.fastq.gz"
        bam    = cell_align / f"{cell}_Aligned.sortedByCoord.out.bam"

        # HISAT2 can read gzipped FASTQs natively on Mac
        cmd = (
            f"hisat2 "
            f"-p {cfg['ALIGN_THREADS']} "
            f"--dta "
            f"-x {cfg['HISAT2_INDEX']} "
            f"-1 {r1p_gz} "
            f"-2 {r2p_gz} "
            f"--no-spliced-alignment "
            f"2>{cell_align}/{cell}_hisat2.log "
            f"| samtools sort -@ 4 -m 2G "
            f"-o {bam}"
        )
        run(cmd, log)

        # Verify BAM has reads
        check = subprocess.run(
            f"samtools view -c {bam}", shell=True, capture_output=True, text=True
        )
        read_count = int(check.stdout.strip()) if check.stdout.strip().isdigit() else 0
        if read_count == 0:
            log.error(f"  BAM for {cell} has 0 reads!")
            bam.unlink(missing_ok=True)
            (cell_align / ".done").unlink(missing_ok=True)
            continue

        run(f"samtools index {bam}", log)
        mark_done(cell_align)
        log.info(f"  [done] {cell} ({read_count:,} reads)")

    # MultiQC on HISAT2 logs
    multiqc_out = Path(cfg["OUTPUT_DIR"]) / "03_aligned_multiqc"
    cmd = f"multiqc {align_dir} --outdir {multiqc_out} --filename alignment_report -f"
    run(cmd, log, check=False)
    log.info(f"  MultiQC alignment report: {multiqc_out}/alignment_report.html")


# ── STEP 4: featureCounts ──────────────────────────────────────────────────────
def step_featurecounts(cells, cfg, log, resume):
    align_dir = Path(cfg["OUTPUT_DIR"]) / "03_aligned"
    count_dir = Path(cfg["OUTPUT_DIR"]) / "04_counts"
    count_dir.mkdir(parents=True, exist_ok=True)
    log.info(f"── STEP 4: featureCounts ──")

    flag_path = count_dir / "all_cells"
    if resume and done_flag(flag_path):
        log.info("  [skip] featureCounts already done")
        return

    # Collect all BAMs
    bams = []
    for cell in cells:
        bam = align_dir / cell / f"{cell}_Aligned.sortedByCoord.out.bam"
        if bam.exists():
            bams.append(str(bam))
        else:
            log.warning(f"  BAM not found for {cell}: {bam}")

    if not bams:
        log.error("No BAM files found. Did alignment complete successfully?")
        sys.exit(1)

    # Decompress GTF if needed
    gtf_gz = Path(cfg["GTF_FILE"])
    gtf = gtf_gz.with_suffix("") if gtf_gz.suffix == ".gz" else gtf_gz
    if not gtf.exists():
        log.info("  Decompressing GTF for featureCounts...")
        run(f"gunzip -k {gtf_gz}", log)

    # BAMs are paired-end aligned (hisat2 -1/-2), so newer featureCounts/subread
    # (>=2.0.3) requires -p just to declare that fact, or it hard-errors with
    # "Paired-end reads were detected in single-end read library". Passing -p
    # alone (without --countReadPairs) keeps the original per-read counting
    # behavior intended here (each mate counted individually, not as one
    # fragment) — see https://support.bioconductor.org/p/9159588/
    paired_flag = "-p"  # SMARTseq: count individual reads, not pairs
    bam_list = " ".join(bams)
    out_counts = count_dir / "raw_counts.txt"

    cmd = (
        f"featureCounts "
        f"-T {cfg['FC_THREADS']} "
        f"-a {gtf} "
        f"-o {out_counts} "
        f"-s {cfg['FC_STRAND']} "
        f"{paired_flag} "
        f"-Q {cfg['FC_MIN_MAPQ']} "
        f"--primary "
        f"-t exon "
        f"-g gene_id "
        f"{bam_list}"
    )
    run(cmd, log)

    # Clean up column names to just cell IDs
    _clean_featurecounts(out_counts, cells, log)
    mark_done(flag_path)
    log.info(f"  Count matrix: {out_counts}")


def _clean_featurecounts(counts_file: Path, cells: list, log):
    """Strip full BAM paths from featureCounts header, extract cell IDs from paths."""
    import re
    cleaned = counts_file.parent / "counts_clean.txt"
    with open(counts_file) as f_in, open(cleaned, "w") as f_out:
        for i, line in enumerate(f_in):
            if i == 0:          # comment line
                f_out.write(line)
                continue
            if i == 1:          # header line — extract cell ID from BAM path
                parts = line.rstrip("\n").split("\t")
                fixed_cols = []
                for col in parts[6:]:
                    # Extract cell ID from path like:
                    # results/03_aligned/B6G7_A1_S1_L001/B6G7_A1_S1_L001_Aligned.sortedByCoord.out.bam
                    # -> B6G7_A1_S1_L001
                    bam_name = col.split("/")[-1]  # get filename
                    cell_id = bam_name.replace("_Aligned.sortedByCoord.out.bam", "")
                    fixed_cols.append(cell_id)
                fixed = parts[:6] + fixed_cols
                f_out.write("\t".join(fixed) + "\n")
                continue
            f_out.write(line)
    log.info(f"  Cleaned count matrix: {cleaned}")


# ── STEP 5: DESeq2 DGE ────────────────────────────────────────────────────────
def step_dge(cfg, log, resume):
    count_dir = Path(cfg["OUTPUT_DIR"]) / "04_counts"
    dge_dir = Path(cfg["OUTPUT_DIR"]) / "05_dge"
    dge_dir.mkdir(parents=True, exist_ok=True)
    log.info("── STEP 5: DESeq2 differential gene expression ──")

    flag_path = dge_dir / "deseq2_done"
    if resume and done_flag(flag_path):
        log.info("  [skip] DESeq2 already done")
        return

    # Write the R script dynamically with config values injected
    r_script = _build_deseq2_script(cfg, count_dir, dge_dir)
    r_script_path = dge_dir / "run_deseq2.R"
    with open(r_script_path, "w") as f:
        f.write(r_script)

    log.info("  Running DESeq2 R script...")
    cmd = f"Rscript {r_script_path} 2>&1 | tee {dge_dir}/deseq2.log"
    run(cmd, log)
    mark_done(flag_path)
    log.info(f"  DGE results: {dge_dir}/")


def _build_deseq2_script(cfg, count_dir, dge_dir) -> str:
    return f"""
# DESeq2 differential expression analysis
# Auto-generated by pipeline.py

suppressPackageStartupMessages({{
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(openxlsx)
}})

cat("=== DESeq2 Analysis ===\\n")

# ── Load count matrix ─────────────────────────────────────────────────────────
counts_file <- "{count_dir}/counts_clean.txt"
cat("Loading counts from:", counts_file, "\\n")

raw <- read.table(counts_file, header=TRUE, skip=1, row.names=1, sep="\\t",
                  check.names=FALSE)
# Keep only expression columns (drop Chr, Start, End, Strand, Length)
expr_cols <- 6:ncol(raw)
count_mat  <- as.matrix(raw[, expr_cols])
storage.mode(count_mat) <- "integer"
cat("Genes:", nrow(count_mat), "| Cells:", ncol(count_mat), "\\n")

# ── Load metadata ─────────────────────────────────────────────────────────────
meta_file <- "{cfg['METADATA_FILE']}"
cat("Loading metadata from:", meta_file, "\\n")
meta <- read.csv(meta_file, stringsAsFactors=FALSE)
rownames(meta) <- meta$cell_id

# Reorder to match count matrix columns
common_cells <- intersect(colnames(count_mat), rownames(meta))
if (length(common_cells) == 0) {{
  stop("No matching cell IDs between count matrix and metadata. \\n",
       "Count matrix cells (first 5): ", paste(head(colnames(count_mat),5), collapse=", "), "\\n",
       "Metadata cell IDs (first 5): ", paste(head(rownames(meta),5), collapse=", "))
}}
cat("Cells with metadata:", length(common_cells), "\\n")
count_mat <- count_mat[, common_cells]
meta      <- meta[common_cells, ]
meta$condition <- relevel(factor(meta$condition), ref="{cfg['REFERENCE_LEVEL']}")

# ── Filter low-count genes ────────────────────────────────────────────────────
keep <- rowSums(count_mat) >= {cfg['MIN_COUNTS']}
cat("Genes after count filter (>={cfg['MIN_COUNTS']} total):", sum(keep), "\\n")
count_mat <- count_mat[keep, ]

# ── DESeq2 ────────────────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = meta,
  design    = {cfg['DESIGN_FORMULA']}
)
dds <- DESeq(dds, test="Wald", fitType="parametric")
cat("DESeq2 complete.\\n")

# ── Results ───────────────────────────────────────────────────────────────────
res     <- results(dds, alpha={cfg['PADJ_THRESHOLD']})
res_df  <- as.data.frame(res) %>%
             rownames_to_column("gene_id") %>%
             arrange(padj)

res_sig <- res_df %>% filter(padj < {cfg['PADJ_THRESHOLD']}, abs(log2FoldChange) >= {cfg['LFC_THRESHOLD']})
cat("Significant DEGs (padj<{cfg['PADJ_THRESHOLD']}, |LFC|>={cfg['LFC_THRESHOLD']}):", nrow(res_sig), "\\n")

# ── Save results ──────────────────────────────────────────────────────────────
out_dir <- "{dge_dir}"

write.csv(res_df, file.path(out_dir, "dge_all_genes.csv"), row.names=FALSE)
write.csv(res_sig, file.path(out_dir, "dge_significant.csv"), row.names=FALSE)

wb <- createWorkbook()
addWorksheet(wb, "All genes")
writeData(wb, "All genes", res_df)
addWorksheet(wb, "Significant DEGs")
writeData(wb, "Significant DEGs", res_sig)
saveWorkbook(wb, file.path(out_dir, "dge_results.xlsx"), overwrite=TRUE)

# ── VST for visualizations ────────────────────────────────────────────────────
vst_data <- vst(dds, blind=FALSE)

# PCA plot
pca_data <- plotPCA(vst_data, intgroup="condition", returnData=TRUE)
pca_var  <- round(100 * attr(pca_data, "percentVar"))
p_pca <- ggplot(pca_data, aes(PC1, PC2, color=condition, label=name)) +
  geom_point(size=3, alpha=0.8) +
  geom_text_repel(size=2.5, max.overlaps=20) +
  xlab(paste0("PC1: ", pca_var[1], "%% variance")) +
  ylab(paste0("PC2: ", pca_var[2], "%% variance")) +
  theme_bw(base_size=12) +
  ggtitle("PCA — VST-normalized expression")
ggsave(file.path(out_dir, "pca_plot.pdf"), p_pca, width=7, height=5)
ggsave(file.path(out_dir, "pca_plot.png"), p_pca, width=7, height=5, dpi=200)

# Volcano plot
volcano_df <- res_df %>% mutate(
  sig = case_when(
    padj < {cfg['PADJ_THRESHOLD']} & log2FoldChange >= {cfg['LFC_THRESHOLD']}  ~ "Up",
    padj < {cfg['PADJ_THRESHOLD']} & log2FoldChange <= -{cfg['LFC_THRESHOLD']} ~ "Down",
    TRUE ~ "NS"
  )
)
top_genes <- volcano_df %>% filter(sig != "NS") %>% slice_min(padj, n=20)

p_volcano <- ggplot(volcano_df, aes(log2FoldChange, -log10(padj), color=sig)) +
  geom_point(size=0.8, alpha=0.6) +
  geom_text_repel(data=top_genes, aes(label=gene_id), size=2.5, max.overlaps=15) +
  scale_color_manual(values=c(Up="firebrick", Down="steelblue", NS="grey70")) +
  geom_vline(xintercept=c(-{cfg['LFC_THRESHOLD']}, {cfg['LFC_THRESHOLD']}), linetype="dashed", color="grey50") +
  geom_hline(yintercept=-log10({cfg['PADJ_THRESHOLD']}), linetype="dashed", color="grey50") +
  theme_bw(base_size=12) +
  ggtitle("Volcano plot") +
  labs(color="Direction")
ggsave(file.path(out_dir, "volcano_plot.pdf"), p_volcano, width=7, height=5)
ggsave(file.path(out_dir, "volcano_plot.png"), p_volcano, width=7, height=5, dpi=200)

# Heatmap — top 50 DEGs
if (nrow(res_sig) >= 2) {{
  top50 <- head(res_sig$gene_id, 50)
  mat_heat <- assay(vst_data)[top50, ]
  mat_heat <- t(scale(t(mat_heat)))
  ann_col  <- data.frame(condition=meta$condition, row.names=rownames(meta))
  pheatmap(mat_heat,
           annotation_col=ann_col,
           show_colnames=(ncol(mat_heat) <= 60),
           fontsize_row=7, fontsize_col=7,
           filename=file.path(out_dir, "heatmap_top50_degs.pdf"),
           width=10, height=10)
  cat("Heatmap saved.\\n")
}}

cat("\\n=== Analysis complete ===\\n")
cat("Outputs in:", out_dir, "\\n")
cat("  dge_all_genes.csv     — full results table\\n")
cat("  dge_significant.csv   — significant DEGs only\\n")
cat("  dge_results.xlsx      — Excel workbook\\n")
cat("  pca_plot.pdf/png      — PCA of cells\\n")
cat("  volcano_plot.pdf/png  — volcano plot\\n")
cat("  heatmap_top50_degs.pdf — expression heatmap\\n")
"""


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="SMARTseq HT scRNA-seq pipeline")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument("--steps", default="qc,trim,align,count,dge",
                        help="Comma-separated steps to run (default: all)")
    parser.add_argument("--resume", action="store_true",
                        help="Skip steps that already have .done flags")
    parser.add_argument("--cells", default=None,
                        help="Comma-separated cell IDs to process (default: all)")
    args = parser.parse_args()

    cfg = load_config(args.config)
    log = setup_logging(Path("logs"))
    steps = [s.strip() for s in args.steps.split(",")]

    log.info("=" * 60)
    log.info("SMARTseq HT Pipeline")
    log.info(f"Config: {args.config}")
    log.info(f"Steps:  {steps}")
    log.info(f"Resume: {args.resume}")
    log.info("=" * 60)

    # Discover cells
    fastq_dir = Path(cfg["FASTQ_DIR"])
    all_cells = discover_cells(fastq_dir, cfg["R1_SUFFIX"], cfg["R2_SUFFIX"], log)

    if args.cells:
        requested = set(args.cells.split(","))
        all_cells = [c for c in all_cells if c in requested]
        log.info(f"Subset requested — processing {len(all_cells)} cell(s).")

    if not all_cells:
        log.error("No cells to process. Check FASTQ_DIR and naming suffixes in config.")
        sys.exit(1)

    t0 = time.time()

    if "qc" in steps:
        step_fastqc(all_cells, cfg, log, args.resume)
    if "trim" in steps:
        step_trim(all_cells, cfg, log, args.resume)
    if "align" in steps:
        step_align(all_cells, cfg, log, args.resume)
    if "count" in steps:
        step_featurecounts(all_cells, cfg, log, args.resume)
    if "dge" in steps:
        step_dge(cfg, log, args.resume)

    elapsed = time.time() - t0
    log.info("=" * 60)
    log.info(f"Pipeline complete in {elapsed/60:.1f} min")
    log.info(f"Results in: {cfg['OUTPUT_DIR']}/")
    log.info("=" * 60)


if __name__ == "__main__":
    main()
