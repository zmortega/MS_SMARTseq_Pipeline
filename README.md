# SMARTseq HT scRNA-seq Pipeline

Full-length, plate-based single-cell RNA-seq pipeline for differential expression,
dimensionality reduction, cluster analysis, and cell type identification.
Optimized for **Takara SMARTseq mRNA HT** kits on **Apple Silicon Mac** with mouse (**GRCm39 / GENCODE vM33**).

---

## Pipeline overview

```
FASTQ (R1/R2 per cell)
    │
    ▼
[1] FastQC          → per-cell read quality reports
    │
    ▼
[2] Trimmomatic     → adapter trimming, quality filtering
    │
    ▼
[3] HISAT2          → splice-aware alignment to GRCm39
    │
    ▼
[4] featureCounts   → gene-level count matrix (all cells → one matrix)
    │
    ▼
[5] per_strain_plots.R
        ├── Metadata/count-matrix alignment (meta reordered to colnames(counts))
        ├── MHCII expression filter (H2-Aa & H2-Ab1, applied globally)
        ├── Per-strain DESeq2 (independent per plate, normalized to that plate's own reference wells)
        │       ├── Volcano + violin plots (1 pair per contrast — contrasts built dynamically
        │       │     from whichever conditions are present, e.g. 3 for MHCIIhi/lo plates, 1 for NODPDL1)
        │       ├── Expression matrix CSV (mean VST + detection stats + upregulation rankings)
        │       ├── Per-strain UMAP (by cluster + by population)
        │       └── Cell identity Excel (CellMarker 2.0 gene set scoring)
        ├── Combined analysis (all MHCII-filtered cells, grouped by strain_group)
        │       ├── Combined UMAP (by Leiden cluster + by strain/population)
        │       ├── Per-strain panels on comprehensive UMAP
        │       ├── Cluster composition bar charts (4 charts, 2x2 PDF)
        │       ├── Cluster marker gene Excel (Wilcoxon rank-sum)
        │       ├── Cluster marker gene heatmap
        │       ├── Cluster violin plots (selected genes)
        │       ├── Cell identity Excel (CellMarker 2.0 gene set scoring)
        │       ├── VAF/VRC correlation Excel (Clarke et al. 2025, GSE292898)
        │       └── MCA cell type correlation Excel (Mouse Cell Atlas)
        └── Combined + Clarke 2025 analysis (your cells + Clarke et al. mouse data)
                ├── Same outputs as combined analysis above
                ├── limma batch correction applied before embedding
                └── Clarke cells treated as additional strain "Clarke2025"
```

---

## Quick start

### 1. Install dependencies

```bash
conda env create -f environment.yml
conda activate smartseq_ht
```

### 2. Build the HISAT2 index (one-time)

```bash
hisat2-build -p 8 reference/GRCm39.primary_assembly.genome.fa reference/hisat2_index/genome
```

### 3. Prepare your metadata

Edit `data/metadata.csv` — one row per cell with columns:
`cell_id`, `strain`, `well`, `condition`, `batch`, `strain_group`

`cell_id` must exactly match the FASTQ stem (everything before `_R1_001.fastq.gz`).

See `data/metadata_example.csv` for format reference.

**Plate layouts currently in `data/metadata.csv`** (564 cells):

| Plate (`strain`) | `strain_group` | `batch` | Layout |
|---|---|---|---|
| B6G7 | B6G7 | batch1 | A1–A12 CD45+ MHCII+ · rest split MHCIIhi/MHCIIlo (42/42) |
| B6MHCIIGFP | B6MHCIIGFP | batch1 | A1–A12 CD45+ MHCII+ · 30 MHCIIhi / 42 MHCIIlo (84 cells) |
| NOD | NOD | batch2 | A1–A12 CD45+ MHCII+ · 42 MHCIIhi / 42 MHCIIlo |
| NOD2 | NOD | batch2 | A1–A12 CD45+ MHCII+ · 42 MHCIIhi / 42 MHCIIlo |
| NODPDL1 | NODPDL1 | batch3 | A1–B6 CD45+ MHCII+ (18) · B7–H12 CD45− MHCII+ (78) |
| NODCD31 | NODCD31 | batch4 | A1–B6 CD45− MHCII+ **CD31−** (18) · B7–H12 CD45− MHCII+ **CD31+** (78) |

NODCD31 has **no CD45+ MHCII+ wells** — it reuses NODPDL1's 18/78 well split but
puts CD31− where the reference block normally sits. See *Per-strain
normalization* below for how that plate is baselined.

### 4. Configure

Edit `config.yaml`:
- Set `FASTQ_DIR` to your folder of `.fastq.gz` files
- Set `METADATA_FILE` to your metadata CSV

### 5. Run alignment and counting

```bash
python pipeline.py --config config.yaml
```

> **`--resume` is safe for `qc`/`trim`/`align`, but NOT for `count`.**
> Those three steps write a `.done` flag per cell, so `--resume` correctly
> reprocesses only new cells. The `count` step writes a *single* global flag at
> `results/04_counts/all_cells/.done` that records only *that* counting
> happened, never *which* cells were counted. Adding a plate and re-running with
> `--resume` therefore skips featureCounts entirely and leaves you with a stale
> `counts_clean.txt` — the new cells align fine, then silently never reach the
> count matrix or any downstream plot.
>
> When adding a plate, clear the flag first:
>
> ```bash
> rm -rf results/04_counts/all_cells
> python pipeline.py --config config.yaml --steps count
> ```
>
> Then verify the column count before going further — should be
> `1 + 5 annotation + n_cells`:
>
> ```bash
> head -2 results/04_counts/counts_clean.txt | tail -1 | awk '{print NF}'
> ```

> **Do not pass `--cells` to the `count` step.** `step_featurecounts()` builds
> its BAM list from whatever cell list it is handed and then overwrites
> `counts_clean.txt`. Restricting to a new plate rebuilds the matrix containing
> *only* that plate and discards every other cell.

### 6. Run all DGE, UMAP, and downstream analyses

```bash
cd ~/Desktop/MS_SMARTseq_Pipeline
Rscript scripts/per_strain_plots.R 2>&1 | tee per_strain_plots.log
```

> **Warning:** This script wipes and rebuilds `results/05_dge/` entirely on every run.
> Upstream folders (01-04, qc_summary) are never touched.

---

## External reference files required

Place these in `reference/` before running:

| File | Source | Used for |
|---|---|---|
| `gencode.vM33.primary_assembly.annotation.gtf` | GENCODE vM33 | Gene symbol mapping |
| `GRCm39.primary_assembly.genome.fa` | GENCODE vM33 | HISAT2 alignment |
| `GSE292898_teyton_don_2025_processed_mouse_raw_counts_matrix.csv.gz` | GEO GSE292898 | Clarke et al. merged analysis |
| `GSE292898_README_sample_annotations.txt` | GEO GSE292898 | Clarke cell metadata |

The Mouse Cell Atlas reference (`ref_MCA.rda`) is downloaded automatically at runtime and cached to `/tmp/`.

---

## File naming convention

```
{STRAIN}_{WELL}_{SAMPLE}_{LANE}_R1_001.fastq.gz
{STRAIN}_{WELL}_{SAMPLE}_{LANE}_R2_001.fastq.gz
```

Examples: `NOD_A1_S1_L002_R1_001.fastq.gz`, `B6G7_H12_S96_L001_R1_001.fastq.gz`

---

## Output structure

```
results/
├── 01_fastqc/
├── 01_fastqc_multiqc/
├── 02_trimmed/
├── 03_aligned/
├── 03_aligned_multiqc/
├── 04_counts/
│   ├── raw_counts.txt
│   └── counts_clean.txt
├── 05_dge/
│   ├── NOD_plots/                            (file count varies with # of contrasts)
│   │   ├── expression_matrix_NOD.csv
│   │   ├── volcano_NOD_*.pdf                 (1 per contrast)
│   │   ├── violins_NOD_*.pdf                 (1 per contrast)
│   │   ├── umap_NOD_by_cluster.pdf
│   │   ├── umap_NOD_by_population.pdf
│   │   ├── cell_identity_NOD_clusters.xlsx
│   │   └── MCA_celltype_correlation_NOD_clusters.xlsx
│   ├── NOD2_plots/                           (same layout, 3 contrasts)
│   ├── B6G7_plots/                           (same layout, 3 contrasts)
│   ├── B6MHCIIGFP_plots/                     (same layout, 3 contrasts)
│   ├── NODPDL1_plots/                        (same layout, 1 contrast — single CD45neg_MHCIIpos condition)
│   ├── NODCD31_plots/                        (same layout, 1 contrast — CD31+ vs CD31−, baselined on CD31−)
│   ├── combined_plots/
│   │   ├── umap_all_by_cluster.pdf
│   │   ├── umap_all_by_strain_population.pdf
│   │   ├── umap_per_strain_on_comprehensive.pdf
│   │   ├── barplots_cluster_composition.pdf
│   │   ├── cluster_marker_genes.xlsx
│   │   ├── heatmap_cluster_markers.pdf
│   │   ├── violin_plots_by_cluster.pdf
│   │   ├── cell_identity_combined_clusters.xlsx
│   │   ├── VAF_VRC_correlation_combined_clusters.xlsx
│   │   └── MCA_celltype_correlation_combined_clusters.xlsx
│   └── CombinedwithVAFPaperPlots/
│       ├── umap_merged_by_cluster.pdf
│       ├── umap_merged_by_strain_population.pdf
│       ├── umap_per_strain_on_merged.pdf
│       ├── barplots_cluster_composition.pdf
│       ├── cluster_marker_genes.xlsx
│       ├── heatmap_cluster_markers.pdf
│       ├── violin_plots_by_cluster.pdf
│       ├── cell_identity_merged_clusters.xlsx
│       ├── VAF_VRC_correlation_merged_clusters.xlsx
│       └── MCA_celltype_correlation_merged_clusters.xlsx
└── qc_summary/
    ├── cell_qc_metrics.csv
    ├── cell_qc_barplots.pdf
    └── metadata_qc_pass.csv
```

---

## per_strain_plots.R — detailed documentation

### Metadata / count matrix alignment

`featureCounts` orders its columns by the BAM order it was handed
(`B6G7_A10, B6G7_A11, B6G7_A12, B6G7_A1, ...`), which is **not** the row order of
`data/metadata.csv` (`B6G7_A1, B6G7_A10, ...`). Any logical mask computed from the
count columns — the MHCII filter, in particular — was previously applied
positionally to `meta`, which silently paired the wrong metadata row with each
cell and dropped cells that had actually passed.

Immediately after load, `meta` is now reordered to `colnames(counts)` and asserted:

```r
meta <- meta[colnames(counts), , drop=FALSE]
stopifnot(identical(rownames(meta), colnames(counts)),
          identical(meta$cell_id,   colnames(counts)))
```

The MHCII mask is named by `cell_id`, and both `counts` and `meta` are subset by
**cell_id character vector**, never by the logical mask, so they cannot
desynchronize downstream. A cell present in `counts` with no metadata row is a
hard error; metadata rows with no count column are reported and dropped.

### MHCII expression filter

Applied globally before any analysis. A cell is retained only if both **H2-Aa**
and **H2-Ab1** exceed a minimum log1p(CPM) threshold (default: log1p(5)).
Cells failing this filter are excluded from all downstream analyses.

```r
MHCII_VST_MIN <- 1.0
MHCII_GENES   <- c("H2-Aa", "H2-Ab1")
```

### Per-strain normalization

Each plate is treated as a fully independent experiment. Size factors are
estimated using only that plate's own **reference wells** — normally the
**CD45pos_MHCIIpos (A1-A12) wells**. No cross-plate pooling occurs.

A plate is not required to sort CD45+ wells. `REF_CONDITION_OVERRIDES` in
`scripts/per_strain_plots.R` sets a per-plate baseline; **NODCD31** is entirely
CD45− MHCII+ (no CD45+ wells at all) and is normalized on its own **CD31−
(A1–B6)** wells instead.

Because of this, "is the reference condition" and "is a CD45+ cell" are no
longer the same test. Normalization and DESeq2 baselining go through
`ref_condition_for(strain)`; the biological CD45 gate (which cells appear in the
CD45− UMAPs and the merged cluster-distribution tab) goes through
`is_cd45neg()`, which reads the condition name. Using the former for the latter
silently drops NODCD31's 18 CD31− cells.

### Per-strain DESeq2

Design: `~ condition`. VST via `varianceStabilizingTransformation()`.

### Per-strain contrasts

Contrasts are no longer hardcoded — `build_contrasts()` generates them per plate from
whichever `condition` values are actually present in `data/metadata.csv`: that plate's
reference condition vs. each non-reference condition, plus all pairwise comparisons
among the non-reference conditions. This reproduces the original fixed 3-contrast
design for NOD/NOD2/B6G7/B6MHCIIGFP and collapses to a single contrast for plates with
one non-reference population:

| Plate(s) | Reference | Contrasts |
|---|---|---|
| NOD, NOD2, B6G7, B6MHCIIGFP | CD45pos_MHCIIpos | CD45+ MHCIIpos vs CD45− MHCIIhi · CD45+ MHCIIpos vs CD45− MHCIIlo · CD45− MHCIIhi vs CD45− MHCIIlo |
| NODPDL1 | CD45pos_MHCIIpos | CD45+ MHCIIpos vs CD45− MHCIIpos |
| NODCD31 | CD45neg_MHCIIpos_CD31neg | CD45− MHCIIpos CD31− vs CD45− MHCIIpos CD31+ |

`build_contrasts()` now hard-errors if a plate's configured reference condition is
absent from its metadata, rather than silently producing zero contrasts.

Each contrast gets one volcano PDF and one violin PDF (top 10 up + top 10 down DEGs,
gene symbols, jittered VST points).

### Combined gene set — adding a plate shrinks it

`common_genes <- Reduce(intersect, ...)` keeps only genes that survive **every**
plate's `rowSums >= 10` filter, so a gene genuinely absent from one plate leaves the
combined matrix for all plates. Adding NODCD31 took the intersection from
**12,755 → 12,055 genes (−700, 5.5%)**.

This is not cosmetic. The combined UMAP picks HVGs from that intersection, so the
gene set change reshuffles Leiden clustering for *every* plate — **cluster numbering
is not comparable across runs where the plate roster changed.** Snapshot
`results/05_dge/` before re-running if you need old cluster IDs (see
`results/05_dge_pre_fix/` for the pre-alignment-fix state).

The most visible casualty is `Ptprc` (CD45): it has **exactly 0 counts across all 96
NODCD31 cells**, because that plate is a pure CD45− sort with no CD45+ reference
block. That is the sort working correctly, not a QC failure.

Consequences for the two violin panels, which both plot the hardcoded `VIOLIN_GENES`:

| Panel | Behavior |
|---|---|
| `combined_plots/violin_plots_by_cluster.pdf` | Sources each gene from the **per-plate** VST matrices (`all_expr_list`), not from `combined_expr`. Ptprc plots from the 5 plates that retain it, and the panel subtitle names the plates where it is absent. |
| `CombinedwithVAFPaperPlots/violin_plots_by_cluster.pdf` | **Skipped for genes missing from `merged_expr`, with the reason logged.** Cannot use the per-plate fallback: `merged_expr` is `limma::removeBatchEffect`-corrected across your cells + Clarke, so an uncorrected per-plate value would put two scales on one axis. |

Both panels previously crashed outright (`pivot_longer(): cols must select at least
one column`) when a gene resolved to nothing. Both now fail soft.

> If the merged Ptprc violin matters for a figure, the only way to recover it is to
> exclude NODCD31 from the combined analysis, which trades away its presence in the
> combined and merged UMAPs.

Combined UMAP filenames embed the actual cell count (`umap_all351_by_cluster.pdf`)
and are derived at runtime — the previously hardcoded `umap_all372_*` names silently
became wrong as soon as the cell count changed.

### Expression matrix columns

| Column | Description |
|---|---|
| `ensembl_id` | Ensembl gene ID |
| `gene_symbol` | Common gene name (GENCODE vM33) |
| `mean_VST_{condition}` | Mean VST expression, one column per condition present on that plate (e.g. `mean_VST_CD45pos_MHCIIpos`, `mean_VST_CD45neg_MHCIIhi`, ...) |
| `n_cells_{condition}` | Number of cells in that condition on that plate |
| `pct_detected_{condition}` | % of those cells with ≥1 raw read for the gene |
| `median_CPM_detected_{condition}` | Median CPM computed **only** across the cells that detected the gene (`NA` if none) |
| `rank_upregulated_{condition}` | Upregulation rank vs the reference condition (1 = most upregulated), one column per non-reference condition (e.g. `rank_upregulated_MHCIIhi`/`MHCIIlo` for most plates, `rank_upregulated_MHCIIpos` for NODPDL1) |

Ranks use combined score `log2FC × -log10(padj)`.

**Reading mean_VST vs. the detection columns.** `mean_VST` is a mean of a
compressive log-like scale, so it cannot distinguish "off in every cell" from
"high in a minority of cells" — a gene expressed strongly in 4 of 30 cells lands
near the VST floor and reads as absent. A **low `mean_VST` with a high
`median_CPM_detected` and low `pct_detected`** is the signature of a bimodal /
subset-expressed gene (e.g. `Cd274` in NODPDL1), not of true absence. Always
check the detection columns before calling a gene negative.

### Per-strain UMAPs

- `umap_{STRAIN}_by_cluster.pdf` — Leiden clusters (kNN graph built on UMAP coordinates for small-n stability, HVG UMAP for visualization)
- `umap_{STRAIN}_by_population.pdf` — colored by CD45+/MHCIIhi/MHCIIlo condition

### Cell identity scoring (per-strain and combined)

Fisher's exact test + Jaccard similarity against CellMarker 2.0 mouse reference
(downloaded at runtime, cached). Falls back to 20-cell-type built-in set if
download fails. When Wilcoxon markers are insufficient, falls back to top
expressed genes.

Excel workbook format (all folders):
- **Summary tab** — top 5 cell type candidates per cluster with rank, cell_type, fisher_pval
- **Per-cluster tabs** — full ranked list of all tested cell types with Fisher p-value, BH-adjusted p-value, Jaccard similarity, odds ratio, and overlapping gene symbols

### Combined UMAP (all strain groups)

**UMAP embedding:** top 2000 HVG PCA → 80% variance PCs → `uwot::umap()`

**Leiden clustering:** all-gene PCA → 80% variance PCs → kNN graph → Leiden

**Per-strain-group panels:** one panel per `strain_group` (not per literal plate — NOD-family
plates share a panel) showing that group's cells on shared comprehensive UMAP coordinates,
colored by comprehensive cluster assignment. Colors and legend sizing scale dynamically with
however many strain groups/conditions are present, rather than a fixed 4×3 palette.

### Cluster marker genes (Wilcoxon rank-sum)

One-vs-rest Wilcoxon on VST matrix. Pre-filter: log2FC ≥ 0.5. BH correction.
Top 100 per cluster in Excel, top 15 per cluster in heatmap.

### Cluster composition bar charts

4 stacked bar charts in one 2×2 PDF: counts and proportions by strain×population
(12 colors) and by strain (4 colors). Segment labels shown for segments ≥5%.

### VAF/VRC correlation (Clarke et al. 2025, GSE292898)

Pearson and Spearman correlation of each cluster's mean VST profile against
mean expression profiles of VAF, VRC, and CD45pos populations from Clarke et al.
mouse pancreatic islet data. CD45neg cells are assigned to VAF or VRC via k-means
clustering scored against known marker genes (Col1a1/Col1a2 for VAF;
Pecam1/Eng/Cdh5 for VRC).

### CD45− MHCII+ cluster distribution sheet

`VAF_VRC_correlation_merged_clusters.xlsx` (CombinedwithVAFPaperPlots) carries an
extra `CD45neg_Cluster_Distribution` sheet answering "where does each CD45−
MHCII+ cell land, cluster-wise, per strain?" The CD45+ MHCII+ wells are
normalization/reference only and are **excluded**; Clarke VAF/VRC cells are kept
as labeled reference rows (`Clarke2025 (ref)`), Clarke CD45pos is excluded.

| Table | Contents |
|---|---|
| 1 | Cell counts by strain, all CD45− subgates pooled |
| 2 | Row percentages by strain (each row sums to 100%) |
| 3 | Cell counts by strain × MHCII subgate (MHCIIhi / MHCIIlo / MHCIIpos / VAF / VRC) |
| 4 | Row percentages by strain × subgate |

Each row's dominant cluster is highlighted green; a bold `All CD45- cells`
column-total row closes every table. Clusters are the merged all-gene PCA +
Leiden clusters.

### Mouse Cell Atlas correlation

Pearson correlation against `ref_MCA` (713 mouse cell types, 8601 genes) from
the clustifyrdata package. Downloaded as a single `.rda` file at runtime,
cached to `/tmp/ref_MCA.rda`. Top 5 matches per cluster in Summary sheet;
full ranked list (713 cell types) in per-cluster sheets.

### CombinedwithVAFPaperPlots analysis

Merges your MHCII-filtered VST data with Clarke et al. 2025 mouse cells
(GSE292898, 117/188 cells passing MHCII filter). Clarke cells are processed
through independent DESeq2 VST then **limma batch correction**
(`removeBatchEffect`) is applied to the merged matrix before embedding.
Clarke cells are labeled as `Clarke2025 CD45pos`, `Clarke2025 VAF`, or
`Clarke2025 VRC`. All combined_plots outputs are reproduced for this
merged dataset. The other 5 output folders are completely unaffected.

---

## Tuning reference

| Parameter | Default | Effect |
|---|---|---|
| `MHCII_VST_MIN` | 1.0 | Minimum VST for MHCII filter |
| `N_HVG` | 2000 | HVGs for UMAP PCA |
| `VAR_THRESHOLD` | 0.80 | Cumulative variance for PC selection |
| `UMAP_N_NEIGHBORS` | 15 | UMAP neighborhood size |
| `UMAP_MIN_DIST` | 0.3 | UMAP point spread |
| `LEIDEN_RESOLUTION` | 0.5 | Combined Leiden granularity |
| `UMAP_SEED` | 42 | Reproducibility seed |
| `MIN_LOG2FC` | 0.5 | Wilcoxon pre-filter threshold |
| `MAX_PADJ` | 0.05 | Adjusted p-value cutoff |
| `TOP_EXCEL` | 100 | Marker genes per cluster in Excel |
| `TOP_HEATMAP` | 15 | Marker genes per cluster in heatmap |
| `MIN_LABEL_PROP` | 0.05 | Minimum bar segment size for label |
| `VIOLIN_GENES` | c("Ptprc") | Genes for cluster violin plots |

---

## Required R packages

Auto-installed by the script if missing:
`igraph`, `uwot`, `BiocManager`, `openxlsx`, `ComplexHeatmap`, `circlize`,
`RColorBrewer`, `limma`, `remotes`

Manual install recommended before first run:
```r
install.packages(c("ggplot2", "ggrepel", "dplyr", "tidyr", "patchwork", "scales"))
BiocManager::install("DESeq2")
```

---

## To run on a different machine

1. Update `base_dir` at the top of `scripts/per_strain_plots.R`
2. Place reference files in `reference/` (see table above)
3. Activate conda environment and run:

```bash
Rscript scripts/per_strain_plots.R 2>&1 | tee per_strain_plots.log
```

---

## Advanced usage

```bash
python pipeline.py --config config.yaml --resume        # resume failed run
python pipeline.py --config config.yaml --steps align,count,dge
python qc_summary.py --config config.yaml               # QC flagging
```

---

## Experiment-specific notes

- **Plates/strains:** NOD, NOD2, B6G7, B6MHCIIGFP, NODPDL1 (`scripts/per_strain_plots.R` discovers
  these dynamically from `data/metadata.csv` — no code changes needed to add NOD3, NOD4, etc.)
- **Strain grouping:** `data/metadata.csv` has both a `strain` column (literal plate, used for
  independent per-plate normalization/DESeq2/output folders) and a `strain_group` column
  (biological grouping used for combined-analysis plots/colors). NOD-family plates (NOD, NOD2,
  NOD3, ...) share `strain_group = NOD` and are shown together by default in combined plots;
  NODPDL1 is biologically distinct and keeps its own group.
- **NOD2:** originally delivered mislabeled as "NODPDL1" — same biological context as NOD (a
  second, independent plate), renamed to NOD2 and relabeled `strain_group = NOD`.
- **NODPDL1 (current):** a distinct plate/biology from NOD. Its 96 wells were originally split
  across two sequencing lanes (L001 + L002) and were concatenated per well into single R1/R2
  FASTQ pairs before running through the pipeline. Condition scheme differs from the other
  plates — no MHCIIhi/MHCIIlo split, just a single CD45- MHCII+ population:
  `CD45pos_MHCIIpos` (A1–A12, B1–B6, used as the normalization reference) and
  `CD45neg_MHCIIpos` (B7–B12, C1–H12).
- **NOD / B6G7 / B6MHCIIGFP / NOD2 conditions:** `CD45pos_MHCIIpos` (A1–A12), `CD45neg_MHCIIhi`
  (B1–E6), `CD45neg_MHCIIlo` (E7–H12)
- **Plate quirks:** B6MHCIIGFP wells H1–H12 empty; NOD2 (formerly mislabeled NODPDL1) layout
  physically flipped but metadata labels are biologically correct
- **Batches:** B6G7 and B6MHCIIGFP = batch1 (L001); NOD and NOD2 = batch2 (L002); NODPDL1 =
  batch3 (L001+L002 merged)
- **Normalization reference:** `CD45pos_MHCIIpos` wells, per plate (well range varies by plate —
  see metadata)
- **Contrasts:** built dynamically per plate from whichever conditions are present (3 contrasts
  for the CD45pos/MHCIIhi/MHCIIlo design, 1 contrast for NODPDL1's simpler 2-condition design)
- **DESeq2 design:** `~ condition` (no batch term; each strain is single-batch)
- **VST method:** `varianceStabilizingTransformation()` (bypasses sample-size check in `vst()`)
- **Aligner:** HISAT2 (switched from STAR due to Apple Silicon RAM constraints)
- **Annotation:** GENCODE vM33 primary assembly GTF
- **featureCounts:** run with `-p` (declares paired-end BAMs, required by subread ≥2.0.3 or it
  hard-errors with "Paired-end reads were detected in single-end read library"). `-p` alone
  (no `--countReadPairs`) preserves the original per-read counting behavior — each mate is
  still counted individually, not as one fragment.
- **QC summary:** `qc_summary.py` now parses HISAT2's `*_hisat2.log` stderr summary instead of
  STAR's `Log.final.out`. HISAT2 doesn't report a "uniquely mapped" figure, so `uniquely_mapped`
  is approximated as `total_reads × overall_alignment_rate` — used only to gate on read depth.
- **Clarke et al. 2025:** Cell Reports 44, 116189. GEO: GSE292898. VAFs are CD45neg MHC-II+ fibroblastic cells from NOD mouse pancreatic islets.
