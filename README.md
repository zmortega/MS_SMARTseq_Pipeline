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
        ├── Combined + Clarke 2025 analysis (your cells + Clarke et al. mouse data)
        │       ├── Same outputs as combined analysis above
        │       ├── Pairwise cluster contrasts vs VAF / VRC (two-sided, up + down)
        │       ├── limma batch correction applied before embedding
        │       └── Clarke cells treated as additional strain "Clarke2025"
        └── NoMHCIIFilter analysis (all sorted cells, MHCII filter bypassed,
                    depth floor only)
                ├── Per-cell lineage typing (curated marker panels)
                ├── UMAP colored by cell type
                ├── Violins of genes of interest within focus cell types
                └── Per-cell expression workbook + all-genes CSVs
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
│       ├── cluster_pairwise_contrasts.xlsx      (merged folder only)
│       ├── heatmap_cluster_markers.pdf
│       ├── violin_plots_by_cluster.pdf
│       ├── cell_identity_merged_clusters.xlsx
│       ├── VAF_VRC_correlation_merged_clusters.xlsx
│       └── MCA_celltype_correlation_merged_clusters.xlsx
│   └── NoMHCIIFilter/
│       ├── umap_all454_by_cell_type.pdf               (cell count is derived)
│       ├── violin_GOI_by_cell_type.pdf
│       ├── GOI_expression_by_cell.xlsx
│       ├── expression_all_genes_by_cell_VST.csv
│       └── expression_all_genes_by_cell_counts.csv
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

### Combined gene set — an all-plate intersection

`common_genes <- Reduce(intersect, ...)` keeps only genes that survive **every**
plate's `rowSums >= 10` filter, so a gene genuinely absent from one plate leaves the
combined matrix for all plates. Current run: per-plate filters keep 11,206–17,375
genes, and the intersection across all six plates is **8,948 genes**. The merged
(your cells + Clarke) matrix intersects that again with Clarke's mappable genes:
**7,772 genes × 469 cells**.

This is not cosmetic. The combined UMAP picks HVGs from that intersection, so the
gene set change reshuffles Leiden clustering for *every* plate — **cluster numbering
is not comparable across runs where the plate roster changed.** Snapshot
`results/05_dge/` before re-running if you need old cluster IDs (see
`results/05_dge_pre_fix/` for the pre-alignment-fix state).

#### The violin-gene filter exemption

The canonical case is `Ptprc` (CD45): it has **exactly 0 counts across all 96
NODCD31 cells**, because that plate is a pure CD45− sort with no CD45+ reference
block. That is the sort working correctly, not a QC failure — but the intersection
rule turned it into a deletion, and one plate's honest zero removed the CD45-purity
QC plot for every plate. Col1a1/Col1a2 hit the same wall (a pure CD31 sort contains
no fibroblasts).

The per-plate filter is therefore now:

```r
keep <- rowSums(s_counts) >= 10 | rownames(s_counts) %in% violin_keep_ens
```

**A zero count is a measurement, not a missing value.** Exempting `VIOLIN_GENES`
keeps real VST values — at the floor, where counts are zero — flowing through to
`combined_expr` and `merged_expr`, so these genes are plottable *as zeros* rather
than absent. All 18 violin genes now resolve in **6/6 plates** and in the merged
matrix; the exemption force-retained 1–7 genes per plate that the filter would
otherwise have dropped, and grew the intersection by 14 genes (8,934 → 8,948).

**Cost:** exempted genes are forced into each plate's DESeq2 object even when
all-zero there. Their own DE statistics on such a plate are meaningless (expect
`NA` padj), and they add ~18 of ~12,000 genes to the multiple-testing burden.
**Do not extend this list to hundreds of genes.**

Both violin panels still fail soft — they previously crashed outright
(`pivot_longer(): cols must select at least one column`) when a gene resolved to
nothing:

| Panel | Behavior when a gene is still unavailable |
|---|---|
| `combined_plots/violin_plots_by_cluster.pdf` | Sources each gene from the **per-plate** VST matrices (`all_expr_list`), not from `combined_expr`. Plots from whichever plates retain the gene; the panel subtitle names the plates where it is absent. |
| `CombinedwithVAFPaperPlots/violin_plots_by_cluster.pdf` | **Skipped for genes missing from `merged_expr`, with the reason logged.** Cannot use the per-plate fallback: `merged_expr` is `limma::removeBatchEffect`-corrected across your cells + Clarke, so an uncorrected per-plate value would put two scales on one axis. |

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

**One-sided and upregulated-only** (`alternative="greater"`, candidates pre-filtered
to `log2FC >= MIN_LOG2FC`). It answers "what marks this cluster against all others
pooled" and **structurally cannot report depletion.** For the complementary question
see the next section.

### Pairwise cluster contrasts vs VAF / VRC (CombinedwithVAFPaperPlots only)

`cluster_pairwise_contrasts.xlsx` answers a different question from
`cluster_marker_genes.xlsx`: for an uncharacterized cluster, what is up **and down**
relative specifically to the VAF cluster and to the VRC cluster?

Depletion is the point. A cluster lacking both Col1a1 and Pecam1 is positive
evidence it is neither fibroblast nor endothelial, and one-vs-rest cannot show that.

| | `cluster_marker_genes.xlsx` | `cluster_pairwise_contrasts.xlsx` |
|---|---|---|
| Comparison | one-vs-rest (all other clusters pooled) | one-vs-one, against VAF and VRC separately |
| Test | Wilcoxon, one-sided (`greater`) | Wilcoxon, two-sided |
| Direction | upregulated only | both, via a `direction` column |
| Matrix | merged batch-corrected VST | merged batch-corrected VST |
| Thresholds | `log2FC ≥ MIN_LOG2FC`, BH `padj < MAX_PADJ` | `\|log2FC\| ≥ MIN_LOG2FC`, BH `padj < MAX_PADJ` |

`log2FC = mean VST(query) − mean VST(reference)`; VST differences *are* log2 fold
changes, so no separate model is fit.

**Reference clusters are re-derived every run**, never hardcoded — the modal merged
Leiden cluster of the Clarke VAF cells and of the Clarke VRC cells respectively.
Leiden IDs are not stable across runs whose plate roster or gene set changed, so a
literal `VAF_CLUSTER <- 3` would silently rot.

Workbook layout: one sheet per contrast named `C{query}_vs_{VAF|VRC}`, plus `Summary`
(cell counts and up/down tallies) and `Notes` (thresholds and the reference-cluster
assignment for that run).

The block skips itself, with a logged reason, if either reference cluster cannot be
located or if VAF and VRC cells land in the *same* cluster — that would mean the
merged clustering does not separate them, and every contrast would be meaningless.
Individual contrasts are skipped when either side has fewer than 3 cells.

> **Check the reference-cluster purity in the log before trusting these sheets.**
> It is printed as e.g. `VAF reference: cluster 1 (23/50 cells, 46.0%)`. The modal
> cluster wins even by a plurality, so a low percentage means the "VAF reference"
> is mostly *not* VAF cells and every `*_vs_VAF` sheet is baselined on a mixed
> population. In the current run VRC is clean (43/46, 93.5%) but **VAF is only
> 46.0%** — the Clarke VAF cells are split across clusters.

### Cluster composition bar charts

4 stacked bar charts in one 2×2 PDF: counts and proportions by strain×population
(12 colors) and by strain (4 colors). Segment labels shown for segments ≥5%.

### VAF/VRC correlation (Clarke et al. 2025, GSE292898)

Two distinct steps that are easy to conflate. **Marker genes are not used to
decide a cluster's identity.**

**1. Cluster → identity match: genome-wide, no marker list.** Each cluster's mean
VST profile is correlated (Pearson *and* Spearman) against the mean profiles of
VAF, VRC, and CD45pos, across **all shared genes — 7,772 in the merged
analysis** (`n_genes` column in the output workbook). Both up- and
down-regulation contribute; nothing is restricted to a marker panel. The
correlations discriminate rather than saturating — in the current merged run
cluster 3 is VAF 0.823 / VRC 0.592 / CD45pos 0.606, and cluster 4 is
VAF 0.610 / VRC 0.877 / CD45pos 0.615 (Pearson, n_genes = 7,772).

**2. Clarke reference labeling: unsupervised split, marker-oriented naming.**
Which *Clarke* cells count as VAF vs VRC is decided in two stages:

- the **split** is unsupervised — PCA on the top 500 variable genes, then k-means
  with k=2, using no marker information at all;
- the **naming** uses the marker panels (Col1a1/Col1a2/Timp3/Spp1/Thy1/Pdpn for
  VAF; Pecam1/Eng/Cdh5/Kdr/Tie1/Vwf for VRC) *only* to decide which of the two
  k-means clusters gets called "VAF". Orientation scores both panels
  (VAF markers minus VRC markers) so a cluster high in both cannot win by
  default.

**Orientation runs in Clarke's own symbol space** (`vaf_cnt_mat2`, before the
Ensembl reindexing and before any intersection with your genes). This is
deliberate: which Clarke cluster is VAF is a Clarke-internal question and must
not depend on which plates you have sequenced.

Three bugs previously lived here, all silent:

0. **Inverted orientation (the per-strain / `combined_plots` block).** The rule was
   `if (c1_vaf_score > c2_vaf_score)` — the VAF panel alone, with the VRC scores
   computed but never used. `Timp3` is strongly endothelial in this dataset (6.967
   in the endothelial cluster vs 0.645 in the fibroblast cluster), so it dominated
   the VAF panel mean and outvoted Col1a1/Col1a2. **The endothelial cluster
   (Pecam1 5.94, Kdr 7.59, Cdh5 5.00) was labeled VAF, and every VAF/VRC
   correlation in both workbooks came out inverted.** The old printout used
   `max()`/`min()` per panel, so it read as self-consistent no matter which way the
   call went and could never reveal the flip.
1. **Symbol vs Ensembl.** The panels are symbols; the matrix had been reindexed
   to Ensembl IDs. `intersect()` returned `character(0)`, which propagated
   `0-row matrix → colMeans → NaN → sc1 > sc2 = NA → ifelse(NA,1,2) = NA →
   cluster == NA → all-NA subscript`, yielding NA vectors for *both* groups.
   `%in%` never matches NA, so **every** CD45neg Clarke cell fell through to
   VRC. The tell in the log was `Clarke VAF cells: 96 | VRC cells: 96` for 96
   total cells.
2. **Intersection leakage.** Routing panels through `sym_to_ens_rev` (built from
   `combined_expr`, the all-plate intersection) meant a marker missing from any
   one plate became unusable. NODCD31 has **0 counts for Col1a1 and Col1a2** — a
   pure CD31 sort contains no fibroblasts — which deleted both canonical VAF
   markers and left a single usable gene per panel.

Guards now in place, in both blocks:

- orientation scores **both** panels and takes the difference (VAF minus VRC), so a
  cluster high in endothelial markers cannot win the VAF label no matter how one
  contaminating gene behaves;
- hard-error below `MIN_ORIENT_MARKERS = 3` usable markers per panel — one marker is
  a coin flip if that gene happens to be bimodal or dropout-prone;
- hard-error on any `NA` score, rather than letting it propagate to an all-NA subscript;
- the log prints **actual per-cluster values** and which cluster won, not `max()`/`min()`;
- an independent **canonical-marker sanity check** (Col1a1/Col1a2 vs Pecam1/Cdh5/Kdr):
  if the VAF cluster is not higher in collagen *and* lower in endothelial markers, it
  emits a warning and prints `*** WARNING: canonical marker check FAILED ***`;
- `stopifnot` that the VAF and VRC id sets partition the CD45neg cells with no NAs.

Current run, all clean:

```
cluster1 (n=100): VAF panel 0.506 | VRC panel 0.178 | diff +0.328  <- VAF
cluster2 (n=64):  VAF panel 1.169 | VRC panel 4.976 | diff -3.807  <- VRC
sanity: collagen VAF 0.739 vs VRC 0.523 | endothelial VAF 0.149 vs VRC 6.178
```

> **Any workbook generated before this fix has VAF and VRC swapped.** Re-run rather
> than reinterpreting old output.

### Is the clustering independent of these labels?

The **cluster assignment step is** unsupervised: HVG → PCA → UMAP → Leiden sees
only expression, and population labels are attached afterward purely for
coloring and tabulation. (Note the UMAP clusters are **Leiden**; k-means appears
only in the Clarke labeling above.)

The **matrix fed into clustering is not fully label-blind**, in three places:

| Step | Label dependence |
|---|---|
| Per-plate size factors | Estimated from that plate's reference wells only, i.e. a population-selected subset |
| `varianceStabilizingTransformation(blind=FALSE)` | Dispersions estimated under `~ condition` |
| `removeBatchEffect(design = ~ condition)` (merged only) | Explicitly *preserves* condition differences so they survive batch correction |

The third is the strongest: for Clarke cells `condition` is literally
`VAF` / `VRC` / `CD45pos_MHCIIpos`, so those labels enter the design matrix that
produces the corrected matrix the merged UMAP and Leiden clustering run on. This
is deliberate and standard — without it `removeBatchEffect` would strip
biological signal confounded with dataset — but it does mean the merged
embedding is **not** independent of the VAF/VRC assignment. A change to those
labels changes the merged clusters, not just the legend.

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

### NoMHCIIFilter analysis

Unsupervised clustering of every sorted cell with the MHCII expression filter
bypassed entirely. Answers "what structure is there across everything I
sorted?", independent of the MHCII gate that defines every other folder.

**"No MHCII filter" is literal, and is the only biological gate removed.** One
quality gate is still applied: a sequencing-depth floor of `NOFILT_MIN_GENES`
(default **500 genes detected**, mirroring `MIN_GENES_DETECTED` in
`config.yaml`). That is a data-quality criterion, not a biological one — no cell
is included or excluded here on the basis of what it expresses. 564 sorted →
**454 analyzed**, 110 removed.

The method is deliberately **identical** to the `combined_plots` UMAP so the two
are directly comparable — per-plate low-count filter (violin-gene exemption
included) → per-plate DESeq2 with size factors from that plate's own reference
wells → VST → intersect genes across plates → top `N_HVG` by variance → PCA →
UMAP for layout, plus a separate all-gene PCA → kNN → Leiden for the cluster
assignment. Same tuning constants, same seed. **The only difference is which
cells go in.**

| | `combined_plots` | `NoMHCIIFilter` |
|---|---|---|
| Gate | MHCII expression | sequencing depth only |
| Cells | 351 | 454 |
| Genes in intersection | 8,948 | 12,049 |
| Leiden clusters | 3 | 3 |

The gene intersection is *larger* here despite the same filter rule: more cells
per plate means more genes clear `rowSums >= 10`.

Cells removed by the depth floor, per plate: B6G7 15, B6MHCIIGFP 20, NOD 23,
NOD2 30, NODCD31 15, NODPDL1 7. The loss is concentrated in the CD45− MHCIIhi /
MHCIIlo subgates and in NODCD31's CD31+ block; a full plate × condition
breakdown is printed to the log every run.

Outputs:

| File | Contents |
|---|---|
| `violin_GOI_by_cell_type.pdf` | Violins of `NF_GENES_OF_INTEREST` within `NF_FOCUS_TYPES`, with n and detection count on each axis label |
| `GOI_expression_by_cell.xlsx` | Genes of interest quantified per cell (6 sheets, below) |
| `expression_all_genes_by_cell_VST.csv` | All genes × all passing cells, VST — for drilling into any other gene |
| `expression_all_genes_by_cell_counts.csv` | Same matrix, raw counts — for checking whether a gene is detected at all |
| `umap_all{N}_by_cell_type.pdf` | UMAP colored by per-cell type; cell count derived at runtime |

Workbook sheets: `GOI_by_cell_long` (one row per gene × cell: VST, raw count,
CPM, detected, cell type), `GOI_genes_x_cells_VST` (literal genes × cells matrix
with a `CELL_TYPE` header row), `Summary_by_cell_type`, `Cell_annotations`,
`Panel_definitions`, `Notes`. The CSV columns are `cell_id`, matching
`Cell_annotations` — take a cell ID from the workbook, look up any gene in the CSV.

#### Per-cell type assignment

**Every cell gets its own label. This does not come from clustering.**
Fibroblasts are ~5% of cells and never form their own Leiden cluster here, so a
cluster-level label cannot produce a fibroblast group at all.

Method, in order:

1. Raw counts for the passing cells — **all genes, deliberately not `nf_expr`**
2. `log1p(CPM)` per cell
3. Z-score each gene across the cells (zero-variance genes drop; ~27,600 remain)
4. Score each panel per cell = mean z-score of that panel's genes
5. Label = highest-scoring panel; `margin` = top minus runner-up

> **Why step 1 says "not `nf_expr`".** `nf_expr` is the cross-plate intersection,
> so one plate's honest zero deletes a gene for every plate. NODCD31 is a pure
> CD31 sort with no fibroblasts and therefore 0 counts for Col1a1/Col1a2 —
> scoring on `nf_expr` resolved the **fibroblast panel to 2 of 10 markers and
> pericyte to 0 of 7**, silently mistyping the exact population the analysis is
> about. Which lineage a cell belongs to is a property of that cell and must not
> depend on which other plates were sequenced. Same failure mode as the Clarke
> VAF/VRC orientation bug above. The script now hard-errors if a panel named in
> `NF_FOCUS_TYPES` resolves to fewer than 3 markers.

##### The marker panels

59 markers across 6 cell types. These are the literal lists in
`NF_LINEAGE_PANELS` at the top of `scripts/per_strain_plots.R` — edit them there.
"Found in last run" is how many resolved to a gene present in the data; the
script prints this every run and hard-errors if a panel named in
`NF_FOCUS_TYPES` falls below 3.

| Cell type | Markers | Found in last run | Genes |
|---|---|---|---|
| `Hematopoietic` | 10 | 10/10 | `Ptprc`, `Coro1a`, `Laptm5`, `Lcp1`, `Fcer1g`, `Ctss`, `Cd52`, `Cd74`, `Arhgdib`, `Cd48` |
| `Endothelial` | 10 | 10/10 | `Pecam1`, `Cdh5`, `Kdr`, `Tie1`, `Vwf`, `Eng`, `Esam`, `Plvap`, `Cldn5`, `Egfl7` |
| `Fibroblast` | 11 | 11/11 | `Col1a1`, `Col1a2`, `Col3a1`, `Col6a1`, `Dcn`, `Lum`, `Pdgfra`, `Postn`, `Fbln1`, `Mgp`, `Serpinf1` |
| `Pericyte_SMC` | 7 | 7/7 | `Acta2`, `Pdgfrb`, `Rgs5`, `Myh11`, `Des`, `Notch3`, `Cspg4` |
| `Endocrine_islet` | 12 | 12/12 | `Chga`, `Chgb`, `Scg2`, `Scg5`, `Ins1`, `Ins2`, `Gcg`, `Sst`, `Ppy`, `Pcsk1n`, `Resp18`, `Pcsk2` |
| `Acinar_ductal` | 9 | 8/9 | `Cela1`, `Ctrb1`, `Prss2`, `Cpa1`, `Amy2a5`, `Krt19`, `Krt18`, `Sox9`, `Spp1` |

`Acinar_ductal` sits at 8/9 because **`Cpa1` has zero counts across all 454
cells**, so z-scoring produces `NaN` and it is dropped. The gene is present in
the count matrix and the annotation — it is simply not expressed in this sort,
which is unsurprising for islet preparations with little exocrine carryover.

That is a harmless shortfall, but it is exactly the kind of thing to check
before trusting a call: a panel silently shrinking is how the `nf_expr` bug
above went unnoticed. Several acinar markers are near-absent here — `Amy2a5`
in 2 of 454 cells, `Prss2` in 3, `Ctrb1` in 10 — so the `Acinar_ductal` label
rests largely on `Krt18` (152 cells) and `Spp1` (105), both of which are broader
than the exocrine compartment. **Treat the 22 acinar/ductal calls with more
caution than the endothelial or hematopoietic ones.**

Panels are deliberately **non-overlapping** — verified: all 59 markers are
unique, no gene appears in two panels, because a shared marker would make the
argmax label meaningless. They are also
deliberately **broad** — these identify lineages, not subtypes.

**Marker provenance — read before publishing.** These panels are conventional
textbook lineage markers, assembled by hand. They are **not** taken from a
database or a citable source. Some overlap the repo's pre-existing fallback
`cell_db` (endothelial, pericyte) but the fibroblast, endocrine and
acinar/ductal panels are largely independent of it, so the repo now holds two
non-identical definitions of "fibroblast". For anything going into a figure,
replace them with a panel from a reference you trust, or cite deliberately.

> **More markers is not automatically better.** The score is a *mean*, so every
> marker carries weight 1/n and **a marker below the panel average drags the
> score down**. Measured on this dataset, adding to the fibroblast panel:
> `Col6a1` z=3.07 raises the mean 2.66→2.70 (kept); `Fap` z=1.38 lowers it to
> 2.55; `S100a4` z=0.29 lowers it to 2.45. `Fap` looks specific (3% elsewhere)
> but is detected in only 26% of fibroblasts — it marks *activated* fibroblasts,
> not fibroblasts. `Vim` is detected in 83% of fibroblasts and 53% of everything
> else, so it barely discriminates. If you want large panels to behave the way
> intuition expects, switch to rank-based scoring (UCell/AUCell-style), which is
> insensitive to panel size and dropout.

> **The fibroblast/pericyte boundary is genuinely fuzzy** — 56 cells have
> Fibroblast as runner-up, the closest pericytes at gaps of 0.025–0.031. Both
> are mesenchymal and pericytes express collagens. Endothelial and hematopoietic
> have 1 ambiguous cell each by comparison. Do not treat the fibroblast/pericyte
> split as sharp regardless of panel.

Cells whose top score beats the runner-up by less than `NF_AMBIGUOUS_MARGIN`
(0.25) are flagged `ambiguous`. **They are still counted in their top type** —
the flag exists so borderline cells can be excluded when it matters, not so they
are silently trusted. There is no "unassigned" class and no threshold on the
absolute score: every cell gets an argmax label regardless of how weak the
evidence is.

#### Legacy cluster-centric outputs

The earlier Leiden-cluster path (cluster UMAP, marker heatmap,
`cluster_lineage_calls.xlsx`, CellMarker and MCA workbooks per cluster) is kept,
not deleted. Set `NF_RUN_CLUSTER_OUTPUTS <- TRUE` to restore it. It is off by
default because its Wilcoxon loops and MCA download dominate this section's
runtime, and the per-cell typing supersedes it.

**Additive by construction.** The section runs last, reads only `counts_all` /
`meta_all` (snapshotted immediately before the MHCII filter overwrites `counts`
/ `meta`), writes only to `results/05_dge/NoMHCIIFilter/`, and assigns nothing
any earlier section reads. Verified: adding it left every prior figure identical
— 8,948 common genes, 351 combined cells, 3 combined clusters, Clarke VAF 50 /
VRC 46, 4 merged clusters, all unchanged.

#### Why the depth floor is not optional here

Removing the MHCII filter also removes the depth screen it was incidentally
performing. Running with no floor at all (`NOFILT_MIN_GENES <- 0`) was tried
first and produced two concrete failures:

1. **A pure artifact cluster.** All 564 cells gave **6** clusters, and cluster 6
   (n=29) sat alone at UMAP1 ≈ 17, far off the main manifold — shallow cells
   clustering by library size rather than by biology. The shallowest cell in the
   dataset detects **51 genes on 101 reads**. With the floor applied the island
   disappears and the structure resolves to 3 clean clusters.
2. **Collapsed size factors.** The per-plate estimator uses genes nonzero in
   *every* reference well, so a single near-empty reference well guts it. With
   no floor: B6G7 251 genes, NODPDL1 41, NODCD31 21, NOD 18, B6MHCIIGFP **4**,
   NOD2 **3**. Normalization at 3 genes is meaningless.

> **The floor improves this but does not fully fix it.** After filtering:
> NOD2 349, NODPDL1 260, B6G7 251 — but B6MHCIIGFP 122, NOD 107, NODCD31 57,
> all still tripping the `*** WARNING ***` that fires below 200. The residual
> cause is structural rather than depth: "nonzero in *every* reference well" is
> a strict criterion when a plate has only 6–15 reference wells left. Treat
> cross-plate comparisons in this folder with corresponding caution. The
> filtered analyses are unaffected — they retain more reference wells because
> the MHCII filter removed the shallow ones for different reasons.

### CombinedwithVAFPaperPlots analysis

Merges your MHCII-filtered VST data with Clarke et al. 2025 mouse cells
(GSE292898, 118/188 cells passing MHCII filter). Clarke cells are processed
through independent DESeq2 VST then **limma batch correction**
(`removeBatchEffect`) is applied to the merged matrix before embedding.
Clarke cells are labeled as `Clarke2025 CD45pos`, `Clarke2025 VAF`, or
`Clarke2025 VRC`. All combined_plots outputs are reproduced for this
merged dataset. The other 5 output folders are completely unaffected.

---

## Which analyses run

`scripts/per_strain_plots.R` has four switches near the top. Turning the first
two off takes a full run from **~15 min to 11 min 25 s** (measured, clean run,
all six plates).

That is a real but modest win, and worth being clear about why: most of the
remaining time is work the flags *cannot* skip — the six per-plate DESeq2/VST
fits are mandatory (see the dependency table below), and the Wilcoxon loops in
the two enabled folders are themselves expensive (NoMHCIIFilter cluster 1 alone
tests ~9,500 candidate genes, and the merged pairwise contrasts run six
two-sided comparisons over ~7,800 genes). Getting substantially below this
needs caching the per-plate VST matrices to disk, not more flags.

```r
RUN_PER_STRAIN_PLOTS <- FALSE   # results/05_dge/<STRAIN>_plots/
RUN_COMBINED_PLOTS   <- FALSE   # results/05_dge/combined_plots/
RUN_VAF_MERGED       <- TRUE    # results/05_dge/CombinedwithVAFPaperPlots/
RUN_NOMHCIIFILTER    <- TRUE    # results/05_dge/NoMHCIIFilter/
```

Nothing is deleted — flip a flag back to `TRUE` to restore that output exactly
as before.

> **`RUN_COMBINED_PLOTS = FALSE` skips the `combined_plots` folder, not the
> combined matrix it is named after.** Several things inside those sections are
> consumed by the analyses that are still on, so they run regardless of the
> flags:
>
> | Still runs | Why |
> |---|---|
> | The per-strain DESeq2/VST loop | Builds `all_expr_list` / `all_meta_list` → `combined_expr` + `combined_meta`, which CombinedwithVAFPaperPlots consumes. Only the loop's plots and spreadsheets are skipped, never its arithmetic. |
> | `combined_meta$strain_condition` | The merged metadata is built from it |
> | `strain_base_colors`, `sc_palette` | The merged UMAP and bar charts extend these palettes |
> | `MIN_LABEL_PROP` | Used by the merged bar charts |
> | `ref_profiles` (Clarke VAF/VRC/CD45pos means) | The merged VAF/VRC correlation reuses them |
> | CellMarker DB, `score_cluster`, `universe_size` | Merged and NoMHCIIFilter both score against them |
> | `sym_to_ens`, the openxlsx/ComplexHeatmap block, `HT_RASTER_DEVICE` | Used throughout |
>
> The gates are placed *around* these, so the shared work happens and only the
> output-writing is skipped. If you gate a new section, check for this class of
> dependency first — the failure mode is a mid-run `object '<x>' not found`
> after several minutes of successful work.

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
| `VIOLIN_GENES` | 18 symbols (see below) | Genes for cluster violin plots — **also exempt from every plate's low-count filter** |
| `MIN_ORIENT_MARKERS` | 3 | Minimum markers per panel before Clarke VAF/VRC orientation will run (hard-errors below this) |
| `NOFILT_MIN_GENES` | 500 | Depth floor for the **NoMHCIIFilter** analysis only — minimum genes detected per cell. Mirrors `MIN_GENES_DETECTED` in `config.yaml` (kept in sync by hand; this script does not parse the YAML). Set to 0 to disable. |
| `NF_GENES_OF_INTEREST` | c("Nod2", "Ciita") | Genes quantified per cell and plotted as violins in **NoMHCIIFilter** |
| `NF_FOCUS_TYPES` | c("Endothelial", "Fibroblast") | Cell types the violins focus on; must match `names(NF_LINEAGE_PANELS)`. A focus panel resolving to <3 markers is a hard error |
| `NF_AMBIGUOUS_MARGIN` | 0.25 | Below this top-vs-runner-up gap a cell is flagged `ambiguous` (still counted in its top type) |
| `NF_LINEAGE_PANELS` | 6 curated panels | Marker panels used for per-cell typing. See the provenance warning above before publishing |
| `NF_RUN_CLUSTER_OUTPUTS` | FALSE | Restore the legacy Leiden-cluster outputs in NoMHCIIFilter |
| `NF_VIOLIN_GENES` | c("Ciita", "Nod2") | Legacy cluster-violin path only |

`VIOLIN_GENES` is defined in the config block at the top of
`scripts/per_strain_plots.R`, not next to the plotting code — the exemption has to
be resolved before the strain loop runs. Symbols, not Ensembl IDs:

| Group | Genes |
|---|---|
| CD45 purity | `Ptprc` |
| VAF panel | `Col1a1`, `Col1a2`, `Timp3`, `Spp1`, `Thy1`, `Pdpn` |
| VAF extras | `S100a4`, `Fn1` |
| VRC panel | `Pecam1`, `Eng`, `Cdh5`, `Kdr`, `Tie1`, `Vwf`, `Esam` |
| Innate sensing / MHCII TF | `Nod2`, `Ciita` |

> **`Timp3` is not a usable VAF marker in this dataset.** It runs ~10× higher in the
> endothelial (VRC) cluster than the fibroblast (VAF) cluster — 6.967 vs 0.645 —
> despite sitting in the VAF panel. Orienting on the VAF panel alone let Timp3
> outvote Col1a1/Col1a2 and invert the entire VAF/VRC call (see below). It is kept
> in the list to be *plotted*, not to be read as fibroblast evidence.

---

## Required R packages

Auto-installed by the script if missing:
`igraph`, `uwot`, `BiocManager`, `openxlsx`, `ComplexHeatmap`, `circlize`,
`RColorBrewer`, `limma`, `remotes`

Manual install required before first run — these are loaded with a bare
`library()` and will hard-fail the run if absent:
```r
install.packages(c("ggplot2", "ggrepel", "dplyr", "tidyr", "patchwork",
                   "scales", "ragg"))
BiocManager::install("DESeq2")
```

**`ragg` and heatmap rendering.** The heatmaps use `use_raster=TRUE`, which needs
a working bitmap device to write a temp PNG. ComplexHeatmap defaults to
`grDevices::png()`, which on macOS routes through cairo/X11 — on a machine
without XQuartz that fails at `draw()` time with an opaque
`unable to open .heatmap_body_*.png` and kills the run *after* most outputs are
already written. `ragg::agg_png` is self-contained, so the script prefers it via
`HT_RASTER_DEVICE` and falls back to the stock device when `ragg` is absent.

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

- **Plates/strains:** NOD, NOD2, B6G7, B6MHCIIGFP, NODPDL1, NODCD31 (`scripts/per_strain_plots.R`
  discovers these dynamically from `data/metadata.csv` — no code changes needed to add NOD3,
  NOD4, etc.)
- **Strain grouping:** `data/metadata.csv` has both a `strain` column (literal plate, used for
  independent per-plate normalization/DESeq2/output folders) and a `strain_group` column
  (biological grouping used for combined-analysis plots/colors). NOD-family plates (NOD, NOD2,
  NOD3, ...) share `strain_group = NOD` and are shown together by default in combined plots;
  NODPDL1 and NODCD31 are biologically distinct and each keep their own group.
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
- **NODCD31:** entirely CD45− MHCII+, with **no CD45pos_MHCIIpos wells at all** — it reuses
  NODPDL1's 18/78 well split but sorts on CD31 instead: `CD45neg_MHCIIpos_CD31neg` (A1–B6,
  used as that plate's normalization reference) and `CD45neg_MHCIIpos_CD31pos` (B7–H12).
  Consequences: `Ptprc` has exactly 0 counts across all 96 cells, as do `Col1a1`/`Col1a2`
  (a pure CD31 sort contains no fibroblasts). All three are real measurements, not QC
  failures — see the violin-gene filter exemption above.
- **Batches:** B6G7 and B6MHCIIGFP = batch1 (L001); NOD and NOD2 = batch2 (L002); NODPDL1 =
  batch3 (L001+L002 merged); NODCD31 = batch4
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
