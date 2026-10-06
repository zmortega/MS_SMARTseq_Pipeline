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
        ├── MHCII transcript gate (H2-Aa OR H2-Ab1 present; sorted MHCII− cells bypass)
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
        ├── ClarkeOnly analysis (Clarke/Don GSE292898 ALONE, no cells of yours)
        │       │     RUNS FIRST — the Combined analysis consumes its output
        │       ├── UMAP by Leiden cluster + UMAP by published sort gate
        │       ├── Cluster marker heatmap + marker workbook
        │       ├── Violins of Ptprc / Col1a1 / Col1a2 / Pecam1 per cluster
        │       └── VAF / VRC / CD45pos reference profiles  ──┐
        │                                                       │
        ├── Combined analysis (your cells + Clarke et al. mouse data)  │
        │       ├── Combined UMAP (by Leiden cluster + by strain/population)
        │       ├── Per-strain panels on the merged UMAP
        │       ├── Per-cell UMAP coordinate CSV (traceable back to plate/condition)
        │       ├── Cluster composition bar charts
        │       ├── Cluster marker genes + heatmap + violins
        │       ├── Pairwise cluster contrasts vs VAF / VRC (two-sided, up + down)
        │       ├── VAF/VRC correlation Excel  <── scored against the reference above
        │       ├── Cell identity + MCA cell type correlation Excel
        │       ├── limma batch correction applied before embedding
        │       └── Clarke cells treated as additional strain "Clarke2025"
        └── NoMHCIIFilter analysis (all sorted cells, MHCII gate bypassed,
                    depth floor only)
                ├── Per-cell lineage typing (curated marker panels)
                ├── UMAP colored by cell type
                ├── Per-cell UMAP coordinate CSV
                ├── Violins of genes of interest within focus cell types
                ├── Per-cell expression workbook + all-genes CSVs
                └── (all embeddings are per-plate batch corrected)
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

**Plate layouts currently in `data/metadata.csv`** (660 cells):

| Plate (`strain`) | `strain_group` | `batch` | Layout |
|---|---|---|---|
| B6G7 | B6G7 | batch1 | A1–A12 CD45+ MHCII+ · rest split MHCIIhi/MHCIIlo (42/42) |
| B6MHCIIGFP | B6MHCIIGFP | batch1 | A1–A12 CD45+ MHCII+ · 30 MHCIIhi / 42 MHCIIlo (84 cells) |
| NOD | NOD | batch2 | A1–A12 CD45+ MHCII+ · 42 MHCIIhi / 42 MHCIIlo |
| NOD2 | NOD | batch2 | A1–A12 CD45+ MHCII+ · 42 MHCIIhi / 42 MHCIIlo |
| NODPDL1 | NODPDL1 | batch3 | A1–B6 CD45+ MHCII+ (18) · B7–H12 CD45− MHCII+ (78) |
| NODCD31 | NODCD31 | batch4 | A1–B6 CD45− MHCII+ **CD31−** (18) · B7–H12 CD45− MHCII+ **CD31+** (78) |
| NODPDL1_2 | NODPDL1 | batch5 | A1–A12 CD45+ MHCII+ PDL1+ (12) · B1–E6 CD45− MHCII+ PDL1+ (42) · E7–H12 CD45− **MHCII−** PDL1+ (42) |

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

> **Trimmomatic heap.** `pipeline.py` invokes `trimmomatic -Xmx8g`. The bioconda
> wrapper hardcodes a 1 GB JVM heap, which OOMs
> (`Trim Stats Collector Exception -> OutOfMemoryError`) on lane-merged FASTQs —
> a 147 MB R1 was enough. Without the flag, most of a lane-merged plate fails.

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
│   ├── ClarkeOnly/                                (GSE292898 alone — builds the VAF/VRC reference)
│   │   ├── clarke_reference_profiles.csv          (the reference the Combined analysis scores against)
│   │   ├── clarke_reference_profiles.rds
│   │   ├── umap_co2_all141_by_cluster.pdf         (cell count is derived)
│   │   ├── umap_clarke_by_population.pdf
│   │   ├── umap_co2_by_sort_gate.pdf
│   │   ├── heatmap_cluster_markers.pdf
│   │   ├── violin_markers_by_cluster.pdf
│   │   └── cluster_marker_genes.xlsx
│   ├── Combined/
│   │   ├── umap_merged_by_cluster.pdf
│   │   ├── umap_merged_by_strain_population.pdf
│   │   ├── umap_per_strain_on_merged.pdf
│   │   ├── umap_coordinates_by_cell.csv           (per-cell UMAP coords + cluster/plate/condition)
│   │   ├── barplots_cluster_composition.pdf
│   │   ├── cluster_marker_genes.xlsx
│   │   ├── cluster_pairwise_contrasts.xlsx        (Combined only)
│   │   ├── heatmap_cluster_markers.pdf
│   │   ├── violin_plots_by_cluster.pdf
│   │   ├── cell_identity_merged_clusters.xlsx
│   │   ├── VAF_VRC_correlation_merged_clusters.xlsx
│   │   └── MCA_celltype_correlation_merged_clusters.xlsx
│   └── NoMHCIIFilter/
│       ├── umap_all514_by_cell_type.pdf               (cell count is derived)
│       ├── umap_coordinates_by_cell.csv               (per-cell UMAP coords + cell type/plate/condition)
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

Applied globally before any analysis. Matches Clarke et al.'s stated criterion —
*"only cells with MHC class II transcript were included in analysis"* — which is
**presence/absence in either chain**, not a quantitative floor in both:

```r
MHCII_GENES      <- c("H2-Aa", "H2-Ab1")
MHCII_MIN_COUNT  <- 1     # raw counts to call a chain "detected"
MHCII_REQUIRE_N  <- 1     # chains required: 1 = either (Clarke), 2 = both
MHCII_BYPASS_SORTED_NEG <- TRUE
```

The test runs on **raw counts**, not VST or CPM. Asking "is the transcript
present" does not need a normalised scale, and normalising first made the
threshold depend on library size.

**Sorted MHCII-negative cells bypass the gate entirely** (`MHCII_BYPASS_SORTED_NEG`).
The gate exists to confirm that cells sorted MHCII-**positive** really carry the
transcript. Applying it to cells deliberately sorted MHCII-negative tests them
against the opposite of their own gate and silently deletes the population you
sorted on purpose. Bypass is keyed on `MHCIIneg` appearing in the condition name.

Current run: **608 of 660 cells pass.** 42 sorted MHCII-negative cells bypass
the gate, 7 of which would otherwise have been dropped.

> **Two bugs this replaced, both silent.**
>
> 1. The filter required **both** chains at `log1p(CPM) >= log1p(5)` — stricter
>    than Clarke on two independent axes at once (AND vs OR, CPM≥5 vs ≥1 count).
>    It discarded 191 of 660 cells that the paper's own criterion retains.
> 2. The constant that *appeared* to control it, `MHCII_VST_MIN <- 1.0`, **was
>    never read by the filter.** The real threshold was a separate hardcoded
>    `log1p(5)` while the console printed "VST 1". Tuning the visible knob did
>    nothing and gave no indication that it had done nothing.
>
> The three constants above are now the only thing the filter consults.

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
combined matrix for all plates. Current run: per-plate filters keep 15,278–18,535
genes, and the intersection across all seven plates is **11,759 genes**. The
merged (your cells + Clarke) matrix intersects that again with Clarke's mappable
genes: **10,218 genes × 755 cells**.

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
than absent. All 18 violin genes now resolve in **7/7 plates** and in the merged
matrix. (The gene-count growth attributable to the exemption was measured at
+14 genes on the earlier 6-plate roster; it has not been re-measured since, as
it requires a deliberate A/B run with the exemption disabled.)

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
| `Combined/violin_plots_by_cluster.pdf` | **Skipped for genes missing from `merged_expr`, with the reason logged.** Cannot use the per-plate fallback: `merged_expr` is `limma::removeBatchEffect`-corrected across your cells + Clarke, so an uncorrected per-plate value would put two scales on one axis. |

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

### Per-cell UMAP coordinates (`umap_coordinates_by_cell.csv`)

Both `Combined/` and `NoMHCIIFilter/` write the embedding itself, not just the
picture: one row per cell with `UMAP1`, `UMAP2`, cluster (or cell type), and the
plate/condition fields.

Without it, a question as basic as *"which sample are those outlying points
from?"* can only be answered by re-running the whole script. With it:

```bash
awk -F, 'NR==1 || ($2>18 && $3<-8)' results/05_dge/Combined/umap_coordinates_by_cell.csv
```

> **This file was silently missing from `NoMHCIIFilter/` on first implementation.**
> The `write.csv` was placed inside the `if (NF_RUN_CLUSTER_OUTPUTS)` block,
> which is `FALSE` by default, so the call never executed and no error was
> raised — the folder simply lacked a file nobody had looked for yet. It now
> sits on the per-cell-type path, which always runs. `nf_df` also carries only
> UMAP/typing columns, so `strain` / `strain_group` / `condition` are joined back
> on from `nf_meta`; without that join the export would have been written, but
> useless for the one question it exists to answer.

### Cluster marker genes (Wilcoxon rank-sum)

One-vs-rest Wilcoxon on VST matrix. Pre-filter: log2FC ≥ 0.5. BH correction.
Top 100 per cluster in Excel, top 15 per cluster in heatmap.

**One-sided and upregulated-only** (`alternative="greater"`, candidates pre-filtered
to `log2FC >= MIN_LOG2FC`). It answers "what marks this cluster against all others
pooled" and **structurally cannot report depletion.** For the complementary question
see the next section.

### Pairwise cluster contrasts vs VAF / VRC (Combined only)

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

> ## ⚠ The "VAF" reference is not fibroblasts. Do not use it as one.
>
> GSE292898 ships **no cell-type labels** — the `cell type` column is empty for
> all 602 samples; only the sort gate is published. Every VAF/VRC/CD45pos
> population in this pipeline is therefore **inferred**, not published.
>
> The reference is now the **ClarkeOnly cluster centroids** (see below), and the
> conclusion survived the change — it got *stronger*. Mean VST in the three
> centroids actually used by the current run:
>
> | Marker class | CD45pos | **VAF** | VRC |
> |---|---|---|---|
> | Immune — `Ptprc` | **7.87** | −2.51 | −3.72 |
> | Endothelial — `Pecam1` / `Cdh5` / `Plvap` | −2.24 / −4.37 / −3.69 | −1.04 / −4.07 / −1.57 | **9.91 / 8.74 / 13.91** |
> | Fibroblast — `Col1a1` / `Pdgfra` / `Dcn` / `Lum` | −3.85 / −4.02 / −3.85 / −3.79 | **−2.86 / −3.34 / −3.03 / −3.25** | −4.09 / −4.37 / −3.97 / −4.09 |
> | Endocrine — `Gcg` / `Ins2` / `Pcsk2` / `Chga` | 0.65 / 6.20 / −3.16 / −2.98 | **9.74 / 9.03 / 7.97 / 7.59** | −0.51 / 6.14 / −3.80 / −3.78 |
>
> **Every fibroblast marker is negative in the VAF centroid. Every endocrine
> marker is strongly positive.** Its top Wilcoxon marker is `Pcsk2`, and `Col1a1`
> is detected in 5 of its 45 cells. `CD45pos` (Ptprc 7.9) and `VRC`
> (Plvap 13.9, Pecam1 9.9) are by contrast cleanly identified and can be read at
> face value — the problem is specific to VAF.
>
> **This is not a clustering artifact.** Across all 95 Clarke CD45− MHCII+
> cells, only **8 carry ≥3 fibroblast markers** (3 carry ≥5), versus 45 carrying
> ≥3 endothelial markers. The cells are not in the dataset. No choice of
> resolution, algorithm, or batch correction recovers a population with ~8
> members. Clarke et al. themselves write that they could not determine a
> convincing embryologic origin for this population.
>
> **How to read a VAF match:** "this cluster resembles Clarke's CD45−
> non-endothelial cells, which are predominantly endocrine." **Not** "this
> cluster is fibroblastic." `cluster_pairwise_contrasts.xlsx` baselines every
> contrast on the same population and inherits the caveat. For fibroblast
> questions use your own cells — NoMHCIIFilter's per-cell typing finds 25–27.

#### Confidence gating on the best match

`best_match` used to be reported with no margin and no floor, so a statistical
tie and a genuine match were formatted identically. A call is now only reported
when it clears **both** tests, and is otherwise blanked to `unresolved` with a
`why_unresolved` reason:

| Constant | Default | Test |
|---|---|---|
| `MIN_CORR_R` | 0.50 | the winning Pearson r must reach this |
| `MIN_CORR_MARGIN` | 0.05 | ...and beat the runner-up by at least this |

The Summary sheet gains `confident`, `runner_up`, `runner_up_r`, `margin`, and
`why_unresolved`; unresolved rows are highlighted. Per-cluster sheets still
carry every raw correlation, so nothing is hidden. The thresholds are round
numbers chosen to separate observed cases, not derived — adjust if they flag
something you consider real.


Two distinct steps that are easy to conflate. **Marker genes are not used to
decide a cluster's identity.**

**1. Cluster → identity match: genome-wide, no marker list.** Each cluster's mean
VST profile is correlated (Pearson *and* Spearman) against the mean profiles of
VAF, VRC, and CD45pos, across **all shared genes — 10,215 in the merged
analysis** (`n_genes` column in the output workbook). Both up- and
down-regulation contribute; nothing is restricted to a marker panel. The
correlations discriminate rather than saturating — in the current merged run:

| Cluster | VAF | VRC | CD45pos | Call |
|---|---|---|---|---|
| 1 | **0.808** | 0.570 | 0.603 | VAF |
| 2 | 0.652 | 0.691 | **0.819** | CD45pos |
| 3 | 0.215 | 0.400 | 0.378 | *unresolved* |
| 4 | 0.626 | **0.850** | 0.633 | VRC |
| 5 | **0.612** | 0.330 | 0.355 | VAF |

(Pearson, n_genes = 10,215.) Cluster 3 fails both confidence tests and is
correctly left uncalled.

**2. Where the reference profiles come from: the ClarkeOnly centroids.**
The reference is built by running Clarke's cells through **this same pipeline,
alone** — the `ClarkeOnly` analysis — and taking its three cluster centroids:

- the **split** is unsupervised: Leiden at `CLARKE_LEIDEN_RESOLUTION` on
  all-gene PCs, using no marker information at all. It independently resolves
  **exactly 3 clusters**, matching the paper's three populations;
- the **naming** uses markers *only* to assign the three labels, strongest
  signal first: highest mean `Ptprc` → CD45pos, then highest `Pecam1` → VRC,
  remainder → VAF. The separation is unambiguous (Ptprc 7.87 vs −2.51/−3.72;
  Pecam1 9.91 vs −2.24/−1.04);
- if ClarkeOnly ever resolves something other than 3 clusters, the extras are
  labelled `other_*`, **excluded** from the reference, and a warning is printed
  — it does not silently fold them into a published population.

The centroids are written to `ClarkeOnly/clarke_reference_profiles.csv` (and
`.rds`), so the reference is inspectable rather than implicit. Switching
`RUN_CLARKEONLY` off falls back to the legacy k-means reconstruction and says so
loudly in the log; the two are on **different scales** (VST vs log1p(CPM)), so r
values are not comparable between them.

> **Two VAF/VRC definitions coexist in this script. Know which one you are reading.**
>
> | Consumer | Definition used |
> |---|---|
> | `VAF_VRC_correlation_*.xlsx` | **ClarkeOnly centroids** (VST) |
> | Clarke cells' `VAF`/`VRC` labels in the merged UMAP and bar charts | k-means (k=2) split |
> | `cluster_pairwise_contrasts.xlsx` baseline cluster | k-means split, via modal cluster |
>
> The k-means path still runs because it assigns a label to *individual* Clarke
> cells, which a centroid cannot do. The two do not perfectly agree: the modal
> merged cluster for k-means "VAF" cells is cluster 1 at only **46.3%**
> (31/67 cells), against 98.3% (57/58) for VRC. Treat the pairwise-contrast VAF
> baseline as the weaker of the two.

**Legacy k-means path (fallback and per-cell labelling).** PCA on the top 500
variable genes, then k-means with k=2; the **naming** uses the marker panels
(Col1a1/Col1a2/Timp3/Spp1/Thy1/Pdpn for VAF; Pecam1/Eng/Cdh5/Kdr/Tie1/Vwf for
VRC) *only* to decide which of the two clusters gets called "VAF". Orientation
scores both panels (VAF markers minus VRC markers) so a cluster high in both
cannot win by default.

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
coloring and tabulation. (Note the UMAP clusters are **Leiden**, as are the
ClarkeOnly clusters the reference profiles come from; k-means appears only in the
per-cell Clarke labeling above.)

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

`VAF_VRC_correlation_merged_clusters.xlsx` (Combined) carries an
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
is included or excluded here on the basis of what it expresses. 660 sorted →
**514 analyzed**, 146 removed.

The method is deliberately **identical** to the `combined_plots` UMAP so the two
are directly comparable — per-plate low-count filter (violin-gene exemption
included) → per-plate DESeq2 with size factors from that plate's own reference
wells → VST → intersect genes across plates → top `N_HVG` by variance → PCA →
UMAP for layout, plus a separate all-gene PCA → kNN → Leiden for the cluster
assignment. Same tuning constants, same seed. **The only difference is which
cells go in.**

| | MHCII-gated path | `NoMHCIIFilter` |
|---|---|---|
| Gate | MHCII transcript present | sequencing depth only |
| Cells | 608 | **514** |
| Genes in intersection | 11,759 | **11,773** |

Note the two gates are **not nested**: the MHCII gate keeps 608 of 660 cells and
the depth floor keeps 514, but neither is a subset of the other — a cell can be
MHCII-positive and shallow, or deep and MHCII-negative. The `combined_plots`
columns previously quoted here came from a run with `RUN_COMBINED_PLOTS <- TRUE`
and an older plate roster; that flag is now `FALSE`, so those figures are not
re-measurable from the current log and have been removed rather than carried
forward.

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

`Acinar_ductal` sits at 8/9 because **`Cpa1` has zero counts across all 514
cells**, so z-scoring produces `NaN` and it is dropped. The gene is present in
the count matrix and the annotation — it is simply not expressed in this sort,
which is unsurprising for islet preparations with little exocrine carryover.

That is a harmless shortfall, but it is exactly the kind of thing to check
before trusting a call: a panel silently shrinking is how the `nf_expr` bug
above went unnoticed. Several acinar markers are near-absent here — `Amy2a5`
in 4 of 514 cells, `Prss2` in 4, `Ctrb1` in 10 — so the `Acinar_ductal` label
rests largely on `Krt18` (166 cells) and `Spp1` (110), both of which are broader
than the exocrine compartment. **Treat the 24 acinar/ductal calls with more
caution than the endothelial or hematopoietic ones.**

Current per-cell type calls (514 cells):

| Cell type | n |
|---|---|
| `Endocrine_islet` | 171 |
| `Endothelial` | 164 |
| `Hematopoietic` | 117 |
| `Fibroblast` | **27** |
| `Acinar_ductal` | 24 |
| `Pericyte_SMC` | 11 |

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
VRC 46, 4 merged clusters, all unchanged. *(Those figures are the state of the
6-plate roster on which the check was run; they are a record of the check, not
current counts. The additivity property they verify is structural and still
holds.)*

#### Why the depth floor is not optional here

Removing the MHCII filter also removes the depth screen it was incidentally
performing. Running with no floor at all (`NOFILT_MIN_GENES <- 0`) was tried
first and produced two concrete failures:

1. **A pure artifact cluster.** All 564 cells (the roster at the time) gave
   **6** clusters, and cluster 6 (n=29) sat alone at UMAP1 ≈ 17, far off the main
   manifold — shallow cells clustering by library size rather than by biology.
   The shallowest cell in the dataset detects **51 genes on 101 reads**. With the
   floor applied the island disappears and the structure resolves cleanly. This
   experiment has not been repeated on the current 7-plate roster; the floor has
   been left on since.
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

### Combined analysis

Merges your MHCII-filtered VST data with Clarke et al. 2025 mouse cells
(GSE292898, 147/188 cells passing the MHCII gate). Clarke cells are processed
through independent DESeq2 VST then **limma batch correction**
(`removeBatchEffect`) is applied to the merged matrix before embedding.
Clarke cells are labeled as `Clarke2025 CD45pos`, `Clarke2025 VAF`, or
`Clarke2025 VRC`. All combined_plots outputs are reproduced for this
merged dataset. The other 5 output folders are completely unaffected.

---

## ClarkeOnly — the Clarke/Don data on its own

`results/05_dge/ClarkeOnly/` analyses GSE292898 **alone**: no cells of yours, no
cross-plate gene intersection, no k-means VAF/VRC reconstruction, no
batch correction (single dataset). It also works in Clarke's own gene-symbol
space rather than mapping through your Ensembl IDs, because mapping first would
silently drop genes absent from your plates.

Outputs: UMAP by Leiden cluster, UMAP by published sort gate, cluster marker
heatmap, marker workbook, and violins of `CLARKE_VIOLIN_GENES`
(`Ptprc`, `Col1a1`, `Col1a2`, `Pecam1`) per cluster with per-cluster detection
counts on the axis labels.

### It is the source of the VAF/VRC reference

ClarkeOnly is **not a side-analysis**. Its three cluster centroids are what the
Combined analysis correlates every cluster against (see
[VAF/VRC correlation](#vafvrc-correlation-clarke-et-al-2025-gse292898)). It runs
*before* the Combined block for exactly that reason.

This is why it is **not frozen to a fixed resolution**: it is clustered by the
same auto-selected procedure as everything else. Pinning it to a hand-picked
resolution that reproduces the paper, then using those centroids as the
yardstick for your own cells, would be circular — it would guarantee the
reference matched the paper's description regardless of what the data said.

### Reproducing the paper's 3 clusters

The paper reports **3 clusters** — fibroblastic (green), vascular (blue), CD45+
(purple) — from **142** retained cells. **The pipeline's defaults now reproduce
this**, after the MHC-II gate was corrected to match Clarke's criterion:

| | This pipeline | Clarke et al. |
|---|---|---|
| MHC-II | **required**, either chain ≥1 count | **required** ("only cells with MHC class II transcript") |
| Cell filter | ≥500 **genes detected** | ≥1000 **total counts** |
| Normalisation | DESeq2 VST | Seurat `LogNormalize` + scale |
| Clustering | Leiden on 80%-variance PCs (auto-selected) | Seurat/Louvain on ~20 PCs |
| **Cells retained** | **141** | **142** |
| **Clusters** | **3** | **3** |

188 cells → 147 pass the MHC-II gate → **141** clear the ≥500-gene depth floor.
Clustering on 20,326 genes resolves 3 clusters with no tuning:

| Cluster | n | Sort gate (neg/pos) | Top marker | `Ptprc` | `Pecam1` | `Col1a1` | Label |
|---|---|---|---|---|---|---|---|
| 1 | 37 | 19 / 18 | `Ctss` | **30/37** | 13/37 | 2/37 | CD45pos (purple) |
| 2 | 45 | 41 / 4 | `Pcsk2` | 7/45 | 20/45 | 5/45 | "VAF" (green) |
| 3 | 59 | 59 / 0 | `Plvap` | 4/59 | **56/59** | 1/59 | VRC (blue) |

Clusters 1 and 3 are unambiguous — `Ptprc` in 81% of cluster 1, `Pecam1` in 95%
of cluster 3, and cluster 3 contains **zero** CD45+ sorted cells.

> **Cluster 2 is not fibroblastic, whatever the paper calls it.** Its top marker
> is `Pcsk2` (endocrine), `Col1a1` appears in 5 of 45 cells (11%), and `Pecam1`
> (20/45) is *more* common than any collagen gene. The label reflects "neither
> immune nor clearly endothelial", not collagen positivity. This is the same
> population the VAF reference is built from — see the warning box above.

> **Deposit/paper discrepancy.** The paper states 114 CD45− and 50 CD45+ cells
> sorted (98 + 44 analysed). The GEO deposit contains **164 CD45neg and 24
> CD45pos**. We recover nearly the right *total* (141 vs 142) but a different
> split. The deposit's composition does not match the figure's description, and
> this is a property of the deposit, not of the filtering.

---

## Per-plate batch correction

`BATCH_CORRECT_BY_PLATE <- TRUE` applies `limma::removeBatchEffect` with **plate
as the batch** to the merged and NoMHCIIFilter embeddings.

**Why it is needed.** Size factors are estimated from each plate's *own*
reference wells, so every plate is scaled to a different baseline and the
normalisation itself imposes a plate-specific offset. Measured on the CD45+
reference cells — nominally the same population across plates, so any separation
between them is technical (measured on the 6-plate roster; the correction itself
now spans 8 batches including Clarke2025):

| | Raw log-CPM | Per-plate VST (what clustering consumes) | After correction |
|---|---|---|---|
| PC1 | 1.3% | **20.5%** | 2.6% |
| PC1–10 combined | 13.3% | **23.3%** | **6.8%** |

The normalisation *amplifies* plate structure into the dominant components.
Measuring this on raw log-CPM instead of VST gives a falsely reassuring answer —
that mistake was made once in this project and led to the wrong conclusion.

The visible symptom was NODPDL1_2 landing **47/47 in a single cluster**, both
subgates together. After correction it spreads 43% / 21% / 34%.

**Condition labels must be harmonised first.** 6 of 9 raw conditions exist on
exactly one plate (all three of NODPDL1_2's, both of NODCD31's), which makes
`~ plate + condition` rank deficient — 13 of 15 columns. `harmonise_condition()`
collapses only the single-plate labels into the nearest label spanning several
plates, restoring full rank (11 of 11). MHCIIhi/MHCIIlo span 4 plates and are
deliberately kept.

> **What this does NOT remove.** The batch offset is one scalar per plate per
> gene applied to every cell on that plate. Differences *within* a plate — CD31+
> vs CD31−, MHCIIhi vs MHCIIlo, VAF vs VRC — shift identically and are fully
> preserved. Only how a plate sits relative to other plates is removed.

`apply_plate_correction()` re-checks the design rank at run time and **skips
correction with a logged message** rather than returning a half-corrected
matrix. Constant genes are held out and re-attached unchanged.

> **Known unexplained failure.** After correction, `prcomp` fails on the real
> matrix with `error code 1 from Lapack routine 'dgesdd'` — while succeeding on
> random matrices of identical dimensions, with every value finite and the range
> unremarkable ([-11, 26]). `run_umap()` falls back to an eigendecomposition of
> the Gram matrix (`x x' = U D^2 U'`, scores `U D`), which is the same PCA via
> `dsyevr` instead of `dgesdd`. The results are sound but the root cause was
> never found. The offending matrix is dumped to `tempdir()` on failure.

---

## Which analyses run

`scripts/per_strain_plots.R` has five switches near the top:

```r
RUN_PER_STRAIN_PLOTS <- FALSE   # results/05_dge/<STRAIN>_plots/
RUN_COMBINED_PLOTS   <- FALSE   # results/05_dge/combined_plots/
RUN_CLARKEONLY       <- TRUE    # results/05_dge/ClarkeOnly/   (runs FIRST)
RUN_VAF_MERGED       <- TRUE    # results/05_dge/Combined/
RUN_NOMHCIIFILTER    <- TRUE    # results/05_dge/NoMHCIIFilter/
```

Turning the first two off is a real but modest win, and it is worth being clear
about why: most of the remaining time is work the flags *cannot* skip — the
seven per-plate DESeq2/VST fits are mandatory (see the dependency table below),
and the Wilcoxon loops in the enabled folders are themselves expensive (the
merged pairwise contrasts run eight two-sided comparisons over ~10,200 genes).
Getting substantially below this needs caching the per-plate VST matrices to
disk, not more flags.

> **`RUN_CLARKEONLY` is not independent of `RUN_VAF_MERGED`.** ClarkeOnly builds
> the VAF/VRC/CD45pos reference profiles that Combined correlates against, which
> is why it runs first. Switching it off does not just drop a folder — it
> silently downgrades Combined to the legacy k-means reference, on a different
> scale (log1p(CPM) rather than VST). The log says so loudly when this happens;
> the workbook's Notes sheet records which reference was actually used.

Nothing is deleted — flip a flag back to `TRUE` to restore that output exactly
as before.

> **`RUN_COMBINED_PLOTS = FALSE` skips the `combined_plots` folder, not the
> combined matrix it is named after.** Several things inside those sections are
> consumed by the analyses that are still on, so they run regardless of the
> flags:
>
> | Still runs | Why |
> |---|---|
> | The per-strain DESeq2/VST loop | Builds `all_expr_list` / `all_meta_list` → `combined_expr` + `combined_meta`, which Combined consumes. Only the loop's plots and spreadsheets are skipped, never its arithmetic. |
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
| `MHCII_GENES` | H2-Aa, H2-Ab1 | MHC-II chains the gate tests |
| `MHCII_MIN_COUNT` | 1 | Raw counts for a chain to count as detected |
| `MHCII_REQUIRE_N` | 1 | Chains required: 1 = either (Clarke), 2 = both |
| `MHCII_BYPASS_SORTED_NEG` | TRUE | Cells sorted MHCII-negative skip the gate entirely |
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
| `BATCH_CORRECT_BY_PLATE` | TRUE | Per-plate `removeBatchEffect` on the merged and NoMHCIIFilter embeddings |
| `MIN_CORR_R` | 0.50 | VAF/VRC: minimum winning Pearson r before a call is reported |
| `MIN_CORR_MARGIN` | 0.05 | VAF/VRC: minimum gap to the runner-up before a call is reported |
| `RUN_CLARKEONLY` | TRUE | `results/05_dge/ClarkeOnly/` — Clarke data analysed alone. **Also builds the VAF/VRC reference the Combined analysis uses**, so switching it off silently downgrades Combined to the legacy k-means reference |
| `CLARKE_LEIDEN_RESOLUTION` | 1.0 | Leiden granularity for ClarkeOnly only; `LEIDEN_RESOLUTION` governs everything else |
| `CLARKE_VIOLIN_GENES` | Ptprc, Col1a1, Col1a2, Pecam1 | Genes given a per-cluster violin in ClarkeOnly |
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
  batch3 (L001+L002 merged); NODCD31 = batch4; NODPDL1_2 = batch5 (L001+L002 merged)
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
