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
        ├── MHCII expression filter (H2-Aa & H2-Ab1, applied globally)
        ├── Per-strain DESeq2 (independent per plate, normalized to A1-A12)
        │       ├── Volcano plots (3 per strain)
        │       ├── Violin plots (3 per strain, top 20 DEGs each)
        │       ├── Expression matrix CSV (mean VST + upregulation rankings)
        │       ├── Per-strain UMAP (by cluster + by population)
        │       └── Cell identity Excel (CellMarker 2.0 gene set scoring)
        ├── Combined analysis (all MHCII-filtered cells, your 4 strains)
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
`cell_id`, `strain`, `condition`, `batch`

See `data/metadata_example.csv` for format reference.

### 4. Configure

Edit `config.yaml`:
- Set `FASTQ_DIR` to your folder of `.fastq.gz` files
- Set `METADATA_FILE` to your metadata CSV

### 5. Run alignment and counting

```bash
python pipeline.py --config config.yaml
```

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
│   ├── NOD_plots/                            (12 files)
│   │   ├── expression_matrix_NOD.csv
│   │   ├── volcano_NOD_*.pdf                 (3 files)
│   │   ├── violins_NOD_*.pdf                 (3 files)
│   │   ├── umap_NOD_by_cluster.pdf
│   │   ├── umap_NOD_by_population.pdf
│   │   └── cell_identity_NOD_clusters.xlsx
│   ├── B6G7_plots/                           (same 12 files)
│   ├── B6MHCIIGFP_plots/                     (same 12 files)
│   ├── NODPDL1_plots/                        (same 12 files)
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
│       ├── VAF_VRC_correlation_merged_clusters.xlsx
│       └── MCA_celltype_correlation_merged_clusters.xlsx
└── qc_summary/
    ├── cell_qc_metrics.csv
    ├── cell_qc_barplots.pdf
    └── metadata_qc_pass.csv
```

---

## per_strain_plots.R — detailed documentation

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
estimated using only the **CD45pos_MHCIIpos (A1-A12) wells** of each plate.
No cross-plate pooling occurs.

### Per-strain DESeq2

Design: `~ condition`. VST via `varianceStabilizingTransformation()`.

### Per-strain contrasts

| Comparison | Volcano | Violins |
|---|---|---|
| CD45+ MHCIIpos vs CD45− MHCIIhi | ✓ | ✓ |
| CD45+ MHCIIpos vs CD45− MHCIIlo | ✓ | ✓ |
| CD45− MHCIIhi vs CD45− MHCIIlo  | ✓ | ✓ |

Each violin PDF: top 10 up + top 10 down DEGs, gene symbols, jittered VST points.

### Expression matrix columns

| Column | Description |
|---|---|
| `ensembl_id` | Ensembl gene ID |
| `gene_symbol` | Common gene name (GENCODE vM33) |
| `mean_VST_CD45pos_MHCIIpos` | Mean VST expression |
| `mean_VST_CD45neg_MHCIIhi` | Mean VST expression |
| `mean_VST_CD45neg_MHCIIlo` | Mean VST expression |
| `rank_upregulated_MHCIIhi` | Upregulation rank vs CD45pos (1 = most upregulated) |
| `rank_upregulated_MHCIIlo` | Upregulation rank vs CD45pos (1 = most upregulated) |

Ranks use combined score `log2FC × -log10(padj)`.

### Per-strain UMAPs

- `umap_{STRAIN}_by_cluster.pdf` — Leiden clusters (all-gene PCA for clustering, HVG UMAP for visualization)
- `umap_{STRAIN}_by_population.pdf` — colored by CD45+/MHCIIhi/MHCIIlo condition

### Cell identity scoring (per-strain and combined)

Fisher's exact test + Jaccard similarity against CellMarker 2.0 mouse reference
(downloaded at runtime, cached). Falls back to 20-cell-type built-in set if
download fails. When Wilcoxon markers are insufficient, falls back to top
expressed genes.

### Combined UMAP (your 4 strains)

**UMAP embedding:** top 2000 HVG PCA → 80% variance PCs → `uwot::umap()`

**Leiden clustering:** all-gene PCA → 80% variance PCs → kNN graph → Leiden

**Per-strain panels:** 4-panel PDF showing each strain's cells on shared
comprehensive UMAP coordinates, colored by comprehensive cluster assignment.

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

- **Strains:** NOD, B6G7, B6MHCIIGFP, NODPDL1
- **Conditions:** `CD45pos_MHCIIpos` (A1–A12), `CD45neg_MHCIIhi` (B1–E6), `CD45neg_MHCIIlo` (E7–H12)
- **Plate quirks:** B6MHCIIGFP wells H1–H12 empty; NODPDL1 layout physically flipped but metadata labels are biologically correct
- **Batches:** B6G7 and B6MHCIIGFP = L001; NOD and NODPDL1 = L002
- **Normalization reference:** CD45pos_MHCIIpos (A1–A12), per plate
- **DESeq2 design:** `~ condition` (no batch term; each strain is single-batch)
- **VST method:** `varianceStabilizingTransformation()` (bypasses sample-size check in `vst()`)
- **Aligner:** HISAT2 (switched from STAR due to Apple Silicon RAM constraints)
- **Annotation:** GENCODE vM33 primary assembly GTF
- **Clarke et al. 2025:** Cell Reports 44, 116189. GEO: GSE292898. VAFs are CD45neg MHC-II+ fibroblastic cells from NOD mouse pancreatic islets.
