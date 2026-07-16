#!/usr/bin/env bash
# ==============================================================================
# setup.sh — Download mm10/GRCm39 reference and build HISAT2 index
# Run once before the pipeline. Takes ~20-30 min on first run.
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

REF_DIR="reference"
HISAT2_INDEX_DIR="reference/hisat2_index"
THREADS=8

echo "=== SMARTseq HT Pipeline — Reference Setup ==="
echo "Working directory: $SCRIPT_DIR"
echo ""

# ── Create directory structure ─────────────────────────────────────────────────
mkdir -p "$REF_DIR" "$HISAT2_INDEX_DIR" data results logs

# ── Download GRCm39 genome (GENCODE M33) ──────────────────────────────────────
GENOME_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/GRCm39.primary_assembly.genome.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.primary_assembly.annotation.gtf.gz"

GENOME_FA_GZ="$REF_DIR/GRCm39.primary_assembly.genome.fa.gz"
GENOME_FA="$REF_DIR/GRCm39.primary_assembly.genome.fa"
GTF_GZ="$REF_DIR/gencode.vM33.primary_assembly.annotation.gtf.gz"
GTF_FILE="$REF_DIR/gencode.vM33.primary_assembly.annotation.gtf"

if [[ ! -f "$GENOME_FA_GZ" ]]; then
    echo "[1/4] Downloading GRCm39 genome FASTA (~800 MB)..."
    curl -L --progress-bar -o "$GENOME_FA_GZ" "$GENOME_URL"
else
    echo "[1/4] Genome FASTA already present, skipping download."
fi

if [[ ! -f "$GTF_GZ" ]]; then
    echo "[2/4] Downloading GENCODE M33 GTF (~30 MB)..."
    curl -L --progress-bar -o "$GTF_GZ" "$GTF_URL"
else
    echo "[2/4] GTF already present, skipping download."
fi

# ── Decompress genome + GTF (HISAT2 build needs uncompressed FASTA) ───────────
if [[ ! -f "$GENOME_FA" ]]; then
    echo "[3/4] Decompressing genome FASTA..."
    gunzip -k "$GENOME_FA_GZ"
else
    echo "[3/4] Decompressed genome FASTA already present."
fi

if [[ ! -f "$GTF_FILE" ]]; then
    echo "[3/4] Decompressing GTF..."
    gunzip -k "$GTF_GZ"
else
    echo "[3/4] Decompressed GTF already present."
fi

# ── Build HISAT2 genome index ──────────────────────────────────────────────────
if [[ ! -f "$HISAT2_INDEX_DIR/GRCm39.1.ht2" ]]; then
    echo "[4/4] Building HISAT2 index (this takes ~15-20 min)..."
    echo "      Threads: $THREADS"
    hisat2-build -p "$THREADS" "$GENOME_FA" "$HISAT2_INDEX_DIR/GRCm39" \
        2>&1 | tee logs/hisat2_index_build.log
    echo "[4/4] HISAT2 index built successfully."
else
    echo "[4/4] HISAT2 index already exists, skipping."
fi

echo ""
echo "=== Setup complete! ==="
echo ""
echo "Next steps:"
echo "  1. Point FASTQ_DIR in config.yaml at your folder of .fastq.gz files."
echo "     Expected naming: {STRAIN}_{WELL}_{SAMPLE}_{LANE}_R1_001.fastq.gz"
echo "                      {STRAIN}_{WELL}_{SAMPLE}_{LANE}_R2_001.fastq.gz"
echo "     Example: NOD_A1_S1_L002_R1_001.fastq.gz"
echo ""
echo "  2. Create your metadata file: data/metadata.csv"
echo "     Required columns: cell_id, strain, condition, batch"
echo "     See metadata_example.csv for format reference."
echo ""
echo "  3. Edit config.yaml if needed (threads, thresholds, design formula)"
echo ""
echo "  4. Run the pipeline:"
echo "     python pipeline.py --config config.yaml"
