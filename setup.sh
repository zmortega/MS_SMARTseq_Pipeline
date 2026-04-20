#!/usr/bin/env bash
# ==============================================================================
# setup.sh — Download mm10/GRCm39 reference and build STAR index
# Run once before the pipeline. Takes ~45 min on first run.
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

REF_DIR="reference"
STAR_INDEX="reference/star_index"
THREADS=8

echo "=== SMARTseq HT Pipeline — Reference Setup ==="
echo "Working directory: $SCRIPT_DIR"
echo ""

# ── Create directory structure ─────────────────────────────────────────────────
mkdir -p "$REF_DIR" "$STAR_INDEX" data/fastq data/metadata results logs

# ── Download GRCm39 genome (GENCODE M33) ──────────────────────────────────────
GENOME_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/GRCm39.primary_assembly.genome.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.primary_assembly.annotation.gtf.gz"

GENOME_FA="$REF_DIR/GRCm39.primary_assembly.genome.fa.gz"
GTF_GZ="$REF_DIR/gencode.vM33.primary_assembly.annotation.gtf.gz"
GTF_FILE="$REF_DIR/gencode.vM33.primary_assembly.annotation.gtf"

if [[ ! -f "$GENOME_FA" ]]; then
    echo "[1/4] Downloading GRCm39 genome FASTA (~800 MB)..."
    curl -L --progress-bar -o "$GENOME_FA" "$GENOME_URL"
else
    echo "[1/4] Genome FASTA already present, skipping download."
fi

if [[ ! -f "$GTF_GZ" ]]; then
    echo "[2/4] Downloading GENCODE M33 GTF (~30 MB)..."
    curl -L --progress-bar -o "$GTF_GZ" "$GTF_URL"
else
    echo "[2/4] GTF already present, skipping download."
fi

# ── Decompress GTF (STAR needs uncompressed GTF) ──────────────────────────────
if [[ ! -f "$GTF_FILE" ]]; then
    echo "[3/4] Decompressing GTF..."
    gunzip -k "$GTF_GZ"
else
    echo "[3/4] Decompressed GTF already present."
fi

# ── Build STAR genome index ────────────────────────────────────────────────────
if [[ ! -f "$STAR_INDEX/SA" ]]; then
    echo "[4/4] Building STAR index (this takes ~30 min and ~30 GB RAM)..."
    echo "      Threads: $THREADS"
    STAR \
        --runMode genomeGenerate \
        --genomeDir "$STAR_INDEX" \
        --genomeFastaFiles <(gunzip -c "$GENOME_FA") \
        --sjdbGTFfile "$GTF_FILE" \
        --sjdbOverhang 149 \
        --runThreadN "$THREADS" \
        --genomeSAindexNbases 14 \
        2>&1 | tee logs/star_index_build.log
    echo "[4/4] STAR index built successfully."
else
    echo "[4/4] STAR index already exists, skipping."
fi

echo ""
echo "=== Setup complete! ==="
echo ""
echo "Next steps:"
echo "  1. Place your .fastq.gz files in: data/fastq/"
echo "     Expected naming: {CELL_ID}_R1.fastq.gz / {CELL_ID}_R2.fastq.gz"
echo ""
echo "  2. Create your metadata file: data/metadata.csv"
echo "     Required columns: cell_id, condition"
echo "     Example:  PLATE1_A01,control"
echo "               PLATE1_B01,treatment"
echo ""
echo "  3. Edit config.yaml if needed (threads, thresholds, design formula)"
echo ""
echo "  4. Run the pipeline:"
echo "     python pipeline.py --config config.yaml"
