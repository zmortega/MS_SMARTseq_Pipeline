#!/usr/bin/env bash
# Removes .done flags for any cell with an empty BAM so pipeline will re-align them
echo "Scanning for empty BAMs..."
count=0
for bam in results/03_aligned/*/*_Aligned.sortedByCoord.out.bam; do
    reads=$(samtools view -c "$bam" 2>/dev/null)
    if [[ "$reads" == "0" ]]; then
        cell_dir=$(dirname "$bam")
        rm -f "$cell_dir/.done"
        rm -f "$bam" "$bam.bai"
        count=$((count + 1))
    fi
done
echo "Cleared $count empty BAM(s). Run pipeline with --resume to re-align."
