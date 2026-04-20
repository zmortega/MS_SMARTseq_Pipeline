#!/usr/bin/env python3
"""Fix counts_clean.txt column names by extracting cell IDs from BAM paths."""
from pathlib import Path

counts_file = Path("results/04_counts/raw_counts.txt")
cleaned = Path("results/04_counts/counts_clean.txt")

with open(counts_file) as f_in, open(cleaned, "w") as f_out:
    for i, line in enumerate(f_in):
        if i == 0:
            f_out.write(line)
            continue
        if i == 1:
            parts = line.rstrip("\n").split("\t")
            fixed_cols = []
            for col in parts[6:]:
                bam_name = col.split("/")[-1]
                cell_id = bam_name.replace("_Aligned.sortedByCoord.out.bam", "")
                fixed_cols.append(cell_id)
            fixed = parts[:6] + fixed_cols
            f_out.write("\t".join(fixed) + "\n")
            continue
        f_out.write(line)

print("Done! Verifying first few column names:")
with open(cleaned) as f:
    for i, line in enumerate(f):
        if i == 1:
            cols = line.split("\t")
            print("Cols 1-5:", cols[0])
            print("First 5 cell cols:", cols[6:11])
            break
