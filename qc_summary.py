#!/usr/bin/env python3
"""
qc_summary.py — Parse STAR logs and featureCounts output to flag low-quality cells.
Run after the pipeline completes (or after step 4).

Usage:
    python qc_summary.py --config config.yaml
"""

import argparse
import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import yaml


def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def parse_star_log(log_file: Path) -> dict:
    """Extract key mapping stats from STAR Log.final.out"""
    stats = {}
    patterns = {
        "total_reads":       r"Number of input reads \|\s+([\d]+)",
        "uniquely_mapped":   r"Uniquely mapped reads number \|\s+([\d]+)",
        "mapping_rate":      r"Uniquely mapped reads % \|\s+([\d.]+)%",
        "multi_mapped":      r"Number of reads mapped to multiple loci \|\s+([\d]+)",
        "unmapped_too_short":r"% of reads unmapped: too short \|\s+([\d.]+)%",
    }
    text = log_file.read_text()
    for key, pat in patterns.items():
        m = re.search(pat, text)
        if m:
            val = m.group(1)
            stats[key] = float(val) if "." in val or key == "mapping_rate" else int(val)
    return stats


def parse_featurecounts(counts_file: Path) -> pd.Series:
    """Return per-cell total assigned counts."""
    df = pd.read_table(counts_file, comment="#", index_col=0)
    expr = df.iloc[:, 5:]
    genes_detected = (expr > 0).sum()
    total_counts = expr.sum()
    return genes_detected, total_counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    cfg = load_config(args.config)

    out_dir = Path(cfg["OUTPUT_DIR"])
    align_dir = out_dir / "03_aligned"
    count_dir = out_dir / "04_counts"
    qc_out = out_dir / "qc_summary"
    qc_out.mkdir(parents=True, exist_ok=True)

    # Thresholds from config
    min_genes = cfg.get("MIN_GENES_DETECTED", 500)
    max_genes = cfg.get("MAX_GENES_DETECTED", 8000)
    min_reads = cfg.get("MIN_READS_MAPPED", 500000)
    min_rate  = cfg.get("MIN_MAPPING_RATE", 0.60)

    # ── Parse STAR logs ────────────────────────────────────────────────────────
    records = []
    for star_log in sorted(align_dir.glob("*/*Log.final.out")):
        cell_id = star_log.parent.name
        stats = parse_star_log(star_log)
        stats["cell_id"] = cell_id
        records.append(stats)

    if not records:
        print(f"No STAR logs found in {align_dir}. Run alignment first.")
        return

    df_star = pd.DataFrame(records).set_index("cell_id")

    # ── Parse featureCounts ────────────────────────────────────────────────────
    counts_clean = count_dir / "counts_clean.txt"
    if counts_clean.exists():
        genes_detected, total_counts = parse_featurecounts(counts_clean)
        df_star["genes_detected"] = genes_detected
        df_star["total_assigned_counts"] = total_counts
    else:
        print("Warning: counts_clean.txt not found — skipping gene/count metrics.")
        df_star["genes_detected"] = pd.NA
        df_star["total_assigned_counts"] = pd.NA

    # ── QC flags ───────────────────────────────────────────────────────────────
    df_star["mapping_rate_frac"] = df_star["mapping_rate"] / 100
    df_star["flag_low_reads"]   = df_star["uniquely_mapped"] < min_reads
    df_star["flag_low_rate"]    = df_star["mapping_rate_frac"] < min_rate
    df_star["flag_low_genes"]   = df_star["genes_detected"]   < min_genes
    df_star["flag_high_genes"]  = df_star["genes_detected"]   > max_genes
    df_star["qc_fail"] = (
        df_star["flag_low_reads"] |
        df_star["flag_low_rate"]  |
        df_star["flag_low_genes"] |
        df_star["flag_high_genes"]
    )

    # ── Save table ─────────────────────────────────────────────────────────────
    df_star.reset_index().to_csv(qc_out / "cell_qc_metrics.csv", index=False)
    fail_cells = df_star[df_star["qc_fail"]].index.tolist()

    print(f"\n=== QC Summary ===")
    print(f"Total cells:   {len(df_star)}")
    print(f"QC pass:       {(~df_star['qc_fail']).sum()}")
    print(f"QC fail:       {df_star['qc_fail'].sum()}")
    if fail_cells:
        print(f"Flagged cells: {', '.join(fail_cells)}")

    # ── Plots ──────────────────────────────────────────────────────────────────
    colors = ["#d9534f" if f else "#5bc0de" for f in df_star["qc_fail"]]

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    fig.suptitle("Per-cell QC metrics", fontsize=13)

    # 1. Uniquely mapped reads
    ax = axes[0]
    ax.bar(range(len(df_star)), df_star["uniquely_mapped"] / 1e6, color=colors)
    ax.axhline(min_reads / 1e6, color="grey", linestyle="--", linewidth=1)
    ax.set_xlabel("Cell")
    ax.set_ylabel("Uniquely mapped reads (M)")
    ax.set_title("Mapped reads")
    ax.set_xticks([])

    # 2. Mapping rate
    ax = axes[1]
    ax.bar(range(len(df_star)), df_star["mapping_rate"], color=colors)
    ax.axhline(min_rate * 100, color="grey", linestyle="--", linewidth=1)
    ax.set_xlabel("Cell")
    ax.set_ylabel("Mapping rate (%)")
    ax.set_title("Mapping rate")
    ax.set_xticks([])

    # 3. Genes detected
    ax = axes[2]
    if df_star["genes_detected"].notna().any():
        ax.bar(range(len(df_star)), df_star["genes_detected"], color=colors)
        ax.axhline(min_genes, color="grey", linestyle="--", linewidth=1, label=f"Min ({min_genes})")
        ax.axhline(max_genes, color="grey", linestyle="-.", linewidth=1, label=f"Max ({max_genes})")
        ax.set_ylabel("Genes detected")
        ax.set_title("Genes detected")
        ax.set_xticks([])
        ax.legend(fontsize=8)

    pass_patch = mpatches.Patch(color="#5bc0de", label="QC pass")
    fail_patch = mpatches.Patch(color="#d9534f", label="QC fail")
    fig.legend(handles=[pass_patch, fail_patch], loc="lower right", fontsize=9)

    plt.tight_layout()
    plt.savefig(qc_out / "cell_qc_barplots.pdf", bbox_inches="tight")
    plt.savefig(qc_out / "cell_qc_barplots.png", bbox_inches="tight", dpi=150)
    plt.close()

    print(f"\nQC outputs:")
    print(f"  {qc_out}/cell_qc_metrics.csv")
    print(f"  {qc_out}/cell_qc_barplots.pdf")

    # ── Write metadata file excluding QC-fail cells ────────────────────────────
    meta_file = Path(cfg["METADATA_FILE"])
    if meta_file.exists():
        meta = pd.read_csv(meta_file)
        meta_pass = meta[~meta["cell_id"].isin(fail_cells)]
        meta_pass_file = qc_out / "metadata_qc_pass.csv"
        meta_pass.to_csv(meta_pass_file, index=False)
        print(f"  {meta_pass_file}  (metadata with QC-fail cells removed)")
        print(f"\nTo run DGE on QC-passing cells only, update config.yaml:")
        print(f"  METADATA_FILE: {meta_pass_file}")


if __name__ == "__main__":
    main()
