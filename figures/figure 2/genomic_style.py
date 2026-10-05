import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


def natural_key(text):
    """Build a natural sort key for chromosome and scaffold identifiers."""
    text = str(text)
    return [
        int(tok) if tok.isdigit() else tok.lower() for tok in re.split("(\\d+)", text)
    ]


def add_genome_offsets(gene_pos_df, chromosome_order=None):
    """Place gene midpoints on concatenated sequences in the specified order."""
    if chromosome_order is None:
        chromosome_order = sorted(
            gene_pos_df["seqid"].dropna().unique(), key=natural_key
        )
    chrom_lengths = (
        gene_pos_df[["seqid", "chrom_length"]]
        .drop_duplicates()
        .set_index("seqid")["chrom_length"]
        .to_dict()
    )
    offset_rows = []
    current_offset = 0
    for chrom in chromosome_order:
        if chrom not in chrom_lengths:
            continue
        length = chrom_lengths[chrom]
        offset_rows.append(
            {
                "seqid": chrom,
                "offset": current_offset,
                "chrom_start": current_offset,
                "chrom_end": current_offset + length,
                "chrom_mid": current_offset + length / 2,
            }
        )
        current_offset += length
    offset_df = pd.DataFrame(offset_rows)
    offset_map = offset_df.set_index("seqid")["offset"].to_dict()
    out = gene_pos_df.copy()
    out["x"] = out["midpoint"] + out["seqid"].map(offset_map)
    return (out, offset_df)


def plot_one_module(
    df,
    offset_df,
    value_col,
    outdir,
    threshold=None,
    figsize=(14, 5),
    point_size=15,
    alpha=0.7,
    chromosome_label_map=None,
    highlight_col=None,
    highlight_label="iModulon genes",
):
    plot_df = df.dropna(subset=[value_col, "x"]).copy()
    plot_df[value_col] = plot_df[value_col].abs()
    fig, ax = plt.subplots(figsize=figsize)
    if highlight_col is not None and highlight_col in plot_df.columns:
        bg = plot_df[~plot_df[highlight_col].fillna(False)]
        hl = plot_df[plot_df[highlight_col].fillna(False)]
        ax.scatter(
            bg["x"],
            bg[value_col],
            s=point_size,
            alpha=alpha,
            color="0.6",
            edgecolors="none",
            zorder=2,
        )
        if not hl.empty:
            ax.scatter(
                hl["x"],
                hl[value_col],
                s=point_size + 15,
                alpha=0.95,
                color="tab:red",
                edgecolors="none",
                zorder=3,
                label=highlight_label,
            )
    else:
        ax.scatter(
            plot_df["x"],
            plot_df[value_col],
            s=point_size,
            alpha=alpha,
            color="0.6",
            edgecolors="none",
            zorder=2,
        )
    for _, row in offset_df.iterrows():
        ax.axvline(row["chrom_start"], color="lightgray", lw=1, zorder=0)
    if not offset_df.empty:
        ax.axvline(offset_df["chrom_end"].iloc[-1], color="lightgray", lw=1, zorder=0)
    ax.axhline(0, color="gray", lw=0.8, zorder=1)
    if threshold is not None:
        ax.axhline(threshold, color="black", ls="--", lw=1, zorder=1)
        ax.axhline(-threshold, color="black", ls="--", lw=1, zorder=1)
    xticks = offset_df["chrom_mid"].tolist()
    xlabels = offset_df["seqid"].tolist()
    if chromosome_label_map is not None:
        xlabels = [chromosome_label_map.get(x, x) for x in xlabels]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, rotation=30)
    ax.set_xlabel("Concatenated genomic position (gene midpoint; bp)")
    ax.set_ylabel("Absolute gene weight")
    ax.set_title({25: "Histone (iM25)", 8: "Ephrin Locus (iM8)"}.get(int(value_col), f"iM{value_col}"))
    plt.tight_layout()
    safe_name = str(value_col).replace("/", "_").replace("\\", "_").replace(" ", "_")
    outpath = os.path.join(outdir, f"figure2{ {25: chr(70), 8: chr(71)}.get(int(value_col), chr(88))}_chromosomal_position_iM{safe_name}.png")
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return outpath
