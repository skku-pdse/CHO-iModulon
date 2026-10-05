"""
Figure 2E: Volcano plot of (iModulon, chromosome) enrichments.

Purpose
-------
Demonstrates that iM11 (shown in panel C as a Manhattan plot) is not an
isolated cherry-picked example. Across all 105 iModulons x 12 chromosomes,
several (iM, chromosome) pairs show statistically significant enrichment,
identifying a category of "locus iModulons" whose gene composition is
spatially clustered on the genome rather than diffuse.

Statistical model
-----------------
For each (iModulon i, chromosome c) pair:
  N = total unique genes assigned to any iModulon (background pool)
  K_c = unique genes on chromosome c (in background)
  n_i = unique genes in iModulon i
  k_{i,c} = unique genes in iModulon i AND on chromosome c

  expected_{i,c} = n_i * K_c / N
  fold_{i,c} = k_{i,c} / expected_{i,c}
  p_{i,c} = P(X >= k_{i,c}) under hypergeom(N, K_c, n_i)   [one-sided, greater]

Multiple-testing correction: Benjamini-Hochberg across all (i, c) pairs tested.

Filtering
---------
- iModulons with n_i < MIN_GENES_PER_IM excluded (too few genes for enrichment).
- (i, c) pairs where k_{i,c} == 0 excluded from the volcano (depletion not of
  interest here; mentioned in caption).

Inputs
------
- iModulon_analysis_final_claude.xlsx, sheet 'iModulon_gene'
  Required columns: 'iModulon', 'Gene', 'Chromosome'
- chrom.sizes (only used to define chromosome ORDER and label mapping;
  NW_023276806.1 + NW_023276807.1 already merged to Chromosome 1 in xlsx,
  scaffolds collapsed to 'Z' in xlsx).

Output
------
- figure2E_chromosome_enrichment.png / .pdf
- chr_enrichment_results.csv (full table for supplement)
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import hypergeom
try:
    from adjustText import adjust_text
    HAS_ADJUSTTEXT = True
except ImportError:
    HAS_ADJUSTTEXT = False
    print("[WARN] adjustText not installed -> label positions may overlap. "
          "pip install adjustText to enable auto-placement.")


# =========================
# 1. Paths
# =========================
from figure2_paths import DATA, RESULTS
PANEL_DIR = DATA
IMODULON_FILE = DATA / "discovery.xlsx"
OUTDIR = RESULTS / "lineage"
OUTDIR.mkdir(parents=True, exist_ok=True)

# =========================
# 2. Configuration
# =========================
# Chromosome order for plotting / coloring (autosomes -> X -> unplaced scaffold bucket)
# Internal label "Z" in source xlsx is renamed to "U" (unplaced scaffold) for display.
CHR_ORDER = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 'X', 'Z']
CHR_DISPLAY = {'Z': 'U'}  # internal -> display label
def disp(c):
    return CHR_DISPLAY.get(c, c)

# Filtering
MIN_GENES_PER_IM = 5     # Skip iModulons with < this many genes
MIN_OBS = 1              # Plot only (iM, chr) pairs with at least this many obs

# Significance thresholds (also used as guide lines)
FDR_CUTOFF = 0.05
LOG2FC_CUTOFF = 1.0      # i.e. fold-enrichment >= 2

# Top-hit labeling
N_LABELS = 12            # Label this many top significant points
ALWAYS_LABEL = [(11, '2')]  # Always label these (i, chr) — matches Panel C example

# Plot tuning
FIGSIZE = (9, 5.6)
DOT_SIZE_NS = 15
DOT_SIZE_SIG = 50
DOT_ALPHA = 0.85
LABEL_FONTSIZE = 18
AXIS_FONTSIZE = 18
TITLE_FONTSIZE = 11
TICK_FONTSIZE = 18

# Broken-axis settings
# Lower panel = "main" range, upper panel = outlier range.
# Auto-detected at runtime, but with these as fallback bounds.
BREAK_GAP_FRAC = 0.04    # gap between two panels (fraction of figure height)
LOWER_HEIGHT_FRAC = 0.78 # bottom panel takes this much of plot area
OUTLIER_THRESHOLD = 20   # -log10(p_adj) above this -> "outlier" panel


# =========================
# 3. Load data
# =========================
def load_imodulon_gene(path):
    df = pd.read_csv(DATA.parents[1] / "discovery iModulon/data/discovery_gene_membership.csv")
    df = df[['iModulon', 'Gene', 'Chromosome']].copy()
    df['Chromosome'] = df['Chromosome'].astype(str).replace({'0': 'Z', 'U': 'Z'})
    df['iModulon'] = df['iModulon'].astype(int)
    # Drop rows with missing values
    df = df.dropna(subset=['Gene', 'Chromosome'])
    df = df[df['Chromosome'].isin(CHR_ORDER)].copy()
    return df


def load_imodulon_names(path):
    try:
        s = pd.read_excel(path, sheet_name='iModulon Summary')
        s['iModulon'] = s['iModulon'].astype(int)
        return dict(zip(s['iModulon'], s['Name']))
    except Exception as e:
        print(f"[WARN] Could not load iModulon names: {e}")
        return {}


# =========================
# 4. Compute enrichment
# =========================
def compute_chr_enrichment(df, min_genes_per_im=MIN_GENES_PER_IM):
    """
    For each (iModulon, chromosome) pair compute hypergeometric one-sided p-value
    for enrichment (greater).

    Background N = unique genes across all iModulons.
    """
    # Unique-gene-level (gene appears once per iM in the table; collapse just in case)
    df_u = df.drop_duplicates(subset=['iModulon', 'Gene']).copy()

    # Background pool: unique genes ever assigned to any iModulon
    bg = df_u.drop_duplicates(subset=['Gene'])
    N = len(bg)
    K_per_chr = bg['Chromosome'].value_counts().to_dict()
    print(f"[INFO] Background N = {N} unique genes")
    print(f"[INFO] Genes per chromosome: { {c: K_per_chr.get(c, 0) for c in CHR_ORDER} }")

    rows = []
    for im, sub in df_u.groupby('iModulon'):
        n_i = sub['Gene'].nunique()
        if n_i < min_genes_per_im:
            continue
        for c in CHR_ORDER:
            K_c = K_per_chr.get(c, 0)
            if K_c == 0:
                continue
            k_ic = (sub['Chromosome'] == c).sum()
            expected = n_i * K_c / N
            fold = (k_ic / expected) if expected > 0 else np.nan
            # One-sided hypergeometric (greater): P(X >= k)
            p_greater = hypergeom.sf(k_ic - 1, N, K_c, n_i) if k_ic > 0 else 1.0
            rows.append({
                'iModulon': im,
                'Chromosome': c,
                'n_iModulon': n_i,
                'K_chromosome': K_c,
                'N_background': N,
                'observed': int(k_ic),
                'expected': expected,
                'log2_fold': np.log2(fold) if fold and fold > 0 else np.nan,
                'p_value': p_greater,
            })
    res = pd.DataFrame(rows)

    # BH-FDR across all tested (iM, chr) pairs
    res = res.sort_values('p_value').reset_index(drop=True)
    m = len(res)
    res['rank'] = np.arange(1, m + 1)
    res['p_adj_BH'] = (res['p_value'] * m / res['rank']).clip(upper=1.0)
    # Enforce monotonic non-decreasing p_adj from smallest p upward
    res['p_adj_BH'] = res['p_adj_BH'][::-1].cummin()[::-1]
    res = res.drop(columns=['rank']).sort_values(['iModulon', 'Chromosome']).reset_index(drop=True)
    return res


# =========================
# 5. Plot volcano (broken y-axis + adjustText)
# =========================
def _draw_break_marks(ax_top, ax_bot, d=0.012):
    gap = 0.03
    # """Add diagonal cut marks at the break between two stacked axes."""
    # kw = dict(transform=ax_top.transAxes, color='k', clip_on=False, lw=1.5)
    # ax_top.plot((-d, +d), (-d, +d), **kw)
    # ax_top.plot((-d, +d),(-0.08 - d, -0.08 + d),  **kw)
    #
    # kw = dict(transform=ax_bot.transAxes, color='k', clip_on=False, lw=1.5)
    # ax_bot.plot((1 -d, 1 +d), (- d, + d), **kw)
    # ax_bot.plot((1 -d, 1 +d),(-0.08 - d, -0.08 + d),  **kw)


def plot_volcano(res, name_map, out_png, out_pdf,
                 fdr=FDR_CUTOFF, log2fc=LOG2FC_CUTOFF, n_labels=N_LABELS,
                 outlier_threshold=OUTLIER_THRESHOLD):
    # Filter to plotted points
    plot_df = res[(res['observed'] >= MIN_OBS) & res['log2_fold'].notna()].copy()
    plot_df['neg_log10_padj'] = -np.log10(plot_df['p_adj_BH'].clip(lower=1e-300))
    plot_df['significant'] = (plot_df['p_adj_BH'] < fdr) & (plot_df['log2_fold'] >= log2fc)

    # Color palette
    cmap = plt.get_cmap('tab20')
    chr_color = {c: cmap(i / max(len(CHR_ORDER) - 1, 1)) for i, c in enumerate(CHR_ORDER)}

    # Detect outliers needing top panel
    outlier_mask = plot_df['neg_log10_padj'] > outlier_threshold
    has_outliers = outlier_mask.any()

    # Determine y-ranges
    y_main_max = max(outlier_threshold * 0.6,
                     plot_df.loc[~outlier_mask, 'neg_log10_padj'].max() * 1.10)
    if has_outliers:
        y_top_min = plot_df.loc[outlier_mask, 'neg_log10_padj'].min() * 0.95
        y_top_max = plot_df.loc[outlier_mask, 'neg_log10_padj'].max() * 1.05
    else:
        y_top_min = y_top_max = None

    # Figure / axes layout
    fig = plt.figure(figsize=FIGSIZE)
    if has_outliers:
        # Two stacked panels with shared x
        gs = fig.add_gridspec(2, 1,
                              height_ratios=[1 - LOWER_HEIGHT_FRAC, LOWER_HEIGHT_FRAC],
                              hspace=BREAK_GAP_FRAC * 2)
        ax_top = fig.add_subplot(gs[0, 0])
        ax_bot = fig.add_subplot(gs[1, 0], sharex=ax_top)
        axes = [ax_top, ax_bot]
    else:
        ax_bot = fig.add_subplot(1, 1, 1)
        ax_top = None
        axes = [ax_bot]

    # Plot points on both axes (matplotlib will clip outside ylim)
    def _scatter_all(ax):
        # Non-significant: light gray
        ns = plot_df[~plot_df['significant']]
        ax.scatter(ns['log2_fold'], ns['neg_log10_padj'],
                   s=DOT_SIZE_NS, c='lightgray', alpha=0.5, edgecolors='none',
                   zorder=1)
        # Significant: colored by chromosome
        for c in CHR_ORDER:
            sig_c = plot_df[plot_df['significant'] & (plot_df['Chromosome'] == c)]
            if len(sig_c):
                ax.scatter(sig_c['log2_fold'], sig_c['neg_log10_padj'],
                           s=DOT_SIZE_SIG, c=[chr_color[c]],
                           alpha=DOT_ALPHA, edgecolors='black', linewidths=0.5,
                           zorder=3)

    for ax in axes:
        _scatter_all(ax)

    # Threshold guide lines (only on bottom panel)
    ax_bot.axhline(-np.log10(fdr), color='gray', ls='--', lw=1.5, alpha=0.7,
                   zorder=2)
    ax_bot.axvline(log2fc, color='gray', ls='--', lw=1.5, alpha=0.7, zorder=2)
    if has_outliers:
        ax_top.axvline(log2fc, color='gray', ls='--', lw=1.5, alpha=0.7,
                       zorder=2)

    # Y-limits
    ax_bot.set_ylim(0, y_main_max)
    if has_outliers:
        ax_top.set_ylim(y_top_min, y_top_max)
        # Hide spines between
        ax_bot.spines['top'].set_visible(False)
        ax_top.spines['bottom'].set_visible(False)
        ax_top.tick_params(labelbottom=False, bottom=False)
        _draw_break_marks(ax_top, ax_bot)
    else:
        ax_bot.spines['top'].set_visible(False)
    # ax_bot.spines['right'].set_visible(False)
    # if has_outliers:
    #     ax_top.spines['right'].set_visible(False)
    #
    # # ---- Labeling with adjustText ----
    # # Pick label set: top N significant (by p_adj) + ALWAYS_LABEL
    # sig = plot_df[plot_df['significant']].copy()
    # sig = sig.sort_values('p_adj_BH').head(n_labels)
    # label_points = sig.copy()
    # for im, ch in ALWAYS_LABEL:
    #     forced = plot_df[(plot_df['iModulon'] == im) & (plot_df['Chromosome'] == ch)]
    #     label_points = pd.concat([label_points, forced]).drop_duplicates(
    #         subset=['iModulon', 'Chromosome'])
    #
    # # Build labels and assign each to the correct axis based on y-value
    # texts_top, texts_bot = [], []
    # for _, row in label_points.iterrows():
    #     nm = name_map.get(int(row['iModulon']), '')
    #     nm = (nm[:18] + '…') if len(nm) > 19 else nm
    #     text = f"iM{int(row['iModulon'])}{(' ' + nm) if nm else ''} @Chr{disp(row['Chromosome'])}"
    #     x, y = row['log2_fold'], row['neg_log10_padj']
    #     if has_outliers and y > outlier_threshold:
    #         t = ax_top.text(x, y, text,
    #                         fontsize=LABEL_FONTSIZE, ha='left', va='bottom',
    #                         bbox=dict(boxstyle='round,pad=0.2', fc='white',
    #                                   ec='gray', lw=0.4, alpha=0.85))
    #         texts_top.append(t)
    #     else:
    #         t = ax_bot.text(x, y, text,
    #                         fontsize=LABEL_FONTSIZE, ha='left', va='bottom',
    #                         bbox=dict(boxstyle='round,pad=0.2', fc='white',
    #                                   ec='gray', lw=0.4, alpha=0.85))
    #         texts_bot.append(t)
    #
    # if HAS_ADJUSTTEXT:
    #     if texts_bot:
    #         adjust_text(
    #             texts_bot, ax=ax_bot,
    #             expand_points=(1.4, 1.6), expand_text=(1.2, 1.4),
    #             arrowprops=dict(arrowstyle='-', color='gray',
    #                             lw=0.5, alpha=0.7),
    #             only_move={'points': 'xy', 'text': 'xy'},
    #         )
    #     if has_outliers and texts_top:
    #         adjust_text(
    #             texts_top, ax=ax_top,
    #             expand_points=(1.4, 1.6), expand_text=(1.2, 1.4),
    #             arrowprops=dict(arrowstyle='-', color='gray',
    #                             lw=0.5, alpha=0.7),
    #             only_move={'points': 'xy', 'text': 'xy'},
    #         )
    #
    # ---- Axis labels / title / legend ----
    ax_bot.set_xlabel(r'log$_2$(observed / expected)',
                      fontsize=AXIS_FONTSIZE)
    # Single shared y-label centered between panels
    fig.text(0.02, 0.5, r'$-$log$_{10}$(BH-adjusted $p$)',
             rotation=90, va='center', ha='center', fontsize=AXIS_FONTSIZE)
    # fig.suptitle(
    #     f'Chromosomal enrichment of iModulon gene composition\n'
    #     f'(hypergeometric one-sided; dashed: FDR={fdr}, fold≥{2**log2fc:.0f}; '
    #     f'U = unplaced scaffold bucket)',
    #     fontsize=TITLE_FONTSIZE, y=0.995)

    # # Legend on bottom panel only
    # legend_handles = []
    # legend_labels = []
    # for c in CHR_ORDER:
    #     n_sig_c = int(plot_df[plot_df['significant'] & (plot_df['Chromosome'] == c)].shape[0])
    #     if n_sig_c > 0:
    #         legend_handles.append(
    #             plt.Line2D([0], [0], marker='o', color='none',
    #                        markerfacecolor=chr_color[c], markeredgecolor='black',
    #                        markeredgewidth=0.5, markersize=7))
    #         legend_labels.append(f'Chr {disp(c)} (n={n_sig_c})')
    # legend_handles.insert(0, plt.Line2D([0], [0], marker='o', color='none',
    #                                      markerfacecolor='lightgray',
    #                                      markeredgecolor='none', markersize=5))
    # legend_labels.insert(0, 'ns')
    # # ax_bot.legend(legend_handles, legend_labels,
    # #               loc='upper left', bbox_to_anchor=(1.02, 1.0),
    # #               fontsize=8, frameon=False, title='Chromosome', title_fontsize=9)

    for ax in axes:
        ax.grid(True, ls=':', lw=0.4, alpha=0.5)
        ax.tick_params(axis='both', labelsize=TICK_FONTSIZE)
    #
    # # Counts annotation
    # n_sig = int(plot_df['significant'].sum())
    # n_im_sig = plot_df.loc[plot_df['significant'], 'iModulon'].nunique()
    # ax_bot.text(0.98, 0.02,
    #             f"{n_sig} sig. (iM, chr) pairs\nfrom {n_im_sig} iModulons",
    #             transform=ax_bot.transAxes, ha='right', va='bottom',
    #             fontsize=9, bbox=dict(boxstyle='round,pad=0.3',
    #                                    fc='white', ec='gray', lw=0.4, alpha=0.85))

    plt.subplots_adjust(left=0.12, right=0.78, top=0.90, bottom=0.10)
    fig.savefig(out_png, dpi=400, bbox_inches='tight')
    fig.savefig(out_pdf, bbox_inches='tight')
    print(f"[SAVED] {out_png}")
    print(f"[SAVED] {out_pdf}")
    return fig


# =========================
# 6. Main
# =========================
if __name__ == '__main__':
    print('[INFO] Loading inputs ...')
    df = load_imodulon_gene(IMODULON_FILE)
    print(f'[INFO] iModulon_gene: {len(df)} rows, '
          f'{df["iModulon"].nunique()} iModulons, '
          f'{df["Gene"].nunique()} unique genes')

    name_map = load_imodulon_names(IMODULON_FILE)

    print('[INFO] Computing chromosomal enrichments ...')
    res = compute_chr_enrichment(df)

    # Save full table for supplement
    out_csv = OUTDIR / 'chr_enrichment_results.csv'
    res.to_csv(out_csv, index=False)
    print(f'[SAVED] {out_csv}')

    # Quick summary
    sig_mask = (res['p_adj_BH'] < FDR_CUTOFF) & (res['log2_fold'] >= LOG2FC_CUTOFF)
    print(f'[INFO] Significant (FDR<{FDR_CUTOFF}, fold>=2): {sig_mask.sum()} '
          f'(iM, chr) pairs from {res.loc[sig_mask, "iModulon"].nunique()} iModulons')
    print('[INFO] Top 10 enrichments:')
    print(res.sort_values('p_adj_BH').head(10).to_string(index=False))

    print('[INFO] Building volcano plot ...')
    out_png = OUTDIR / 'figure2E_chromosome_enrichment.png'
    out_pdf = OUTDIR / 'figure2E_chromosome_enrichment.pdf'
    plot_volcano(res, name_map, out_png, out_pdf)
