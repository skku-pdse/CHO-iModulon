import numpy as np

BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
PURPLE = "#6F4AA8"


def clean(ax):
    ax.spines[["top", "right"]].set_visible(False)


def letter(ax, s):
    ax.text(-0.16, 1.08, s, transform=ax.transAxes, fontweight="bold", fontsize=12)


def panelmatch(ax, d, l, DM, EM, dg, eg, t):
    """Plot one matched pair using its shared and cohort-specific members."""
    z = t[t.discovery_iModulon.eq(d)].iloc[0]
    e = int(z.expanded_iModulon)
    ds = set(dg[dg.iModulon.eq(d)].Gene)
    es = set(eg[eg.iModulon.eq(e)].Gene)
    genes = DM.index
    xx = DM[d].to_numpy()
    yy = EM[e].to_numpy() * np.sign(z.weight_r)
    xx = xx / np.std(xx)
    yy = yy / np.std(yy)
    masks = [
        ~genes.isin(ds | es),
        genes.isin(ds - es),
        genes.isin(es - ds),
        genes.isin(ds & es),
    ]
    for mask, c, sz in zip(masks, ["#D9D9D9", ORANGE, GREEN, BLUE], [2, 11, 11, 13]):
        ax.scatter(
            xx[mask],
            yy[mask],
            s=sz,
            c=c,
            alpha=0.65,
            edgecolors="none",
            rasterized=True,
        )
    ax.axhline(0, c="#ddd", lw=0.5)
    ax.axvline(0, c="#ddd", lw=0.5)
    ax.set(
        xlabel="Discovery gene weight / SD",
        ylabel="Expanded gene weight / SD\n(sign-aligned)",
        title=f"{z.discovery_name}\nd-iM{d} → e-iM{e}",
    )
    ax.text(
        0.03,
        0.97,
        f"|r| = {z.weight_abs_r:.3f}\nMembership J = {z.membership_J:.3f}\nShared = {int(z.shared_genes)} ({int(z.discovery_n)} → {int(z.expanded_n)})",
        va="top",
        transform=ax.transAxes,
        fontsize=6.8,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85},
    )
    clean(ax)
    letter(ax, l)
    if d == 39:
        hist = ["LOC100753502", "LOC100754104", "LOC107979446", "LOC100753804"]
        mask = genes.isin(hist)
        assert set(hist) <= ds - es
        ax.scatter(
            xx[mask],
            yy[mask],
            marker="D",
            s=27,
            facecolors="none",
            edgecolors=PURPLE,
            lw=1,
        )
