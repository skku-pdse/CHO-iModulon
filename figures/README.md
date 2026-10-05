# Figure scripts

```bash
python figures/generate_manuscript_figures.py "figure 1" --output-dir "outputs/figures/figure 1"
```

Use `--format pdf` for vector output. The default format is PNG.

| Figure | Panels and scripts |
|---|---|
| 1 | B/C/E: `module_overview.py`; D: `explained_variance.py`; F: `metadata_distribution.py`; G: `representative_heatmap.py` |
| 2 | A/B: `lineage_violin.py`; C/D: `lineage_effects.py` + `lineage_pairwise_heatmap.py`; E: `chromosome_enrichment.py`; F/G: `genomic_panels.py` |
| 3 | A: `temperature_activity.py`; C: `tds_genes.py`; D: `antiviral_activity.py`; F: `isg1_genes.py` |
| 4 | B/C: `projected_lineage.py`; D: `projected_temperature.py`; E/F: `projected_atf4.py` |
| 5 | A–G: `make_figure5.py` and panel helpers |
| Supplementary | `within_project_lineage.py` |

Figures use the supplied publication matrices and annotations. Figures 1–4 are exported as data panels; schematic panels (1A, 3B/E, 4A) are not generated. Figure 1D uses the supplied reconstruction-EV table.
