# Expanded motif analysis

Rebuild summaries from archived outputs (from the repository root):

```bash
python "expanded iModulon/motif analysis/summarize_expanded_motif_analysis.py" --output-dir "outputs/expanded/motifs"
```

Use `--results-dir PATH` for another numeric-module directory. MEME XML and TOMTOM TSV are sufficient; TOMTOM XML supplies TF names, otherwise target IDs are retained. AME is summarized independently. FIMO, CentriMo and SpaMo are not required.

## Rerun searches

MEME Suite 5.5.5 was used for the original analysis ([installation](https://meme-suite.org/meme/doc/install.html)). Make `meme`, `tomtom`, `ame`, and `fasta-get-markov` available on PATH. The supplied JASPAR 2024 CORE non-redundant databases (879 vertebrate and 224 mouse motifs) are selected automatically:

```bash
python "expanded iModulon/motif analysis/run_expanded_motif_analysis.py" --output-dir "outputs/expanded/motif_rerun" --jobs 1 --threads 8
```

Use `--vertebrate-db` and `--mouse-db` to override the supplied databases. Use a new output directory. `--modules 0 1` selects modules; `--dry-run` prepares inputs and commands only. Concurrent CPU demand is approximately jobs times threads.

The archive supplies original promoter FASTAs and BED intervals from CriGri-PICRH-1.0 gene boundaries (500 bp upstream, 100 bp downstream). Controls contain other member genes in the cohort union. MEME uses ZOOPS, 10 motifs, widths 6-40, both strands, and a first-order background. TOMTOM uses Pearson distance separately for each database. AME uses the vertebrate database, average scoring, Fisher's method and the supplied controls.

## Results

- `*_meme_tomtom_summary.csv`: publication-style columns; MEME E-value < 0.01. Separate mouse and vertebrate CSVs are also generated. Each module/motif/database/TF retains its lowest TOMTOM q-value; cross-database hits remain separate.
- `*_tomtom_detailed_hits.csv`: target motif IDs retained, duplicate target hits collapsed by lowest q-value. HIGH confidence means the row's own TOMTOM q-value < 0.05.
- `*_meme_motif_inventory.csv`: all MEME motifs, including motifs without TOMTOM hits.
- `*_ame_summary.csv`: AME statistics and TP/FP percentages. AME_Significance preserves the original adjusted-p-value rule; AME_Evalue_significant identifies E-value < 0.05.
- `*_motif_coverage.csv`: missing outputs distinguished from completed searches with no hits.

Total_Sites uses optional archived FIMO counts (unique sequence IDs with q-value < 0.01). This field stays blank for new MEME/TOMTOM-only searches; MEME site counts are not substituted. Database-specific processing avoids the original multiplication of rows when joining both databases. TOMTOM similarity does not establish direct TF binding.

Expanded iModulon 106 contains one gene and is excluded from motif analysis. Its coverage status is `excluded_single_gene`; it is also skipped during reruns.
