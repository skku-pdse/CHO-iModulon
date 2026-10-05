# ICA

`TPM matrix/` contains expression inputs; `code/` contains ICA and gene-membership scripts. All commands below run from the repository root.

## TPM inputs

| File | Samples |
|---|---:|
| `TPM matrix/discovery_tpm.tsv.gz` | 778 |
| `TPM matrix/expanded_tpm.tsv.gz` | 1,030 |

These are gzip-compressed copies of the final TPM files. Rows are genes; two header rows identify project and sample. Scripts read them directly without extraction. Checksums are recorded in `tpm_input_manifest.json`.

## Run ICA

```bash
python reproduce_analysis.py --steps ica --dataset discovery
python reproduce_analysis.py --steps ica --dataset expanded
```

The workflow retains genes with TPM > 1 in any sample, applies log2(TPM + 1), and standardizes each sample across genes. Defaults are 200 Picard components, 256 seeds, and DBSCAN with eps 0.1 and min_samples 129. It saves gene weights, signed-RMS activities, cluster statistics, and logs under `outputs/<dataset>/ica/`.

## Individual stages

```bash
python ICA/code/run_robust_ica.py --dataset discovery --tpm "ICA/TPM matrix/discovery_tpm.tsv.gz" --output-dir outputs/discovery/ica --stage ica
python ICA/code/run_robust_ica.py --dataset discovery --tpm "ICA/TPM matrix/discovery_tpm.tsv.gz" --output-dir outputs/discovery/ica --stage distances
python ICA/code/run_robust_ica.py --dataset discovery --tpm "ICA/TPM matrix/discovery_tpm.tsv.gz" --output-dir outputs/discovery/ica --stage cluster
```

Completed seed matrices and similarity blocks are reused when rerunning the same configuration. A full fit requires substantial memory, disk space, and runtime. Numerical results and component numbering can vary across environments.

## Membership for a new fit

```bash
python ICA/code/determine_thresholds.py outputs/discovery/ica/discovery_ica_gene_weights.csv --output-dir outputs/discovery/ica --prefix discovery_refit_membership
```

Genes are selected at the last qualifying gap among all adjacent differences in sorted absolute weights; a gap qualifies at >= 1% of the maximum absolute weight. This produces membership and threshold tables, without assigning biological annotations.

Published annotations use the supplied matrices; newly fitted components must be matched before reusing their labels.

## Data availability

Sample-level TPM, activity, and metadata for NISTCHO_Characterisation (also named R17704 or Project_BOKU) are omitted from the general analysis inputs because the associated manuscript is in preparation. Figure-specific inputs and plotted results are retained. Supplied models and aggregate results were obtained from the full research dataset; reruns using the released subset may differ.
