# Projection

Project external expression onto the fixed discovery gene-weight matrix.

| File in `data/` | Contents |
|---|---|
| `validation_tpm.tsv.gz` | TPM input with 258 samples; project and sample header rows |
| `projected_sample_metadata.csv` | Metadata selecting the 255 projected samples |
| `projected_activity_matrix.csv` | 105 iModulons × 255 samples |
| `validation_tpm_manifest.json` | Input checksum and sample count |

From the repository root:

```bash
python reproduce_analysis.py --steps projection
```

`project.py` selects samples by Project/SampleID, aligns genes to the discovery basis, log-transforms and standardizes TPM within each sample, and solves the least-squares projection. Results are saved under `outputs/projection/`.

## Data availability

Sample-level TPM, activity, and metadata for NISTCHO_Characterisation (also named R17704 or Project_BOKU) are omitted from the general analysis inputs because the associated manuscript is in preparation. Figure-specific inputs and plotted results are retained. Supplied models and aggregate results were obtained from the full research dataset; reruns using the released subset may differ.
