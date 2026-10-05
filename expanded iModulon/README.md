# Expanded iModulons

148 iModulons identified from 1,155 samples; sample-level inputs are released for 1,030 samples.

| Folder | Contents |
|---|---|
| `data/` | Gene weights, activities, membership, names, and sample metadata |
| `functional enrichment/` | STRING functional enrichment code and results |
| `genomic/` | Chromosomal enrichment code and results |
| `conditional association/` | Associations between activities and categorical sample metadata |
| [motif analysis/](motif%20analysis/expanded_motif_README.md) | Promoter inputs, MEME/TOMTOM/AME code, and summaries |
| `iModulon comparison/` | Discovery–expanded module matching, activity comparisons, and result tables |

## Data files

Files in `data/` use the `expanded_` prefix.

| Filename suffix | Contents |
|---|---|
| `gene_weight_matrix.csv` | Genes in rows; iModulons in columns |
| `activity_matrix.csv` | iModulons in rows; samples in columns, with Project/SampleID headers |
| `gene_membership.csv` | Member genes, loadings, ranks, and genomic locations |
| `imodulon_names.csv` | Module names and annotations |
| `sample_metadata.csv` | Sample identifiers, cell lines, and experimental conditions |
| `string_protein_mapping.csv` | Gene-to-STRING protein identifiers |

Match samples using both `Project` and `SampleID`.

`expanded_association_media_labels.csv` retains the media-group labels used for the reported association analysis.

From the repository root, run selected analyses with:

```bash
python reproduce_analysis.py --dataset expanded --steps membership genomic association motifs
```

## Data availability

Sample-level TPM, activity, and metadata for NISTCHO_Characterisation (also named R17704 or Project_BOKU) are omitted from the general analysis inputs because the associated manuscript is in preparation. Figure-specific inputs and plotted results are retained. Supplied models and aggregate results were obtained from the full research dataset; reruns using the released subset may differ.
