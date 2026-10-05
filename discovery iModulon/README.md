# Discovery iModulons

105 iModulons identified from 778 samples.

| Folder | Contents |
|---|---|
| `data/` | Gene weights, activities, membership, names, and sample metadata |
| `functional enrichment/` | STRING functional enrichment code and results |
| `genomic/` | Chromosomal enrichment code and results |
| `conditional association/` | Associations between activities and categorical sample metadata |
| [motif analysis/](motif%20analysis/discovery_motif_README.md) | Promoter inputs, MEME/TOMTOM/AME code, and summaries |

## Data files

Files in `data/` use the `discovery_` prefix.

| Filename suffix | Contents |
|---|---|
| `gene_weight_matrix.csv` | Genes in rows; iModulons in columns |
| `activity_matrix.csv` | iModulons in rows; samples in columns, with Project/SampleID headers |
| `gene_membership.csv` | Member genes, loadings, ranks, and genomic locations |
| `imodulon_names.csv` | Module names and annotations |
| `sample_metadata.csv` | Sample identifiers, cell lines, and experimental conditions |
| `string_protein_mapping.csv` | Gene-to-STRING protein identifiers |

Match samples using both `Project` and `SampleID`.

From the repository root, run selected analyses with:

```bash
python reproduce_analysis.py --dataset discovery --steps membership genomic association motifs
```
