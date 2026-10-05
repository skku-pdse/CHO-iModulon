# CHO-iModulon

Data and code for ICA-based analysis of CHO cell transcriptomes: 105 discovery iModulons from 778 samples and 148 expanded-compendium iModulons from 1,155 samples.

## Data availability

Sample-level TPM, activity, and metadata for NISTCHO_Characterisation (also named R17704 or Project_BOKU) are omitted from the general analysis inputs because the associated manuscript is in preparation. Figure-specific inputs and plotted results are retained. Supplied models and aggregate results were obtained from the full research dataset; reruns using the released subset may differ.

## Install

Tested with Python 3.12. From a terminal:

```bash
git clone https://github.com/skku-pdse/CHO-iModulon.git
cd CHO-iModulon
python -m venv .venv
```

Activate the environment with `.venv\Scripts\activate` on Windows or `source .venv/bin/activate` on macOS/Linux, then install:

```bash
python -m pip install -r requirements.txt
```

## Reproduce the analyses and figures

```bash
python reproduce_analysis.py
```

This checks the inputs, recalculates gene membership, chromosome enrichment, metadata associations, external-sample projection, compendium comparisons, and MEME/TOMTOM/AME summaries, then generates Figures 1–5 data panels and the supplementary plot. It uses the supplied iModulon matrices and saves results and execution logs under `outputs/`, preserving the packaged data.

The run also writes `reference_comparison.csv`. Expanded association uses the original media-group labels retained alongside the curated metadata.

To run selected steps:

```bash
python reproduce_analysis.py --steps association --dataset discovery
python reproduce_analysis.py --steps figures
```

STRING results are supplied for offline use. To query STRING again (requires internet):

```bash
python reproduce_analysis.py --steps string
```

## Contents

| Directory | Contents |
|---|---|
| [ICA/](ICA/README.md) | Final discovery/expanded TPM inputs, seeds, and ICA code |
| [discovery iModulon/](discovery%20iModulon/README.md) | Discovery M/A matrices, metadata, annotated memberships, and statistical results |
| [expanded iModulon/](expanded%20iModulon/README.md) | Expanded matrices, annotations, results, and compendium comparison |
| [projection/](projection/README.md) | External TPM, sample selection, projected activities, and projection code |
| `figures/` | Figure inputs, scripts, and PNG panels; see [panel mapping](figures/README.md) |
| `scripts/` | Shared motif analysis code |
| `validation/` | Input checks and comparison with reference results |
| `databases/` | JASPAR 2024 vertebrate and mouse motif databases |
| [tests/](tests/README.md) | Small numerical and execution tests |

Gene-weight matrices have genes in rows and iModulons in columns; activity matrices have iModulons in rows and samples in columns. Samples are matched by the combined `Project` and `SampleID` identifiers. Membership tables link `Gene` identifiers to `iModulon` numbers.

The workflow starts from quantified TPM. Raw-read processing and manual schematic artwork are outside this package.
