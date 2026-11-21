# DepMap data: download instructions and provenance

This directory documents the DepMap data used for the real-data experiments
and explains how to obtain the complete gene-effect matrix for analysing
other genes. The data source is **DepMap 24Q4 Public, version 1**.

## Source and citation

Broad DepMap (2024). *DepMap 24Q4 Public*. Figshare+. Dataset.
[DOI: 10.25452/figshare.plus.27993248.v1](https://doi.org/10.25452/figshare.plus.27993248.v1).

- [Release page](https://plus.figshare.com/articles/dataset/DepMap_24Q4_Public/27993248)
- [Version 1 file metadata](https://api.figshare.com/v2/articles/27993248/versions/1)
- [DepMap portal](https://depmap.org/portal/)

Use this release to reproduce the supplied experiments. Later releases may
contain different scores, models or annotations.

## Required files

| File | Contents | Download size |
|---|---|---:|
| `CRISPRGeneEffect.csv` | Chronos gene-effect scores: rows identify cell-line models and columns identify genes. | 428,678,699 bytes |
| `Model.csv` | Cell-line metadata, including `ModelID`, `PatientID`, recorded `Sex` and `OncotreeLineage`. | 645,696 bytes |

The first column of the gene-effect file contains model identifiers; some
CSV readers may give an unnamed identifier column a generated name. Gene
columns use names such as `RPL38 (6169)` and `RPS23 (6228)`. Join this
identifier to `ModelID` in the metadata.

Each score describes the effect of knocking out a gene in a cancer cell-line
model. Negative scores indicate reduced cell growth. Rows represent models,
not independent clinical treatment responses; more than one model may
share a `PatientID`.

## Download

Download these two individual files, rather than the complete release:

- [CRISPRGeneEffect.csv](https://ndownloader.figshare.com/files/51064667)
- [Model.csv](https://ndownloader.figshare.com/files/51065297)

Save both in this root-level `data/` directory. Alternatively, run the
following commands **from the repository root**:

```bash
mkdir -p data

curl --fail --location --retry 3 \
  https://ndownloader.figshare.com/files/51064667 \
  --output data/CRISPRGeneEffect.csv

curl --fail --location --retry 3 \
  https://ndownloader.figshare.com/files/51065297 \
  --output data/Model.csv
```

These commands replace existing files with the same names. Keep the full
downloaded matrix separate from any single-gene extracts.

The release records the following MD5 checksums:

| File | MD5 |
|---|---|
| `CRISPRGeneEffect.csv` | `6edf7ade09b9b34199210b559d4745d3` |
| `Model.csv` | `675210d17675f3517b0ce39a3c274f16` |

On macOS, check the downloaded files with:

```bash
md5 data/CRISPRGeneEffect.csv data/Model.csv
```

On Linux, use `md5sum` with the same paths. These checksums refer to the
original full files, not to gene subsets or files rewritten by CSV software.

## Data already supplied with the experiments

The repository includes small, single-gene extracts in:

- [`experiments/depmap/RPL38/data/`](../experiments/depmap/RPL38/data/)
- [`experiments/depmap/RPS23/data/`](../experiments/depmap/RPS23/data/)

Their files named `CRISPRGeneEffect.csv` contain only the selected gene and
model identifiers, retaining the 1,178 model rows from the original matrix.
They are sufficient for their corresponding experiments, but cannot be
used to analyse a gene that is absent from the extract.

The full `Model.csv` is also included in each experiment directory.
Release references and extraction details are recorded in the accompanying
provenance files. Cohort filtering is applied by the analysis code, not
pre-applied to these input CSVs.

Small extracts may also be stored in this root-level directory using
explicit names, for example `CRISPRGeneEffect_RPL38.csv` and
`CRISPRGeneEffect_RPS23.csv`. Reserve `data/CRISPRGeneEffect.csv` for the
full downloaded matrix. Document any newly added extract's gene columns,
row selection and source release.

## Using the root-level data

The existing experiment scripts default to their own `data/` directories.
To use the full root-level downloads instead, supply explicit paths.
For example, from the repository root:

```bash
python -u experiments/depmap/RPL38/run_analysis.py \
  --gene-effect data/CRISPRGeneEffect.csv \
  --models data/Model.csv \
  --gene RPL38 \
  --sex Female \
  --min-source-n 10 \
  --min-target-n 10 \
  --permutations 100 \
  --alpha 0.05 \
  --example-target Myeloid \
  --out experiments/depmap/RPL38/results_reproduced
```

Install the experiment's dependencies first, as described in its README.
Use a fresh output directory, or use the runner's `--overwrite` option
when intentionally replacing an existing completed run.

For another gene in the full matrix, change `--gene` to its symbol. Also
choose the appropriate recorded sex, an eligible `--example-target` lineage
and a separate output directory. Check that figure titles and descriptions
match the new experiment. The loader reads the identifier and selected
gene column; a separate full-matrix download for each gene is unnecessary.

The analysis filters unusable scores and metadata, applies the selected
sex filter, retains one model per recorded patient and enforces the stated
lineage sample-size thresholds. Target lineages are held out from prior
fitting. Eligibility and the Gaussian working assumptions must still be
assessed for each analysis; availability of a gene does not establish that
its data satisfy those assumptions.

The example specifies a 5% significance level explicitly. Always record
the chosen level and other settings rather than relying on runner defaults.

## Repository storage

The full gene-effect matrix is approximately 429 MB and exceeds GitHub's
[100 MiB limit for ordinary Git files](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github).
Download it locally using the instructions above. Commit this README,
metadata and small clearly named extracts as needed.

Add this line to the repository's root `.gitignore` to exclude only the
full root-level matrix:

```gitignore
/data/CRISPRGeneEffect.csv
```

## Licence and attribution

The DepMap 24Q4 Public release is distributed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Credit Broad DepMap, cite the version-specific dataset DOI and identify any
extraction or preprocessing when redistributing derived files. The
repository's software licence does not replace the data's original licence.
