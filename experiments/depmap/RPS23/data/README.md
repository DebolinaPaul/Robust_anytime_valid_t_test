# Data provenance

Source: Broad DepMap, **DepMap 24Q4 Public**, version 1.

- Release: https://plus.figshare.com/articles/dataset/DepMap_24Q4_Public/27993248
- DOI: https://doi.org/10.25452/figshare.plus.27993248.v1
- Version metadata: https://api.figshare.com/v2/articles/27993248/versions/1
- Full `CRISPRGeneEffect.csv`: https://ndownloader.figshare.com/files/51064667
- `Model.csv`: https://ndownloader.figshare.com/files/51065297

The full gene-effect matrix was downloaded from the release. Its size is
428,678,699 bytes and its MD5 is `6edf7ade09b9b34199210b559d4745d3`, matching
the release metadata. This package contains an extract with the model
identifier and **RPS23 (6228)** score column, retaining all 1,178 model rows.
The extract was written using pandas; the original release contains many
additional gene columns.

`Model.csv` is the 2,105-row metadata file included in the supplied RPL38
package. Its MD5 is `675210d17675f3517b0ce39a3c274f16`, matching the release.
The model identifier joins the scores to PatientID, recorded Sex and
OncotreeLineage. The CSV is retained in full so that the same filtering code
can be used.

No source/target filtering is pre-applied to these CSVs. `cohort.py` applies
the sex, metadata and donor filters; the runner then excludes lineages with
fewer than ten observations from fitting and testing. Input SHA256 hashes
are recorded in `results/settings.json`.
