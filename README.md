# Resistance-associated substitution profiles outperform whole-genome sequences for machine learning prediction of hepatitis C treatment failure

Data and analysis code for the manuscript of the same name.

162 hepatitis C virus whole-genome sequences, retrieved from NCBI GenBank and
contributed by 146 patients, were used to compare five representations of the
viral genome as inputs to thirteen machine learning classifiers predicting
sustained virological response after direct-acting antiviral therapy.

## Layout

```
data/       input data
code/       analysis notebooks and scripts
figures/    R scripts that generate the manuscript figures
results/    model outputs, including all five sensitivity analyses
```

## Data

| File | Contents |
|---|---|
| `data/data_RAS.xlsx` | 162 sequences with geno2pheno[HCV] resistance annotations, subtype assignment, regimen and outcome; the input to the RAS feature set |
| `data/full_gene_and_3_regions.xlsx` | aligned whole-genome and target-gene sequences; the input to the four sequence-based feature sets |
| `data/master_162.csv` | one row per sequence: accession, genotype, subtype, regimen, outcome |
| `data/genbank_metadata.csv` | GenBank record metadata for all 162 accessions, retrieved with `code/fetch_genbank_metadata.py` |
| `data/accession_provenance.csv` | each accession mapped to its source study and patient identifier; the basis of the patient-level analysis |

No participant-level clinical or demographic data accompany deposited sequences,
so none are present here.

## Reproducing the analysis

```bash
pip install -r requirements.txt
```

**Primary analysis.** Run the notebooks in `code/`. `Predictive_model_ras.ipynb`
covers the resistance feature set, `Predictive_model_3seq.ipynb` the target-gene
sets and `Predictive_model_full.ipynb` the whole-genome sets. Each writes into a
timestamped folder under `results/`.

**Sensitivity analysis 5, patient-level partitioning.** The 162 sequences come
from 146 patients; sixteen contributed two samples each, so a partition drawn
over sequences can place both samples of one patient on opposite sides of the
train/test boundary.

```bash
python code/fetch_genbank_metadata.py       # provenance; writes genbank_metadata.csv
python code/rerun_grouped.py                # partition drawn over patients
python code/rerun_ungrouped_multiseed.py    # sequence-level baseline, same seeds
```

`rerun_grouped.py` does not modify the notebook. It executes the notebook's own
function definitions and substitutes only the splitting function they call, so
every other step is identical to the primary analysis.

**Statistical analyses.** `code/p2_analysis.py` reproduces the subtype and
position-level comparisons; `code/firth.py` is the penalized logistic regression
used where complete separation makes maximum-likelihood estimates diverge.

## Figures

```r
install.packages(c("ggplot2","dplyr","tidyr","readxl","ragg","scales"))
setwd("figures")
source("Figure2_performance.R")   # A heatmap, B WG vs 3TG, C best AUC, D ROC
source("Figure3_SHAP.R")          # A importance, B shift plot, C1/C2, D covariation
source("Figure4_force.R")         # four individual prediction decompositions
```

Each panel is written three times: a vector PDF for the journal, a 600 dpi PNG
for figure assembly, and a 600 dpi uncompressed TIFF. Panel letters are added
during assembly and are deliberately not drawn by the scripts.

## Results

`results/model_runs/` holds the per-experiment outputs of the primary analysis:
metrics, bootstrap distributions, SHAP values and learning curves for every
model and feature set.

`results/sensitivity/` holds the patient-level partitioning analysis, the
matched sequence-level rerun over the same five seeds, and the comparison
between them.

## Environment

Analyses were run under Python 3.11 to 3.13 with the package versions pinned in
`requirements.txt`, and figures under R 4.4. Results are stable across these
versions to the third decimal, with the exception noted in the manuscript for
the quadratic discriminant model under patient-level partitioning.

## Citation

See `CITATION.cff`. Please cite the article rather than this repository alone.

## Licence

Code is released under the MIT licence (`LICENSE`). Data are released under
CC0 1.0; the underlying sequences remain publicly available from GenBank under
the accessions listed in `data/master_162.csv`.
