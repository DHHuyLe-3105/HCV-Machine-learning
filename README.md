# Resistance-associated substitution profiles outperform whole-genome sequences for machine learning prediction of hepatitis C treatment failure

Data and analysis code for the study comparing five representations of the
hepatitis C virus genome as inputs to machine learning models predicting
sustained virological response after direct-acting antiviral therapy.

162 whole-genome sequences from 146 patients, retrieved from NCBI GenBank,
spanning 27 subtypes across genotypes 1, 2, 4 and 6.

## Layout

data/ input data and provenance
code/ analysis notebooks and scripts
figures/ R scripts that generate the manuscript figures
sensitivity/ outputs of the patient-level partitioning analysis


## data/

| File | Contents |
|---|---|
| `data_RAS.xlsx` | resistance annotations, subtype, regimen and outcome for all 162 sequences; the input to the RAS feature set |
| `full_gene_and_3_regions.xlsx` | aligned whole-genome and target-gene sequences; the input to the four sequence-based feature sets |
| `master_162.csv` | one row per sequence: accession, genotype, subtype, regimen, outcome |
| `genbank_metadata.csv` | GenBank record metadata for every accession |
| `accession_provenance.csv` | each accession mapped to its source study and patient identifier |

No participant-level clinical or demographic data accompany deposited sequences,
so none are present here.

## code/

`Predictive_model_ras.ipynb`, `Predictive_model_3seq.ipynb` and
`Predictive_model_full.ipynb` run the primary analysis for the resistance,
target-gene and whole-genome feature sets. `model_utils.py` holds the shared
helpers.

`fetch_genbank_metadata.py` retrieves the provenance table.
`rerun_grouped.py` repeats the pipeline with the train/test partition drawn over
patients rather than sequences; `rerun_ungrouped_multiseed.py` produces the
matched sequence-level baseline over the same seeds. Neither modifies the
notebooks: they execute the notebook's own function definitions and substitute
only the splitting function.

`firth.py` and `p2_analysis.py` reproduce the penalized regression and the
subtype and position-level comparisons.

```bash
pip install -r requirements.txt
python code/rerun_grouped.py
```

## figures/

```r
install.packages(c("ggplot2","dplyr","tidyr","readxl","ragg","scales"))
setwd("figures")
source("Figure2_performance.R")
source("Figure3_SHAP.R")
source("Figure4_force.R")
```

Each panel is written as a vector PDF, a 600 dpi PNG for figure assembly, and a
600 dpi TIFF. Panel letters are added during assembly and are not drawn by the
scripts. `Figure_Data_All.xlsx` and the three `shap_gbm_test_*.csv` files are the
inputs.

## sensitivity/

The 162 sequences come from 146 patients; sixteen contributed two samples each,
so a partition drawn over sequences can place both samples of one patient on
opposite sides of the train/test boundary. This folder holds the analysis that
redraws the partition over patients, the matched sequence-level rerun over the
same five seeds, and the comparison between them.

## Licence

Code under the MIT licence. Data under CC0 1.0; the underlying sequences remain
publicly available from GenBank under the accessions listed in
`data/master_162.csv`.

## Citation

See `CITATION.cff`. Please cite the article rather than this repository alone.

