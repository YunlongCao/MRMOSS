# MRMOSS

MRMOSS is an R package for MR-MOSS, a Mendelian randomization method for one exposure and many correlated outcomes. Using GWAS summary statistics, it estimates the causal effects on all outcomes in a joint model and tests them with likelihood ratio tests: for each outcome, for any prespecified set of outcomes, such as the diseases in a disease domain, and for all outcomes together.

## Installation

```r
# install.packages("remotes")
remotes::install_github("YunlongCao/MRMOSS")
```

R 4.1 or later and a C++ compiler are required.

## Input

MR-MOSS needs two inputs:

- a table with one row per instrument and the columns `SNP`, `effect_allele`, `other_allele`, `exposure_sample_size`, `outcome_sample_size`, `beta_exposure` and `beta_outcome_<name>` for each outcome. The associations are aligned to the effect allele and standardized as the z-score divided by the square root of the GWAS sample size.
- the correlation matrix of the outcomes, which can be estimated with `estimate_outcome_correlation()` from the z-scores of variants that are not associated with any of the outcomes.

## Example

```r
library(MRMOSS)

path <- system.file("extdata", "positive_control", package = "MRMOSS")
data <- read_mrmoss_summary_stats(file.path(path, "LPA_5e-08.tsv.gz"))
R <- read_outcome_correlation(file.path(path, "outcome_correlation.tsv"))

fit <- fit_mrmoss(data, R)
outcome_lrt(fit)
global_lrt(fit)
subset_lrt(fit, list(cardiac = c("CAD", "HF", "AF", "AS")))
```

## Positive control analysis

`inst/extdata/positive_control` contains the data of the positive control analysis in the MR-MOSS paper. The exposures are LDL cholesterol (`LDL`), systolic blood pressure (`SBP`), lifetime smoking (`SMOKING`), HbA1c (`HBA1C`) and lipoprotein(a) (`LPA`), with instruments selected at P < 5 × 10<sup>-8</sup> and at P < 5 × 10<sup>-7</sup>. The outcomes are nine diseases from FinnGen: coronary artery disease (`CAD`), ischemic stroke (`IS`), peripheral artery disease (`PAD`), abdominal aortic aneurysm (`AAA`), heart failure (`HF`), atrial fibrillation (`AF`), aortic stenosis (`AS`), chronic kidney disease (`CKD`) and type 2 diabetes (`T2D`). The outcome sample size is the mean of the nine FinnGen sample sizes.

The following script fits the ten combinations of exposure and threshold and prints the global P values next to the published values in Supplementary Table 7. It takes about ten minutes.

```r
source(system.file("scripts", "reproduce_positive_control.R", package = "MRMOSS"))
```

## License

MIT
