# Online GWAS data with TwoSampleMR

This tutorial shows how to start from online GWAS records available through
TwoSampleMR/OpenGWAS and prepare the inputs needed by MR-MOSS.

This route is useful when a user wants to try MR-MOSS without manually
downloading GWAS Catalog or provider-specific summary-statistic files. The main
idea is:

```text
OpenGWAS IDs
  -> TwoSampleMR::extract_instruments()
  -> TwoSampleMR::extract_outcome_data()
  -> TwoSampleMR::harmonise_data()
  -> MRMOSS::twosamplemr_to_mrmoss_input()
  -> MRMOSS::fit_mrmoss()
```

MR-MOSS is still responsible only for the final multi-outcome likelihood. The
instrument selection, online lookup and allele harmonization steps are handled
by TwoSampleMR.

## Requirements

Install MR-MOSS and TwoSampleMR:

```r
install.packages("remotes")
remotes::install_github("YunlongCao/MRMOSS")
remotes::install_github("MRCIEU/TwoSampleMR")
```

TwoSampleMR queries the OpenGWAS service. Some datasets or high-volume queries
may require OpenGWAS authentication. Follow the current TwoSampleMR/OpenGWAS
instructions if your query asks for a JWT token or fails because of access
limits.

## Choose GWAS IDs

Use TwoSampleMR/OpenGWAS to identify one exposure GWAS and two or more outcome
GWAS records. For example:

```r
library(TwoSampleMR)

ao <- available_outcomes()
subset(ao, grepl("body mass index", trait, ignore.case = TRUE))[1:5, ]
subset(ao, grepl("coronary heart disease", trait, ignore.case = TRUE))[1:5, ]
```

For the template below, replace the example IDs with records that match your
scientific question, ancestry, genome build and data-use requirements.

## Extract and harmonize online data

```r
library(TwoSampleMR)
library(MRMOSS)

exposure_id <- "ieu-a-2"  # Example: BMI in many OpenGWAS tutorials.
outcome_ids <- c(
  CHD = "ieu-a-7",
  T2D = "ieu-a-26"
)
outcome_name_map <- stats::setNames(names(outcome_ids), unname(outcome_ids))

exposure_dat <- extract_instruments(
  outcomes = exposure_id,
  p1 = 5e-8,
  clump = TRUE
)

outcome_dat <- extract_outcome_data(
  snps = exposure_dat$SNP,
  outcomes = unname(outcome_ids)
)

harmonised <- harmonise_data(
  exposure_dat = exposure_dat,
  outcome_dat = outcome_dat,
  action = 2
)
```

The harmonized data frame is still in long TwoSampleMR format: one row per
SNP--outcome pair. Convert it to the wide MR-MOSS input:

```r
summary_wide <- twosamplemr_to_mrmoss_input(
  harmonised,
  outcome_names = outcome_name_map,
  output_path = "summary_stats.tsv"
)

attr(summary_wide, "outcome_map")
dim(summary_wide)
```

The output file contains one row per complete-case instrument and columns such
as:

```text
SNP
effect_allele
other_allele
beta_exposure
se_exposure
exposure_sample_size
outcome_sample_size
beta_outcome_CHD
se_outcome_CHD
beta_outcome_T2D
se_outcome_T2D
```

## Estimate an outcome-correlation matrix

MR-MOSS needs a working correlation matrix across outcomes. The package includes
a default candidate SNP panel so users do not need to create a separate null-SNP
file before trying the TwoSampleMR workflow.

Use the built-in panel, exclude the exposure instruments, query the same outcome
GWAS records and estimate \(R\):

```r
candidate_snps <- mrmoss_null_snp_panel(
  n = 500,
  exclude = exposure_dat$SNP
)

null_outcome_dat <- extract_outcome_data(
  snps = candidate_snps,
  outcomes = unname(outcome_ids),
  proxies = FALSE
)

R <- twosamplemr_to_outcome_correlation(
  null_outcome_dat,
  outcome_names = outcome_name_map,
  pval_threshold = 1e-5,
  near_pd = TRUE,
  output_path = "outcome_correlation.tsv"
)
```

This estimates the correlation from complete-case outcome z-scores,
`beta.outcome / se.outcome`, after removing rows with outcome association
`P <= 1e-5`. The example uses `n = 500` and `proxies = FALSE` to keep the
online OpenGWAS query short. The built-in SNP panel is a convenience
common-variant panel for workflow setup and small analyses. For a final
manuscript analysis, you may increase `n`, replace `candidate_snps` with a
larger analysis-specific candidate set, or estimate the matrix from local
downloaded summary statistics.

## Fit MR-MOSS

```r
summary_stats <- read_mrmoss_summary_stats("summary_stats.tsv")
R <- read_outcome_correlation("outcome_correlation.tsv")

check_mrmoss_inputs(
  gamma_hat = summary_stats$gamma_hat,
  Gamma_hat = summary_stats$Gamma_hat,
  R = R,
  n1 = summary_stats$n1,
  n2 = summary_stats$n2
)

fit <- fit_mrmoss(summary_stats = summary_stats, R = R)

outcome_results <- outcome_lrt(fit)
global_results <- global_lrt(fit)

outcome_results
global_results
```

If you have grouped outcomes, create a domain map:

```r
domain_map <- data.frame(
  outcome = c("CHD", "T2D"),
  domain = c("cardiometabolic", "cardiometabolic")
)
write.table(domain_map, "domain_map.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

domain_results <- domain_lrt(fit, "domain_map.tsv")
domain_results
```

## Full script

A runnable template is provided at:

```text
scripts/06_twosamplemr_online_mrmoss_template.R
```

Edit the GWAS IDs and output directory at the top of the script before running
it. OpenGWAS availability, authentication and rate limits are controlled by
TwoSampleMR/OpenGWAS, not by MR-MOSS.

## Notes for real analyses

- Confirm that exposure and outcome GWAS records use compatible ancestry and
  sample definitions.
- Record the exposure-instrument threshold and clumping settings used by
  `extract_instruments()`.
- Review harmonization messages from `harmonise_data()`, especially palindromic
  variants and dropped SNPs.
- The MR-MOSS reader uses scalar exposure and outcome sample sizes. If sample
  sizes vary across SNPs or outcomes, document the effective sample-size choice.
- The outcome-correlation matrix is part of the model input. The built-in
  candidate SNP panel is convenient for onboarding, but final analyses should
  document how candidate variants and putatively-null filters were chosen.
