# Real GWAS Tutorial: From Summary Statistics To MR-MOSS Results

This tutorial is for users who already have access to real GWAS summary
statistics and want to run MR-MOSS on one exposure and a correlated outcome
panel. It follows a manuscript-style example:

- exposure: smoking initiation;
- outcomes: age-related macular degeneration outcomes, for example `AMD`,
  `AMD_dry`, and `AMD_wet`.

The repository does not include these real GWAS files. Use files from the
manuscript reproducibility archive if they are redistributable, or obtain the
source GWAS summary statistics from the original providers and respect their
data-use terms.

## 1. Create a local analysis folder

Work outside the software repository so large or restricted files are not
accidentally committed.

```text
my_mrmoss_analysis/
  raw_gwas/
    smoking_initiation.csv
    AMD.csv
    AMD_dry.csv
    AMD_wet.csv
  prepared/
  results/
```

Do not commit `raw_gwas/`, `prepared/`, `results/`, `*.rds`, compressed GWAS
files, PLINK files, or controlled-access files to a public GitHub repository.

## 2. Standardize each GWAS file

Before MR-MOSS, each exposure or outcome file used by the preparation template
should have these columns:

```text
SNP,effect_allele,other_allele,beta,se,sample_size
```

Optional exposure column:

```text
p_value
```

Important conventions:

- `beta` must be on the scale you intend to analyze. For binary traits, this is
  usually log odds ratio.
- `effect_allele` is the allele corresponding to the sign of `beta`.
- `SNP` IDs should use one genome build and one naming convention across files.
- Remove duplicate SNP IDs before harmonization or keep one pre-specified row.

If your source files use other names, rename them before using the template. For
example:

```r
gwas <- read.csv("raw_gwas/smoking_initiation.csv", check.names = FALSE)

standard <- data.frame(
  SNP = gwas$rsid,
  effect_allele = gwas$EA,
  other_allele = gwas$NEA,
  beta = gwas$beta,
  se = gwas$se,
  sample_size = gwas$N,
  p_value = gwas$p,
  stringsAsFactors = FALSE
)

write.csv(standard, "raw_gwas/smoking_initiation_standard.csv",
          row.names = FALSE)
```

Repeat this for each outcome file.

## 3. Select and LD-clump exposure instruments

MR-MOSS assumes the final rows are approximately independent instruments. It
does not perform LD clumping.

Use your preferred clumping tool, for example PLINK, TwoSampleMR, or another
validated pipeline. Record:

- exposure p-value threshold;
- LD reference panel and ancestry;
- genome build;
- LD `r2` threshold;
- clumping window size;
- how duplicated SNPs and proxies were handled.

After clumping, save the exposure file as:

```text
raw_gwas/smoking_initiation_clumped.csv
```

with the standardized columns above.

## 4. Create the wide MR-MOSS summary-statistic file

Copy the repository template into your analysis folder:

```sh
cp /path/to/MRMOSS/scripts/05_prepare_wide_input_template.R .
```

Edit the paths at the top:

```r
exposure_path <- "raw_gwas/smoking_initiation_clumped.csv"
outcome_paths <- c(
  AMD = "raw_gwas/AMD_standard.csv",
  AMD_dry = "raw_gwas/AMD_dry_standard.csv",
  AMD_wet = "raw_gwas/AMD_wet_standard.csv"
)
output_path <- "prepared/harmonized_summary_stats.csv"
```

Then run:

```sh
Rscript 05_prepare_wide_input_template.R
```

The template:

- matches outcome rows to the clumped exposure SNPs;
- keeps rows with the same allele orientation;
- flips outcome beta signs when alleles are reversed;
- never flips standard errors;
- drops unresolved mismatches and, by default, palindromic A/T or C/G rows;
- writes one wide row per retained instrument.

Review the resulting CSV manually. Its header should look like:

```text
SNP,effect_allele,other_allele,beta_exposure,se_exposure,exposure_sample_size,outcome_sample_size,beta_outcome_AMD,se_outcome_AMD,beta_outcome_AMD_dry,se_outcome_AMD_dry,beta_outcome_AMD_wet,se_outcome_AMD_wet
```

Outcome names parsed by MR-MOSS will be `AMD`, `AMD_dry`, and `AMD_wet`.

## 5. Build the outcome correlation matrix

MR-MOSS needs a working correlation matrix across outcomes. A practical approach
is to estimate it from null outcome-GWAS z-scores:

1. Start from variants shared across the outcome GWAS files.
2. Exclude the exposure instrument SNPs and nearby instrument loci.
3. Exclude variants strongly associated with one or more outcomes.
4. Compute `z = beta / se` for each outcome.
5. Keep enough nonconstant null variants to estimate a stable correlation.

Create:

```text
prepared/null_outcome_z_scores.csv
```

with this shape:

```text
SNP,AMD,AMD_dry,AMD_wet
rsNull1,0.10,-0.21,0.05
rsNull2,-0.33,-0.18,-0.27
```

Then estimate and save the matrix:

```r
library(MRMOSS)

null_z <- read.csv("prepared/null_outcome_z_scores.csv", check.names = FALSE)
rownames(null_z) <- null_z$SNP
null_z$SNP <- NULL

R <- estimate_outcome_correlation(null_z, near_pd = TRUE)

write.csv(
  data.frame(outcome = rownames(R), R, check.names = FALSE),
  "prepared/outcome_correlation.csv",
  row.names = FALSE
)
```

The saved file should look like:

```text
outcome,AMD,AMD_dry,AMD_wet
AMD,1,0.62,0.58
AMD_dry,0.62,1,0.41
AMD_wet,0.58,0.41,1
```

`near_pd = TRUE` is a numerical repair, not a scientific validation. If the
estimate needs heavy repair, revisit the null variant set or outcome sources.

## 6. Optional domain map

If you want grouped tests, create:

```text
prepared/domain_map.csv
```

Example:

```text
outcome,domain
AMD,overall_amd
AMD_dry,amd_subtype
AMD_wet,amd_subtype
```

The required columns are lowercase `outcome` and `domain`. Outcome values must
match the names parsed from the wide summary-statistic file.

## 7. Run MR-MOSS

Copy and edit the real-data analysis template:

```sh
cp /path/to/MRMOSS/scripts/04_template_real_gwas_analysis.R .
```

Set:

```r
summary_path <- "prepared/harmonized_summary_stats.csv"
cor_path <- "prepared/outcome_correlation.csv"
domain_path <- "prepared/domain_map.csv"
results_dir <- "results"
```

Run:

```sh
Rscript 04_template_real_gwas_analysis.R
```

The script writes:

```text
results/mrmoss_outcome_lrt.csv
results/mrmoss_global_lrt.csv
results/mrmoss_domain_lrt.csv
results/mrmoss_fit.rds
```

The `mrmoss_fit.rds` object can contain harmonized derived GWAS inputs. Keep it
local unless redistribution is allowed.

## 8. Interpret the result tables

`mrmoss_outcome_lrt.csv` has one row per outcome:

- `beta_hat`: fitted outcome-specific causal-effect estimate;
- `lrt`: likelihood-ratio statistic for that outcome;
- `p_value`: outcome-specific p-value;
- `minus_log10_p`: capped `-log10(P)` for plotting;
- `censored`: whether the reported `-log10(P)` hit the cap.

Add multiple-testing summaries if needed:

```r
out <- read.csv("results/mrmoss_outcome_lrt.csv")
out$fdr_bh <- p.adjust(out$p_value, method = "BH")
write.csv(out, "results/mrmoss_outcome_lrt_with_fdr.csv", row.names = FALSE)
```

Use the global and domain tests as panel-level summaries. They should not
replace outcome-specific inspection or MR sensitivity analyses.

## 9. Report the analysis

For a manuscript or report, document:

- source GWAS names, ancestry, genome build, sample sizes, and access terms;
- instrument selection threshold and LD-clumping settings;
- allele harmonization rules and palindromic SNP handling;
- final number of retained instruments;
- outcome correlation matrix source and null-variant filtering rules;
- scalar `n1` and `n2` strategy;
- MR-MOSS package version and GitHub release tag;
- whether harmonized inputs can be shared or require provider access.

This is the minimum information another user needs to understand and audit the
real-data analysis.
