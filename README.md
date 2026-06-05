# MR-MOSS

MR-MOSS is an R package for multi-outcome Mendelian randomization (MR) with
GWAS summary statistics. It fits one exposure against a correlated panel of
outcomes and reports:

- outcome-specific likelihood-ratio tests;
- optional prespecified subset/domain likelihood-ratio tests;
- one all-outcome global likelihood-ratio test.

MR-MOSS starts after instrument selection, LD clumping and allele harmonization.
The two workflows below are the intended first user-facing entry points.

## Install

```r
install.packages("remotes")
remotes::install_github("YunlongCao/MRMOSS")
```

MR-MOSS requires R 4.1 or later and a working C++ toolchain.

## 1. Run a built-in manuscript example

This is the fastest way to see what MR-MOSS expects and returns. The package
ships processed examples built from the MR-MOSS manuscript Supplementary Data:
negative controls, lipid cardiovascular controls and the MVP72 chronic disease
panel. These are not toy simulations. Each example already contains a
precomputed complete-case IV set, a harmonized exposure-by-multiple-outcome
summary-statistic matrix and a working outcome-correlation matrix.

### 1.1 List the available examples

```r
library(MRMOSS)

examples <- list_mrmoss_examples()
head(examples[, c(
  "example_id", "analysis", "exposure", "iv_threshold",
  "n_instruments", "n_outcomes"
)])
```

Filter by analysis type:

```r
list_mrmoss_examples(analysis = "cvd_positive_control")
list_mrmoss_examples(analysis = "mvp72", max_rows = 10)
list_mrmoss_examples(analysis = "negative_control", max_rows = 10)
```

A good first run is ApoB against the five cardiovascular positive-control
outcomes:

```r
example_id <- "cvd__Apolipoprotein_B_levels__p5e_08"
```

### 1.2 Load the prepared MR-MOSS input

```r
dat <- load_mrmoss_example(example_id)
dat
```

The loaded object exposes the actual input contract used by MR-MOSS:

```r
# SNP-exposure association vector
length(dat$summary_stats$gamma_hat)

# SNP-by-outcome association matrix
dim(dat$summary_stats$Gamma_hat)
dat$summary_stats$Gamma_hat[1:3, , drop = FALSE]

# Working outcome-correlation matrix
dat$R

# Optional domain or subset map
dat$domain_map
```

The same files can also be read directly:

```r
summary_stats <- read_mrmoss_summary_stats(dat$paths$summary_stats)
R <- read_outcome_correlation(dat$paths$outcome_correlation)
```

### 1.3 Fit MR-MOSS and inspect key results

For a concise first run, use the helper:

```r
res <- run_mrmoss_example(example_id)
```

By default this prints a console summary and writes no files. The returned
object contains the main result tables:

```r
res$global_lrt
head(res$outcome_lrt)
head(res$domain_lrt)
```

The helper is just a wrapper around the core MR-MOSS functions:

```r
fit <- fit_mrmoss(summary_stats = summary_stats, R = R)

global_lrt(fit)
outcome_lrt(fit)
domain_lrt(fit, dat$domain_map)
```

To export tables and a small HTML report, supply an output directory:

```r
res <- run_mrmoss_example(
  example_id,
  out_dir = "mrmoss_example_apob_cvd",
  open_report = FALSE
)
```

This writes:

```text
mrmoss_example_apob_cvd/
  input_qc.tsv
  outcome_lrt.tsv
  domain_lrt.tsv
  global_lrt.tsv
  mrmoss_fit.rds
  mrmoss_report.html
```

## 2. Start from OpenGWAS IDs with TwoSampleMR

Many MR users begin with `TwoSampleMR`, which can query OpenGWAS/IEU GWAS
records online. MR-MOSS can use this workflow by converting harmonized
TwoSampleMR output into MR-MOSS input.

Install TwoSampleMR if needed:

```r
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
```

Some OpenGWAS queries require current OpenGWAS authentication. If TwoSampleMR
asks for a JWT token or reports access/rate-limit errors, follow the OpenGWAS
login instructions and rerun the same code.

### 2.1 Choose one exposure and multiple outcomes

Search OpenGWAS records from R:

```r
library(TwoSampleMR)

ao <- available_outcomes()
subset(ao, grepl("body mass index", trait, ignore.case = TRUE))[1:5, ]
subset(ao, grepl("coronary heart disease", trait, ignore.case = TRUE))[1:5, ]
```

Set one exposure GWAS ID and at least two outcome GWAS IDs:

```r
exposure_id <- "ieu-a-2"
outcome_ids <- c("ieu-a-7", "ieu-a-26")
outcome_labels <- c("ieu-a-7" = "CHD", "ieu-a-26" = "T2D")
```

Replace these IDs with records appropriate for your analysis.

### 2.2 Extract instruments and harmonize alleles

```r
library(TwoSampleMR)
library(MRMOSS)

exposure_dat <- extract_instruments(
  outcomes = exposure_id,
  p1 = 5e-8,
  clump = TRUE
)

outcome_dat <- extract_outcome_data(
  snps = exposure_dat$SNP,
  outcomes = outcome_ids
)

harmonised <- harmonise_data(
  exposure_dat = exposure_dat,
  outcome_dat = outcome_dat,
  action = 2
)
```

`harmonised` is a long-format table with one row per SNP--outcome pair.

### 2.3 Convert the harmonized table to MR-MOSS input

```r
summary_wide <- twosamplemr_to_mrmoss_input(
  harmonised,
  outcome_names = outcome_labels,
  output_path = "summary_stats.tsv"
)

attr(summary_wide, "outcome_map")
dim(summary_wide)
```

The saved `summary_stats.tsv` has one complete-case instrument per row and
outcome-specific columns such as:

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

### 2.4 Estimate the outcome-correlation matrix

MR-MOSS also needs a working correlation matrix across the outcomes. Standard
TwoSampleMR/OpenGWAS calls retrieve associations for requested SNPs; they do
not generally return all genome-wide SNPs for a GWAS ID through this workflow.
To keep the first online workflow self-contained, MR-MOSS provides a built-in
candidate rsID panel. The user supplies only OpenGWAS IDs.

```r
candidate_snps <- mrmoss_null_snp_panel(
  n = 500,
  exclude = exposure_dat$SNP
)

null_outcome_dat <- extract_outcome_data(
  snps = candidate_snps,
  outcomes = outcome_ids,
  proxies = FALSE
)

R <- twosamplemr_to_outcome_correlation(
  null_outcome_dat,
  outcome_names = outcome_labels,
  pval_threshold = 1e-5,
  near_pd = TRUE,
  output_path = "outcome_correlation.tsv"
)
```

The example uses `n = 500` and `proxies = FALSE` so the online OpenGWAS query
usually finishes quickly. The `pval_threshold = 1e-5` filter is applied after
the SNP associations are returned; it removes variants strongly associated with
any fitted outcome before estimating the correlation matrix. For a final
manuscript-scale analysis, increase `n`, use a larger analysis-specific
candidate set, or estimate \(R\) from downloaded full GWAS summary statistics.
The code above is designed so that a new user does not need to prepare a
separate `null_snps.txt` file.

### 2.5 Fit MR-MOSS

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

If you have prespecified outcome domains or subsets:

```r
domain_map <- data.frame(
  outcome = c("CHD", "T2D"),
  domain = c("cardiometabolic", "cardiometabolic")
)

domain_results <- domain_lrt(fit, domain_map)
domain_results
```

Optional exports:

```r
write.table(outcome_results, "mrmoss_outcome_lrt.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(global_results, "mrmoss_global_lrt.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(domain_results, "mrmoss_domain_lrt.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(fit, "mrmoss_fit.rds")
```

## 3. Use downloaded GWAS Catalog files: ApoB and cardiovascular outcomes

This example starts from actual GWAS Catalog harmonized summary-statistic files
used in the MR-MOSS manuscript. The exposure is Apolipoprotein B levels
(`GCST90474314`). The outcome panel contains three MVP cardiovascular outcomes:
ischemic heart disease (`GCST90475929`), coronary atherosclerosis
(`GCST90475936`) and myocardial infarction (`GCST90480130`).

Unlike the TwoSampleMR online workflow, this route downloads the GWAS files once
and then works locally. It is the better pattern for manuscript-scale analyses
because estimating the outcome-correlation matrix no longer requires repeated
remote SNP queries.

The complete runnable version of this workflow is provided in
`scripts/07_gwascatalog_local_mrmoss_example.R`.

### 3.1 Download the GWAS Catalog harmonized files

The files below are GWAS Catalog harmonized `.tsv.gz` files. They already have
consistent harmonized variant fields such as `rsid`, `variant_id`,
`effect_allele`, `other_allele`, `p_value` and either `beta` or `odds_ratio`.

```r
dir.create("raw_gwas", showWarnings = FALSE)
dir.create("prepared", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

gwas_files <- data.frame(
  label = c("ApoB", "IHD", "CAD", "MI"),
  role = c("exposure", "outcome", "outcome", "outcome"),
  accession = c("GCST90474314", "GCST90475929", "GCST90475936", "GCST90480130"),
  sample_size_fallback = c(429941, 417274, 424341, 432053),
  url = c(
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90474001-GCST90475000/GCST90474314/harmonised/GCST90474314.h.tsv.gz",
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475929/harmonised/GCST90475929.h.tsv.gz",
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475936/harmonised/GCST90475936.h.tsv.gz",
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90480001-GCST90481000/GCST90480130/harmonised/GCST90480130.h.tsv.gz"
  ),
  stringsAsFactors = FALSE
)
gwas_files$path <- file.path("raw_gwas", paste0(gwas_files$accession, ".h.tsv.gz"))

for (i in seq_len(nrow(gwas_files))) {
  if (!file.exists(gwas_files$path[i])) {
    download.file(gwas_files$url[i], gwas_files$path[i], mode = "wb")
  }
}
```

These files are large. Use `data.table::fread()` to read them efficiently:

```r
install.packages("data.table")  # run once if needed
library(data.table)
library(MRMOSS)
```

Inspect one header:

```r
names(fread(gwas_files$path[gwas_files$label == "ApoB"], nrows = 0))
```

For this GWAS Catalog example, the ApoB exposure file contains `beta` and
`standard_error`. The cardiovascular outcome files contain `odds_ratio`,
confidence interval columns and sometimes missing `standard_error`; the code
below converts odds ratios to log odds ratios and derives standard errors from
the confidence interval when needed.

### 3.2 Convert the GWAS Catalog files to a simple local format

```r
standardize_gwascatalog_harmonized <- function(path, label, sample_size_fallback = NA_real_) {
  x <- data.table::fread(path, data.table = FALSE)

  se <- if ("standard_error" %in% names(x)) {
    suppressWarnings(as.numeric(x$standard_error))
  } else {
    rep(NA_real_, nrow(x))
  }
  if ("ci_upper" %in% names(x) && "ci_lower" %in% names(x)) {
    from_ci <- (log(as.numeric(x$ci_upper)) - log(as.numeric(x$ci_lower))) / (2 * 1.96)
    missing_se <- !is.finite(se)
    se[missing_se] <- from_ci[missing_se]
  }

  beta <- if ("beta" %in% names(x)) {
    as.numeric(x$beta)
  } else if ("odds_ratio" %in% names(x)) {
    log(as.numeric(x$odds_ratio))
  } else {
    stop("No beta or odds_ratio column found in ", path)
  }

  n <- if ("n" %in% names(x)) {
    suppressWarnings(as.numeric(x$n))
  } else {
    rep(sample_size_fallback, nrow(x))
  }
  n[!is.finite(n) | n <= 0] <- sample_size_fallback

  snp <- if ("rsid" %in% names(x)) as.character(x$rsid) else as.character(x$variant_id)

  out <- data.frame(
    SNP = snp,
    effect_allele = toupper(as.character(x$effect_allele)),
    other_allele = toupper(as.character(x$other_allele)),
    beta = beta,
    se = se,
    p_value = as.numeric(x$p_value),
    sample_size = n,
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$SNP) & nzchar(out$SNP), ]
  out <- out[!duplicated(out$SNP), ]
  out <- out[is.finite(out$beta) & is.finite(out$se) & out$se > 0, ]

  out_path <- file.path("prepared", paste0(label, "_standard.tsv"))
  write.table(out, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
  out_path
}

standard_paths <- setNames(
  mapply(
    standardize_gwascatalog_harmonized,
    path = gwas_files$path,
    label = gwas_files$label,
    sample_size_fallback = gwas_files$sample_size_fallback,
    SIMPLIFY = TRUE
  ),
  gwas_files$label
)
```

### 3.3 Select and LD-clump exposure instruments

MR-MOSS assumes the final instruments are approximately independent. For this
example, select ApoB instruments at \(P \leq 5\times10^{-8}\), then LD-clump
them. If you use TwoSampleMR clumping, the call below may use the OpenGWAS LD
service; alternatively, clump locally with PLINK and a matched ancestry
reference panel.

```r
exposure <- data.table::fread(standard_paths["ApoB"], data.table = FALSE)
iv_candidates <- exposure[exposure$p_value <= 5e-8, ]

if (requireNamespace("TwoSampleMR", quietly = TRUE)) {
  clump_input <- data.frame(
    SNP = iv_candidates$SNP,
    pval = iv_candidates$p_value,
    id.exposure = "ApoB",
    stringsAsFactors = FALSE
  )
  clumped <- TwoSampleMR::clump_data(
    clump_input,
    clump_kb = 10000,
    clump_r2 = 0.001,
    pop = "EUR"
  )
  exposure_iv <- iv_candidates[match(clumped$SNP, iv_candidates$SNP), ]
} else {
  stop("Install TwoSampleMR for this example, or replace this block with local PLINK clumping.")
}

write.table(exposure_iv, "prepared/ApoB_iv.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
```

### 3.4 Harmonize ApoB instruments against the three outcomes

This step writes the wide `summary_stats.tsv` input expected by
`read_mrmoss_summary_stats()`.

```r
is_palindromic <- function(a1, a2) {
  paste0(toupper(a1), toupper(a2)) %in% c("AT", "TA", "CG", "GC")
}

make_mrmoss_summary_stats <- function(exposure_path, outcome_paths, output_path,
                                      drop_palindromic = TRUE) {
  exposure <- data.table::fread(exposure_path, data.table = FALSE)
  if (drop_palindromic) {
    exposure <- exposure[!is_palindromic(exposure$effect_allele,
                                         exposure$other_allele), ]
  }

  wide <- data.frame(
    SNP = exposure$SNP,
    effect_allele = exposure$effect_allele,
    other_allele = exposure$other_allele,
    beta_exposure = exposure$beta,
    se_exposure = exposure$se,
    exposure_sample_size = exposure$sample_size,
    stringsAsFactors = FALSE
  )
  outcome_sample_sizes <- numeric()

  for (outcome_name in names(outcome_paths)) {
    outcome <- data.table::fread(outcome_paths[[outcome_name]], data.table = FALSE)
    if (drop_palindromic) {
      outcome <- outcome[!is_palindromic(outcome$effect_allele,
                                         outcome$other_allele), ]
    }
    outcome <- outcome[match(wide$SNP, outcome$SNP), ]

    same <- wide$effect_allele == outcome$effect_allele &
      wide$other_allele == outcome$other_allele
    flipped <- wide$effect_allele == outcome$other_allele &
      wide$other_allele == outcome$effect_allele
    same[is.na(same)] <- FALSE
    flipped[is.na(flipped)] <- FALSE

    beta <- rep(NA_real_, nrow(wide))
    se <- rep(NA_real_, nrow(wide))
    beta[same] <- outcome$beta[same]
    beta[flipped] <- -outcome$beta[flipped]
    se[same | flipped] <- outcome$se[same | flipped]

    wide[[paste0("beta_outcome_", outcome_name)]] <- beta
    wide[[paste0("se_outcome_", outcome_name)]] <- se
    outcome_sample_sizes <- c(
      outcome_sample_sizes,
      unique(outcome$sample_size[is.finite(outcome$sample_size)])
    )
  }

  outcome_cols <- grep("^beta_outcome_|^se_outcome_", names(wide), value = TRUE)
  wide <- wide[stats::complete.cases(wide[, outcome_cols, drop = FALSE]), ]
  wide$outcome_sample_size <- round(mean(outcome_sample_sizes, na.rm = TRUE))
  wide <- wide[, c(
    "SNP", "effect_allele", "other_allele",
    "beta_exposure", "se_exposure",
    "exposure_sample_size", "outcome_sample_size",
    outcome_cols
  )]
  if (!nrow(wide)) stop("No complete-case instruments remained after harmonization.")

  write.table(wide, output_path, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(wide)
}

outcome_paths <- standard_paths[c("IHD", "CAD", "MI")]

summary_wide <- make_mrmoss_summary_stats(
  exposure_path = "prepared/ApoB_iv.tsv",
  outcome_paths = outcome_paths,
  output_path = "prepared/summary_stats.tsv"
)
```

### 3.5 Estimate the outcome-correlation matrix from local null variants

Use variants present in all three outcome files, remove ApoB instruments and
retain variants with \(P>10^{-5}\) for every fitted outcome.

```r
estimate_local_R <- function(outcome_paths, exposure_iv_path, output_path,
                             pval_threshold = 1e-5,
                             max_snps = 50000) {
  exposure_iv <- data.table::fread(exposure_iv_path, data.table = FALSE)
  outcomes <- lapply(outcome_paths, data.table::fread, data.table = FALSE)

  common <- Reduce(intersect, lapply(outcomes, function(x) x$SNP))
  common <- setdiff(common, exposure_iv$SNP)
  if (length(common) > max_snps) common <- common[seq_len(max_snps)]

  z <- matrix(NA_real_, nrow = length(common), ncol = length(outcomes),
              dimnames = list(common, names(outcomes)))
  keep <- rep(TRUE, length(common))

  for (k in seq_along(outcomes)) {
    x <- outcomes[[k]]
    x <- x[match(common, x$SNP), ]
    z[, k] <- x$beta / x$se
    keep <- keep & is.finite(z[, k]) & is.finite(x$p_value) &
      x$p_value > pval_threshold
  }

  z <- z[keep, , drop = FALSE]
  if (nrow(z) < 2) stop("Too few complete-case putatively-null variants.")

  R <- estimate_outcome_correlation(z, near_pd = TRUE)
  write.table(data.frame(outcome = rownames(R), R, check.names = FALSE),
              output_path, sep = "\t", row.names = FALSE, quote = FALSE)
  R
}

R <- estimate_local_R(
  outcome_paths = outcome_paths,
  exposure_iv_path = "prepared/ApoB_iv.tsv",
  output_path = "prepared/outcome_correlation.tsv",
  pval_threshold = 1e-5
)
```

### 3.6 Fit MR-MOSS and export results

```r
summary_stats <- read_mrmoss_summary_stats("prepared/summary_stats.tsv")
R <- read_outcome_correlation("prepared/outcome_correlation.tsv")

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

domain_map <- data.frame(
  outcome = c("IHD", "CAD", "MI"),
  domain = c("cardiovascular", "cardiovascular", "cardiovascular")
)
domain_results <- domain_lrt(fit, domain_map)

write.table(outcome_results, "results/mrmoss_outcome_lrt.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(global_results, "results/mrmoss_global_lrt.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(domain_results, "results/mrmoss_domain_lrt.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(fit, "results/mrmoss_fit.rds")
```
