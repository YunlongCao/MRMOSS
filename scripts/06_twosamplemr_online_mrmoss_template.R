#!/usr/bin/env Rscript

# Template: online GWAS lookup through TwoSampleMR, followed by MR-MOSS.
#
# Edit the configuration block below before running. TwoSampleMR/OpenGWAS
# controls dataset availability, authentication and rate limits.

library(MRMOSS)

if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  stop(
    "TwoSampleMR is required for this template. Install it with:\n",
    "  install.packages('remotes')\n",
    "  remotes::install_github('MRCIEU/TwoSampleMR')"
  )
}

# -------------------------------------------------------------------------
# 1. Configuration
# -------------------------------------------------------------------------

# Replace these example IDs with GWAS records for your analysis. Use
# TwoSampleMR::available_outcomes() to search OpenGWAS from R.
exposure_id <- "ieu-a-2"
outcome_ids <- c(
  CHD = "ieu-a-7",
  T2D = "ieu-a-26"
)

instrument_p_threshold <- 5e-8
# Keep the default small enough for a quick OpenGWAS run. Increase this value,
# or estimate R from downloaded full summary statistics, for final analyses.
correlation_candidate_snps <- 500
results_dir <- "mrmoss_twosamplemr_results"

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
outcome_name_map <- stats::setNames(names(outcome_ids), unname(outcome_ids))

# -------------------------------------------------------------------------
# 2. Fetch instruments and outcome associations from OpenGWAS
# -------------------------------------------------------------------------

message("Extracting exposure instruments from OpenGWAS: ", exposure_id)
exposure_dat <- TwoSampleMR::extract_instruments(
  outcomes = exposure_id,
  p1 = instrument_p_threshold,
  clump = TRUE
)
if (!nrow(exposure_dat)) stop("No exposure instruments were returned")

message("Extracting outcome associations for ", length(outcome_ids), " outcomes")
outcome_dat <- TwoSampleMR::extract_outcome_data(
  snps = exposure_dat$SNP,
  outcomes = unname(outcome_ids)
)
if (!nrow(outcome_dat)) stop("No outcome associations were returned")

message("Harmonizing alleles with TwoSampleMR")
harmonised <- TwoSampleMR::harmonise_data(
  exposure_dat = exposure_dat,
  outcome_dat = outcome_dat,
  action = 2
)

# -------------------------------------------------------------------------
# 3. Convert harmonized TwoSampleMR output to MR-MOSS wide input
# -------------------------------------------------------------------------

summary_path <- file.path(results_dir, "summary_stats.tsv")
summary_wide <- twosamplemr_to_mrmoss_input(
  harmonised,
  outcome_names = outcome_name_map,
  output_path = summary_path
)

message("MR-MOSS complete-case instruments: ", nrow(summary_wide))
message("Outcome mapping:")
print(attr(summary_wide, "outcome_map"))

# -------------------------------------------------------------------------
# 4. Estimate the working outcome-correlation matrix
# -------------------------------------------------------------------------

cor_path <- file.path(results_dir, "outcome_correlation.tsv")

message("Estimating outcome correlation from the built-in candidate SNP panel")
candidate_snps <- mrmoss_null_snp_panel(
  n = correlation_candidate_snps,
  exclude = exposure_dat$SNP
)
null_outcome_dat <- TwoSampleMR::extract_outcome_data(
  snps = candidate_snps,
  outcomes = unname(outcome_ids),
  proxies = FALSE
)
R <- twosamplemr_to_outcome_correlation(
  null_outcome_dat,
  outcome_names = outcome_name_map,
  pval_threshold = 1e-5,
  near_pd = TRUE,
  output_path = cor_path
)

domain_path <- file.path(results_dir, "domain_map.tsv")
domain_map <- data.frame(
  outcome = names(outcome_ids),
  domain = "example_domain",
  stringsAsFactors = FALSE
)
write.table(domain_map, domain_path, sep = "\t", row.names = FALSE, quote = FALSE)

# -------------------------------------------------------------------------
# 5. Fit MR-MOSS and export results
# -------------------------------------------------------------------------

summary_stats <- read_mrmoss_summary_stats(summary_path)
R <- read_outcome_correlation(cor_path)
checked <- check_mrmoss_inputs(
  gamma_hat = summary_stats$gamma_hat,
  Gamma_hat = summary_stats$Gamma_hat,
  R = R,
  n1 = summary_stats$n1,
  n2 = summary_stats$n2
)
message("Outcome-correlation min eigenvalue: ", signif(checked$min_eigenvalue, 4))

fit <- fit_mrmoss(summary_stats = summary_stats, R = R)
outcome_results <- outcome_lrt(fit)
global_results <- global_lrt(fit)
domain_results <- domain_lrt(fit, domain_path)

write.table(outcome_results, file.path(results_dir, "mrmoss_outcome_lrt.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(global_results, file.path(results_dir, "mrmoss_global_lrt.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(domain_results, file.path(results_dir, "mrmoss_domain_lrt.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(fit, file.path(results_dir, "mrmoss_fit.rds"))

message("MR-MOSS results written to: ", normalizePath(results_dir))
print(global_results)
