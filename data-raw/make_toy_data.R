# Regenerate the public toy data bundled in inst/extdata.
#
# Run this script from the package root:
#   Rscript data-raw/make_toy_data.R
#
# The toy data are synthetic. They are designed to exercise the package input
# readers, correlation-matrix checks, model fit, and LRT helpers. They are not
# intended to mimic a specific GWAS dataset.

set.seed(20260604)

n_snps <- 80L
n_outcomes <- 3L
n1 <- 100000L
n2 <- 90000L
true_beta <- c(0.16, -0.08, 0.06)
outcomes <- paste0("outcome_", seq_len(n_outcomes))

R <- matrix(
  c(
    1.00, 0.35, 0.10,
    0.35, 1.00, 0.25,
    0.10, 0.25, 1.00
  ),
  nrow = n_outcomes,
  byrow = TRUE,
  dimnames = list(outcomes, outcomes)
)

alleles <- data.frame(
  effect_allele = rep(c("A", "C", "G", "T"), length.out = n_snps),
  other_allele = rep(c("G", "T", "A", "C"), length.out = n_snps),
  stringsAsFactors = FALSE
)

gamma_hat <- rnorm(n_snps, mean = 0.04, sd = 0.010)
se_exposure <- runif(n_snps, min = 0.004, max = 0.008)

Gamma_hat <- vapply(seq_len(n_outcomes), function(j) {
  gamma_hat * true_beta[j] + rnorm(n_snps, mean = 0, sd = 0.006)
}, numeric(n_snps))
colnames(Gamma_hat) <- outcomes

se_outcome <- matrix(
  runif(n_snps * n_outcomes, min = 0.006, max = 0.010),
  nrow = n_snps,
  ncol = n_outcomes,
  dimnames = list(NULL, outcomes)
)

toy <- data.frame(
  SNP = paste0("rsToy", seq_len(n_snps)),
  alleles,
  beta_exposure = gamma_hat,
  se_exposure = se_exposure,
  exposure_sample_size = n1,
  outcome_sample_size = n2,
  stringsAsFactors = FALSE
)

for (outcome in outcomes) {
  toy[[paste0("beta_", outcome)]] <- Gamma_hat[, outcome]
  toy[[paste0("se_", outcome)]] <- se_outcome[, outcome]
}

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
write.csv(toy, "inst/extdata/toy_summary_stats.csv", row.names = FALSE, quote = FALSE)

R_csv <- data.frame(outcome = rownames(R), R, check.names = FALSE)
write.csv(R_csv, "inst/extdata/toy_outcome_correlation.csv", row.names = FALSE, quote = FALSE)

domain_map <- data.frame(
  outcome = outcomes,
  domain = c("domain_a", "domain_a", "domain_b"),
  stringsAsFactors = FALSE
)
write.csv(domain_map, "inst/extdata/toy_domain_map.csv", row.names = FALSE, quote = FALSE)

writeLines(
  c(
    "summary_stats: inst/extdata/toy_summary_stats.csv",
    "outcome_correlation: inst/extdata/toy_outcome_correlation.csv",
    "domain_map: inst/extdata/toy_domain_map.csv",
    "maxiter: 500",
    "rd: 1",
    "seed: 20260604",
    "n_snps: 80",
    "n_outcomes: 3",
    "exposure_sample_size: 100000",
    "outcome_sample_size: 90000",
    "true_beta:",
    "  - 0.16",
    "  - -0.08",
    "  - 0.06"
  ),
  "inst/extdata/toy_config.yml"
)
