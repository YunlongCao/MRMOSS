#!/usr/bin/env Rscript

# Template: downloaded GWAS Catalog harmonized summary statistics followed by
# MR-MOSS. This example uses one manuscript positive-control exposure, ApoB
# (GCST90474314), and three MVP cardiovascular outcomes.

library(MRMOSS)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Install data.table first: install.packages('data.table')")
}
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  stop("Install TwoSampleMR for LD clumping, or replace the clumping block with local PLINK clumping.")
}

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
    message("Downloading ", gwas_files$accession[i])
    download.file(gwas_files$url[i], gwas_files$path[i], mode = "wb")
  }
}

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

message("Selecting and LD-clumping ApoB instruments")
exposure <- data.table::fread(standard_paths["ApoB"], data.table = FALSE)
iv_candidates <- exposure[exposure$p_value <= 5e-8, ]
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
write.table(exposure_iv, "prepared/ApoB_iv.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

is_palindromic <- function(a1, a2) {
  paste0(toupper(a1), toupper(a2)) %in% c("AT", "TA", "CG", "GC")
}

make_mrmoss_summary_stats <- function(exposure_path, outcome_paths, output_path,
                                      drop_palindromic = TRUE) {
  exposure <- data.table::fread(exposure_path, data.table = FALSE)
  if (drop_palindromic) {
    exposure <- exposure[!is_palindromic(exposure$effect_allele, exposure$other_allele), ]
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
      outcome <- outcome[!is_palindromic(outcome$effect_allele, outcome$other_allele), ]
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
make_mrmoss_summary_stats(
  exposure_path = "prepared/ApoB_iv.tsv",
  outcome_paths = outcome_paths,
  output_path = "prepared/summary_stats.tsv"
)

estimate_local_R <- function(outcome_paths, exposure_iv_path, output_path,
                             pval_threshold = 1e-5, max_snps = 50000) {
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
    keep <- keep & is.finite(z[, k]) & is.finite(x$p_value) & x$p_value > pval_threshold
  }

  z <- z[keep, , drop = FALSE]
  if (nrow(z) < 2) stop("Too few complete-case putatively-null variants.")

  R <- estimate_outcome_correlation(z, near_pd = TRUE)
  write.table(data.frame(outcome = rownames(R), R, check.names = FALSE),
              output_path, sep = "\t", row.names = FALSE, quote = FALSE)
  R
}

estimate_local_R(
  outcome_paths = outcome_paths,
  exposure_iv_path = "prepared/ApoB_iv.tsv",
  output_path = "prepared/outcome_correlation.tsv"
)

summary_stats <- read_mrmoss_summary_stats("prepared/summary_stats.tsv")
R <- read_outcome_correlation("prepared/outcome_correlation.tsv")
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
domain_results <- domain_lrt(
  fit,
  data.frame(
    outcome = c("IHD", "CAD", "MI"),
    domain = "cardiovascular",
    stringsAsFactors = FALSE
  )
)

write.table(outcome_results, "results/mrmoss_outcome_lrt.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(global_results, "results/mrmoss_global_lrt.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(domain_results, "results/mrmoss_domain_lrt.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(fit, "results/mrmoss_fit.rds")

print(global_results)
