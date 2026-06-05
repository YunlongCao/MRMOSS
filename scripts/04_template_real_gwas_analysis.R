library(MRMOSS)

# Edit these paths before running.
summary_path <- "path/to/harmonized_summary_stats.csv"
cor_path <- "path/to/outcome_correlation.csv"
domain_path <- "path/to/domain_map.csv"  # Set to NA_character_ if not used.
results_dir <- "mrmoss_results"

# The saved fit object can contain harmonized, derived GWAS inputs. Keep the
# results directory local unless you have confirmed that these derived files can
# be redistributed under the original GWAS data-use terms.

if (!file.exists(summary_path)) {
  stop("Edit summary_path so it points to your harmonized MR-MOSS summary-statistic CSV.")
}
if (!file.exists(cor_path)) {
  stop("Edit cor_path so it points to your outcome correlation CSV.")
}
if (!is.na(domain_path) && !file.exists(domain_path)) {
  stop("Edit domain_path or set it to NA_character_ if you do not have domain groups.")
}

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

is_palindromic <- function(a1, a2) {
  paste0(toupper(a1), toupper(a2)) %in% c("AT", "TA", "CG", "GC")
}

raw_summary <- read.csv(summary_path, check.names = FALSE)
message("Rows in summary file: ", nrow(raw_summary))
message("Duplicated SNP IDs: ", sum(duplicated(raw_summary$SNP)))
bad_alleles <- !toupper(raw_summary$effect_allele) %in% c("A", "C", "G", "T") |
  !toupper(raw_summary$other_allele) %in% c("A", "C", "G", "T")
message("Rows with non-ACGT alleles: ", sum(bad_alleles, na.rm = TRUE))
message("Palindromic A/T or C/G rows: ",
        sum(is_palindromic(raw_summary$effect_allele, raw_summary$other_allele),
            na.rm = TRUE))

summary_stats <- read_mrmoss_summary_stats(summary_path)
R <- read_outcome_correlation(cor_path)

message("Parsed outcomes: ", paste(summary_stats$outcomes, collapse = ", "))
if (length(unique(summary_stats$exposure_sample_size)) > 1L) {
  message("Exposure sample size varies across rows; MR-MOSS will use n1 = ",
          summary_stats$n1)
}
if (length(unique(summary_stats$outcome_sample_size)) > 1L) {
  message("Outcome sample size varies across rows; MR-MOSS will use n2 = ",
          summary_stats$n2)
}

missing_in_R <- setdiff(summary_stats$outcomes, rownames(R))
extra_in_R <- setdiff(rownames(R), summary_stats$outcomes)
if (length(missing_in_R) || length(extra_in_R)) {
  stop("Outcome names differ between summary stats and R. Missing in R: ",
       paste(missing_in_R, collapse = ", "),
       "; extra in R: ", paste(extra_in_R, collapse = ", "))
}

if (!is.na(domain_path)) {
  domain_map <- read.csv(domain_path, check.names = FALSE)
  required_domain_cols <- c("outcome", "domain")
  missing_domain_cols <- setdiff(required_domain_cols, names(domain_map))
  if (length(missing_domain_cols)) {
    stop("Domain map is missing columns: ", paste(missing_domain_cols, collapse = ", "))
  }
  duplicated_domain_rows <- duplicated(domain_map[, required_domain_cols])
  if (any(duplicated_domain_rows)) stop("Domain map contains duplicated outcome/domain rows.")
  unknown_domain_outcomes <- setdiff(domain_map$outcome, summary_stats$outcomes)
  if (length(unknown_domain_outcomes)) {
    stop("Domain map contains outcomes not present in the fit: ",
         paste(unknown_domain_outcomes, collapse = ", "))
  }
  unmapped_outcomes <- setdiff(summary_stats$outcomes, domain_map$outcome)
  if (length(unmapped_outcomes)) {
    message("Outcomes without domain labels: ", paste(unmapped_outcomes, collapse = ", "))
  }
}

checked <- check_mrmoss_inputs(
  gamma_hat = summary_stats$gamma_hat,
  Gamma_hat = summary_stats$Gamma_hat,
  R = R,
  n1 = summary_stats$n1,
  n2 = summary_stats$n2
)

message("Outcome correlation matrix minimum eigenvalue: ",
        signif(checked$min_eigenvalue, 4))

fit <- fit_mrmoss(summary_stats = summary_stats, R = R)

outcome_results <- outcome_lrt(fit)
global_results <- global_lrt(fit)

write.csv(outcome_results, file.path(results_dir, "mrmoss_outcome_lrt.csv"),
          row.names = FALSE)
write.csv(global_results, file.path(results_dir, "mrmoss_global_lrt.csv"),
          row.names = FALSE)

if (!is.na(domain_path)) {
  domain_results <- domain_lrt(fit, domain_path)
  write.csv(domain_results, file.path(results_dir, "mrmoss_domain_lrt.csv"),
            row.names = FALSE)
}

saveRDS(fit, file.path(results_dir, "mrmoss_fit.rds"))
message("MR-MOSS results written to: ", normalizePath(results_dir))
