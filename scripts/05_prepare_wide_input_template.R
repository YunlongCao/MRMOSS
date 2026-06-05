# Template for creating the wide MR-MOSS input CSV from one exposure GWAS file
# and several outcome GWAS files.
#
# This script assumes each input file has already been restricted to the
# intended genome build and that exposure instruments have been selected and
# LD-clumped upstream. It performs only basic allele matching/flipping and
# should be reviewed before use in a manuscript analysis.

exposure_path <- "path/to/exposure_gwas_clumped.csv"
outcome_paths <- c(
  CAD = "path/to/cad_outcome_gwas.csv",
  STROKE = "path/to/stroke_outcome_gwas.csv"
)
output_path <- "harmonized_summary_stats.csv"

# Expected columns in each input file after any user-side renaming:
# SNP,effect_allele,other_allele,beta,se,sample_size
# The exposure file can also contain p_value; if present, p_value_threshold is
# applied after reading.
p_value_threshold <- 5e-8
drop_palindromic <- TRUE

required_cols <- c("SNP", "effect_allele", "other_allele", "beta", "se", "sample_size")

read_gwas <- function(path, label) {
  if (!file.exists(path)) stop("File not found for ", label, ": ", path)
  dat <- read.csv(path, check.names = FALSE)
  missing <- setdiff(required_cols, names(dat))
  if (length(missing)) {
    stop(label, " file is missing columns: ", paste(missing, collapse = ", "))
  }
  dat <- dat[!duplicated(dat$SNP), , drop = FALSE]
  dat$effect_allele <- toupper(dat$effect_allele)
  dat$other_allele <- toupper(dat$other_allele)
  dat$beta <- as.numeric(dat$beta)
  dat$se <- as.numeric(dat$se)
  dat$sample_size <- as.numeric(dat$sample_size)
  dat
}

is_palindromic <- function(a1, a2) {
  pair <- paste0(a1, a2)
  pair %in% c("AT", "TA", "CG", "GC")
}

exposure <- read_gwas(exposure_path, "exposure")
if ("p_value" %in% names(exposure)) {
  exposure <- exposure[as.numeric(exposure$p_value) <= p_value_threshold, , drop = FALSE]
}
if (drop_palindromic) {
  exposure <- exposure[!is_palindromic(exposure$effect_allele, exposure$other_allele), , drop = FALSE]
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

outcome_sample_sizes <- list()

for (outcome_name in names(outcome_paths)) {
  outcome <- read_gwas(outcome_paths[[outcome_name]], outcome_name)
  if (drop_palindromic) {
    outcome <- outcome[!is_palindromic(outcome$effect_allele, outcome$other_allele), , drop = FALSE]
  }

  outcome <- outcome[match(wide$SNP, outcome$SNP), , drop = FALSE]

  same <- wide$effect_allele == outcome$effect_allele &
    wide$other_allele == outcome$other_allele
  flipped <- wide$effect_allele == outcome$other_allele &
    wide$other_allele == outcome$effect_allele
  same[is.na(same)] <- FALSE
  flipped[is.na(flipped)] <- FALSE
  keep <- same | flipped

  beta <- rep(NA_real_, nrow(wide))
  beta[same] <- outcome$beta[same]
  beta[flipped] <- -outcome$beta[flipped]

  se <- outcome$se
  se[!keep] <- NA_real_

  wide[[paste0("beta_outcome_", outcome_name)]] <- beta
  wide[[paste0("se_outcome_", outcome_name)]] <- se
  outcome_sample_sizes[[outcome_name]] <- outcome$sample_size
}

complete_cols <- grep("^beta_outcome_|^se_outcome_", names(wide), value = TRUE)
wide <- wide[stats::complete.cases(wide[, complete_cols, drop = FALSE]), , drop = FALSE]

all_outcome_n <- unlist(outcome_sample_sizes, use.names = FALSE)
effective_outcome_n <- round(mean(all_outcome_n[is.finite(all_outcome_n)], na.rm = TRUE))
wide$outcome_sample_size <- effective_outcome_n

wide <- wide[, c(
  "SNP", "effect_allele", "other_allele",
  "beta_exposure", "se_exposure",
  "exposure_sample_size", "outcome_sample_size",
  complete_cols
)]

if (!nrow(wide)) stop("No SNPs remained after harmonization and complete-case filtering.")

write.csv(wide, output_path, row.names = FALSE, quote = FALSE)
message("Wrote ", nrow(wide), " harmonized SNP rows to ", output_path)
message("Using scalar outcome_sample_size = ", effective_outcome_n,
        ". Confirm this effective sample size is appropriate for your analysis.")
