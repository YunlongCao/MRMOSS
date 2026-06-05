#' Convert harmonized TwoSampleMR output to MR-MOSS wide input
#'
#' Converts the data frame returned by `TwoSampleMR::harmonise_data()` into the
#' wide summary-statistic format expected by [read_mrmoss_summary_stats()]. This
#' helper does not call TwoSampleMR or OpenGWAS directly; it only reshapes an
#' already harmonized data frame.
#'
#' @param harmonised_dat A data frame returned by `TwoSampleMR::harmonise_data()`
#'   for one exposure and multiple outcomes.
#' @param outcome_names Optional character vector giving MR-MOSS outcome labels.
#'   If named, names should match `id.outcome` values. If unnamed, values are
#'   used in the order of outcomes in `harmonised_dat`.
#' @param output_path Optional path to write the wide input file. `.csv`,
#'   `.tsv`, `.txt` and `.gz` paths are supported.
#' @return A data frame in MR-MOSS wide input format. The attribute
#'   `outcome_map` records the mapping from TwoSampleMR outcome IDs to MR-MOSS
#'   outcome labels.
#' @export
twosamplemr_to_mrmoss_input <- function(harmonised_dat, outcome_names = NULL,
                                        output_path = NULL) {
  if (!is.data.frame(harmonised_dat)) stop("harmonised_dat must be a data frame")
  dat <- harmonised_dat
  if ("mr_keep" %in% names(dat)) {
    keep <- dat$mr_keep
    if (is.character(keep)) keep <- toupper(keep) == "TRUE"
    dat <- dat[is.na(keep) | keep, , drop = FALSE]
  }
  required <- c(
    "SNP",
    "beta.exposure", "se.exposure",
    "beta.outcome", "se.outcome",
    "effect_allele.exposure", "other_allele.exposure"
  )
  missing <- setdiff(required, names(dat))
  if (length(missing)) stop("harmonised_dat is missing columns: ", paste(missing, collapse = ", "))
  id_col <- if ("id.outcome" %in% names(dat)) "id.outcome" else if ("outcome" %in% names(dat)) "outcome" else NULL
  if (is.null(id_col)) stop("harmonised_dat must contain id.outcome or outcome")
  dat[[id_col]] <- as.character(dat[[id_col]])
  outcome_ids <- unique(dat[[id_col]])
  if (!length(outcome_ids)) stop("No outcome rows are available after filtering")
  labels <- mrmoss_resolve_outcome_labels(dat, id_col, outcome_ids, outcome_names)

  key <- paste(dat$SNP, dat[[id_col]], sep = "\r")
  if (anyDuplicated(key)) {
    stop("harmonised_dat contains duplicate SNP/outcome rows after filtering")
  }
  snps_by_outcome <- split(as.character(dat$SNP), dat[[id_col]])
  common_snps <- Reduce(intersect, snps_by_outcome[outcome_ids])
  if (!length(common_snps)) stop("No complete-case SNPs are shared across all outcomes")
  first_order <- as.character(dat$SNP[dat[[id_col]] == outcome_ids[[1L]]])
  common_snps <- first_order[first_order %in% common_snps]

  n1_vec <- mrmoss_twosamplemr_sample_size(dat, "exposure")
  n2_vec <- mrmoss_twosamplemr_sample_size(dat, "outcome")
  dat$.__n1 <- n1_vec
  dat$.__n2 <- n2_vec

  base <- dat[match(common_snps, dat$SNP), , drop = FALSE]
  wide <- data.frame(
    SNP = common_snps,
    effect_allele = as.character(base$effect_allele.exposure),
    other_allele = as.character(base$other_allele.exposure),
    beta_exposure = as.numeric(base$beta.exposure),
    se_exposure = as.numeric(base$se.exposure),
    exposure_sample_size = as.numeric(base$.__n1),
    outcome_sample_size = NA_real_,
    stringsAsFactors = FALSE
  )

  n2_by_outcome <- matrix(NA_real_, nrow = length(common_snps), ncol = length(outcome_ids))
  for (k in seq_along(outcome_ids)) {
    id <- outcome_ids[[k]]
    label <- labels[[k]]
    sub <- dat[dat[[id_col]] == id, , drop = FALSE]
    sub <- sub[match(common_snps, sub$SNP), , drop = FALSE]
    if (any(is.na(sub$SNP))) stop("Internal complete-case alignment failed for outcome: ", id)
    wide[[paste0("beta_outcome_", label)]] <- as.numeric(sub$beta.outcome)
    wide[[paste0("se_outcome_", label)]] <- as.numeric(sub$se.outcome)
    n2_by_outcome[, k] <- as.numeric(sub$.__n2)
  }
  wide$outcome_sample_size <- round(rowMeans(n2_by_outcome, na.rm = TRUE))

  numeric_cols <- c("beta_exposure", "se_exposure", "exposure_sample_size",
                    "outcome_sample_size",
                    grep("^beta_outcome_|^se_outcome_", names(wide), value = TRUE))
  bad_numeric <- numeric_cols[!vapply(wide[numeric_cols], function(x) all(is.finite(x)), logical(1))]
  if (length(bad_numeric)) {
    stop("Converted MR-MOSS input contains non-finite values in: ",
         paste(bad_numeric, collapse = ", "))
  }
  if (any(wide$se_exposure <= 0) ||
      any(as.matrix(wide[grep("^se_outcome_", names(wide))]) <= 0)) {
    stop("Converted MR-MOSS input contains non-positive standard errors")
  }
  outcome_map <- data.frame(
    twosamplemr_outcome_id = outcome_ids,
    mrmoss_outcome = labels,
    stringsAsFactors = FALSE
  )
  if ("outcome" %in% names(dat)) {
    outcome_map$twosamplemr_outcome <- vapply(
      outcome_ids,
      function(id) as.character(dat$outcome[dat[[id_col]] == id][[1L]]),
      character(1)
    )
  }
  attr(wide, "outcome_map") <- outcome_map
  attr(wide, "n_input_rows") <- nrow(harmonised_dat)
  attr(wide, "n_retained_rows") <- nrow(dat)
  attr(wide, "n_complete_case_snps") <- nrow(wide)

  if (!is.null(output_path)) mrmoss_write_delimited(wide, output_path)
  wide
}

#' Return the built-in MR-MOSS candidate SNP panel for estimating outcome correlation
#'
#' Returns a built-in panel of common rsID variants that can be queried with
#' `TwoSampleMR::extract_outcome_data()` and then filtered as putatively null
#' variants by [twosamplemr_to_outcome_correlation()]. The panel is intended as a
#' convenient default for workflow setup and small analyses; users can replace it
#' with an analysis-specific null-variant set for final studies.
#'
#' @param n Number of SNPs to return. Use `Inf` to return all packaged SNPs.
#' @param exclude Optional vector of SNP IDs to remove, for example exposure
#'   instruments.
#' @return A character vector of rsIDs.
#' @export
mrmoss_null_snp_panel <- function(n = 5000L, exclude = NULL) {
  path <- system.file("extdata", "default_null_snp_panel.tsv.gz", package = "MRMOSS")
  if (!nzchar(path) || !file.exists(path)) stop("Built-in candidate SNP panel is not installed")
  dat <- mrmoss_read_table(path)
  if (!"SNP" %in% names(dat)) stop("Built-in candidate SNP panel is missing an SNP column")
  snps <- unique(as.character(dat$SNP))
  snps <- snps[nzchar(snps)]
  if (!is.null(exclude)) {
    snps <- setdiff(snps, as.character(exclude))
  }
  if (length(n) != 1L || is.na(n)) stop("n must be a positive integer or Inf")
  if (!is.infinite(n)) {
    n <- suppressWarnings(as.integer(n))
    if (!is.finite(n) || n < 1L) stop("n must be a positive integer or Inf")
    snps <- utils::head(snps, n)
  }
  snps
}

#' Estimate an MR-MOSS outcome-correlation matrix from TwoSampleMR outcome data
#'
#' Converts a data frame returned by `TwoSampleMR::extract_outcome_data()` for a
#' candidate SNP set into outcome z-score columns, filters on outcome
#' association strength when requested, then estimates the working
#' outcome-correlation matrix with
#' [estimate_outcome_correlation()].
#'
#' @param outcome_dat Data frame containing `SNP`, outcome identifier,
#'   `beta.outcome` and `se.outcome` columns.
#' @param outcome_names Optional labels passed to the MR-MOSS correlation matrix.
#'   The same mapping should be used in [twosamplemr_to_mrmoss_input()].
#' @param pval_threshold Optional threshold. If `pval.outcome` is present, rows
#'   with `pval.outcome <= pval_threshold` are excluded before estimating the
#'   matrix.
#' @param near_pd If `TRUE`, allow [estimate_outcome_correlation()] to project a
#'   non-positive-definite estimate to a nearest positive-definite correlation
#'   matrix.
#' @param output_path Optional path to write the correlation matrix.
#' @return A numeric outcome-correlation matrix.
#' @export
twosamplemr_to_outcome_correlation <- function(outcome_dat, outcome_names = NULL,
                                               pval_threshold = NULL,
                                               near_pd = TRUE,
                                               output_path = NULL) {
  if (!is.data.frame(outcome_dat)) stop("outcome_dat must be a data frame")
  dat <- outcome_dat
  required <- c("SNP", "beta.outcome", "se.outcome")
  missing <- setdiff(required, names(dat))
  if (length(missing)) stop("outcome_dat is missing columns: ", paste(missing, collapse = ", "))
  id_col <- if ("id.outcome" %in% names(dat)) "id.outcome" else if ("outcome" %in% names(dat)) "outcome" else NULL
  if (is.null(id_col)) stop("outcome_dat must contain id.outcome or outcome")
  dat[[id_col]] <- as.character(dat[[id_col]])
  if (!is.null(pval_threshold)) {
    if (!"pval.outcome" %in% names(dat)) stop("pval_threshold requires a pval.outcome column")
    dat <- dat[suppressWarnings(as.numeric(dat$pval.outcome)) > pval_threshold, , drop = FALSE]
  }
  dat$beta.outcome <- suppressWarnings(as.numeric(dat$beta.outcome))
  dat$se.outcome <- suppressWarnings(as.numeric(dat$se.outcome))
  dat <- dat[is.finite(dat$beta.outcome) & is.finite(dat$se.outcome) & dat$se.outcome > 0, , drop = FALSE]
  outcome_ids <- unique(dat[[id_col]])
  if (!length(outcome_ids)) stop("No outcome rows are available for correlation estimation")
  labels <- mrmoss_resolve_outcome_labels(dat, id_col, outcome_ids, outcome_names)
  key <- paste(dat$SNP, dat[[id_col]], sep = "\r")
  if (anyDuplicated(key)) stop("outcome_dat contains duplicate SNP/outcome rows")
  snps_by_outcome <- split(as.character(dat$SNP), dat[[id_col]])
  common_snps <- Reduce(intersect, snps_by_outcome[outcome_ids])
  if (length(common_snps) < 2L) stop("At least two complete-case null SNPs are required")
  first_order <- as.character(dat$SNP[dat[[id_col]] == outcome_ids[[1L]]])
  common_snps <- first_order[first_order %in% common_snps]
  z <- matrix(NA_real_, nrow = length(common_snps), ncol = length(outcome_ids),
              dimnames = list(common_snps, labels))
  for (k in seq_along(outcome_ids)) {
    id <- outcome_ids[[k]]
    sub <- dat[dat[[id_col]] == id, , drop = FALSE]
    sub <- sub[match(common_snps, sub$SNP), , drop = FALSE]
    z[, k] <- sub$beta.outcome / sub$se.outcome
  }
  R <- estimate_outcome_correlation(z, near_pd = near_pd)
  if (!is.null(output_path)) {
    mrmoss_write_delimited(data.frame(outcome = rownames(R), R, check.names = FALSE), output_path)
  }
  R
}

mrmoss_resolve_outcome_labels <- function(dat, id_col, outcome_ids, outcome_names = NULL) {
  if (is.null(outcome_names)) {
    raw <- if ("outcome" %in% names(dat)) {
      vapply(outcome_ids, function(id) as.character(dat$outcome[dat[[id_col]] == id][[1L]]), character(1))
    } else {
      outcome_ids
    }
  } else {
    if (!is.character(outcome_names)) stop("outcome_names must be a character vector")
    if (!is.null(names(outcome_names)) && all(outcome_ids %in% names(outcome_names))) {
      raw <- as.character(outcome_names[outcome_ids])
    } else {
      if (length(outcome_names) != length(outcome_ids)) {
        stop("Unnamed outcome_names must have one value per outcome")
      }
      raw <- as.character(outcome_names)
    }
  }
  labels <- gsub("[^A-Za-z0-9]+", "_", raw)
  labels <- gsub("^_+|_+$", "", labels)
  labels[!nzchar(labels)] <- paste0("outcome_", which(!nzchar(labels)))
  labels[grepl("^[0-9]", labels)] <- paste0("outcome_", labels[grepl("^[0-9]", labels)])
  make.unique(labels, sep = "_")
}

mrmoss_twosamplemr_sample_size <- function(dat, role) {
  direct <- paste0("samplesize.", role)
  case_col <- paste0("ncase.", role)
  control_col <- paste0("ncontrol.", role)
  if (direct %in% names(dat)) {
    n <- suppressWarnings(as.numeric(dat[[direct]]))
  } else {
    n <- rep(NA_real_, nrow(dat))
  }
  if (any(!is.finite(n)) && all(c(case_col, control_col) %in% names(dat))) {
    n_case <- suppressWarnings(as.numeric(dat[[case_col]]))
    n_control <- suppressWarnings(as.numeric(dat[[control_col]]))
    from_case_control <- n_case + n_control
    replace <- !is.finite(n) & is.finite(from_case_control)
    n[replace] <- from_case_control[replace]
  }
  if (!all(is.finite(n)) || any(n <= 0)) {
    stop("Could not derive positive sample sizes for ", role,
         ". Add samplesize.", role, " or ncase.", role, "/ncontrol.", role,
         " columns before conversion.")
  }
  n
}

mrmoss_write_delimited <- function(x, output_path) {
  if (!is.character(output_path) || length(output_path) != 1L || !nzchar(output_path)) {
    stop("output_path must be a non-empty scalar path")
  }
  lower <- tolower(output_path)
  is_tab <- grepl("\\.(tsv|tab|txt)(\\.gz)?$", lower)
  sep <- if (is_tab) "\t" else ","
  con <- if (grepl("\\.gz$", lower)) gzfile(output_path, open = "wt") else file(output_path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, con, sep = sep, row.names = FALSE, quote = FALSE, na = "")
  invisible(output_path)
}
