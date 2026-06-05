#' MR-MOSS package
#'
#' MR-MOSS provides R wrappers around the likelihood-based C++ core for
#' correlated-outcome Mendelian randomization with summary statistics.
#'
#' @keywords internal
#' @useDynLib MRMOSS, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom Matrix nearPD
#' @importFrom methods is
#' @importFrom stats cor
#' @importFrom utils read.csv
"_PACKAGE"

#' Read a working outcome correlation matrix
#'
#' Reads a CSV file whose first column contains outcome names and whose remaining
#' columns contain a square outcome-by-outcome correlation matrix.
#'
#' @param path Path to a CSV correlation matrix.
#' @return A numeric matrix with row names taken from the first CSV column.
#' @export
read_outcome_correlation <- function(path) {
  dat <- mrmoss_read_table(path)
  if (ncol(dat) < 2L) stop("R matrix CSV must have an outcome column plus matrix columns")
  rn <- as.character(dat[[1L]])
  if (anyNA(rn) || any(!nzchar(rn))) stop("R matrix row names must be non-empty outcome names")
  if (anyDuplicated(rn)) stop("R matrix row names must be unique")
  R <- as.matrix(dat[, -1L, drop = FALSE])
  storage.mode(R) <- "numeric"
  if (nrow(R) != ncol(R)) stop("R matrix must be square after removing the outcome-name column")
  if (is.null(colnames(R)) || any(!nzchar(colnames(R)))) {
    stop("R matrix columns must be named with outcome names")
  }
  if (anyDuplicated(colnames(R))) stop("R matrix column names must be unique")
  rownames(R) <- rn
  if (!setequal(rownames(R), colnames(R))) {
    stop("R matrix row names and column names must contain the same outcomes")
  }
  R <- R[rownames(R), rownames(R), drop = FALSE]
  if (!all(is.finite(R))) stop("R matrix must contain only finite numeric values")
  R
}

#' Read wide MR-MOSS summary statistics
#'
#' Reads a toy or user-prepared wide summary-statistic CSV with one row per SNP,
#' SNP-exposure association columns, and one beta column per outcome. Outcome
#' columns can be named `beta_outcome_CAD` with matching `se_outcome_CAD`, or
#' `beta_CAD` with matching `se_CAD`.
#'
#' @param path Path to a CSV summary-statistic file.
#' @return A list containing `gamma_hat`, `Gamma_hat`, standard-error vectors
#'   and matrices, sample-size vectors and scalar means, outcome names, and the
#'   raw data frame.
#' @export
read_mrmoss_summary_stats <- function(path) {
  dat <- mrmoss_read_table(path)
  read_cols <- c(
    "SNP", "effect_allele", "other_allele",
    "beta_exposure", "se_exposure",
    "exposure_sample_size", "outcome_sample_size"
  )
  missing <- setdiff(read_cols, names(dat))
  if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))
  outcome_cols <- grep("^beta_outcome_|^beta_", names(dat), value = TRUE)
  outcome_cols <- setdiff(outcome_cols, "beta_exposure")
  if (!length(outcome_cols)) stop("No outcome beta columns found")
  outcomes <- ifelse(
    grepl("^beta_outcome_", outcome_cols),
    sub("^beta_outcome_", "", outcome_cols),
    sub("^beta_", "", outcome_cols)
  )
  outcomes <- ifelse(grepl("^[0-9]+$", outcomes), paste0("outcome_", outcomes), outcomes)
  if (anyDuplicated(outcomes)) stop("Outcome beta columns map to duplicated outcome names")
  se_cols <- sub("^beta_", "se_", outcome_cols)
  missing_se <- setdiff(se_cols, names(dat))
  if (length(missing_se)) {
    stop("Missing required outcome standard-error columns: ", paste(missing_se, collapse = ", "))
  }
  gamma <- suppressWarnings(as.numeric(dat$beta_exposure))
  Gamma <- as.matrix(dat[, outcome_cols, drop = FALSE])
  storage.mode(Gamma) <- "numeric"
  colnames(Gamma) <- outcomes
  se_Gamma <- as.matrix(dat[, se_cols, drop = FALSE])
  storage.mode(se_Gamma) <- "numeric"
  colnames(se_Gamma) <- outcomes
  se_gamma <- suppressWarnings(as.numeric(dat$se_exposure))
  n1_vec <- suppressWarnings(as.numeric(dat$exposure_sample_size))
  n2_vec <- suppressWarnings(as.numeric(dat$outcome_sample_size))
  if (!all(is.finite(gamma)) || !all(is.finite(Gamma))) {
    stop("Exposure and outcome beta columns must contain only finite numeric values")
  }
  if (!all(is.finite(se_gamma)) || !all(is.finite(se_Gamma)) ||
      any(se_gamma <= 0) || any(se_Gamma <= 0)) {
    stop("Exposure and outcome standard-error columns must be finite and positive")
  }
  if (!all(is.finite(n1_vec)) || !all(is.finite(n2_vec)) ||
      any(n1_vec <= 0) || any(n2_vec <= 0)) {
    stop("Sample-size columns must be finite and positive")
  }
  if (length(unique(n1_vec)) > 1L) {
    warning("exposure_sample_size varies across SNPs; using the rounded mean as scalar n1")
  }
  if (length(unique(n2_vec)) > 1L) {
    warning("outcome_sample_size varies across SNPs; using the rounded mean as scalar n2")
  }
  n1 <- round(mean(n1_vec, na.rm = TRUE))
  n2 <- round(mean(n2_vec, na.rm = TRUE))
  list(
    snp = as.character(dat$SNP),
    effect_allele = as.character(dat$effect_allele),
    other_allele = as.character(dat$other_allele),
    gamma_hat = gamma,
    Gamma_hat = Gamma,
    se_gamma = se_gamma,
    se_Gamma = se_Gamma,
    exposure_sample_size = n1_vec,
    outcome_sample_size = n2_vec,
    n1 = n1,
    n2 = n2,
    outcomes = outcomes,
    raw = dat
  )
}

#' Check MR-MOSS model inputs
#'
#' Validates dimensions, numeric finiteness, sample sizes, outcome-name
#' alignment, symmetry, and positive definiteness of the working outcome
#' correlation matrix before calling the C++ fitting core.
#'
#' @param gamma_hat Numeric vector of SNP-exposure association estimates.
#' @param Gamma_hat Numeric matrix of SNP-outcome association estimates, with
#'   SNPs in rows and outcomes in columns.
#' @param R Numeric outcome correlation matrix.
#' @param n1 Exposure GWAS sample size.
#' @param n2 Outcome GWAS sample size.
#' @return A normalized input list, with `R` reordered to match `Gamma_hat`
#'   outcome columns when names are available.
#' @export
check_mrmoss_inputs <- function(gamma_hat, Gamma_hat, R, n1, n2) {
  gamma_hat <- as.numeric(gamma_hat)
  Gamma_hat <- as.matrix(Gamma_hat)
  storage.mode(Gamma_hat) <- "numeric"
  R <- as.matrix(R)
  storage.mode(R) <- "numeric"
  n1 <- suppressWarnings(as.numeric(n1))
  n2 <- suppressWarnings(as.numeric(n2))
  if (length(gamma_hat) != nrow(Gamma_hat)) stop("gamma_hat length must equal rows of Gamma_hat")
  outcomes <- colnames(Gamma_hat)
  if (!is.null(outcomes) && anyDuplicated(outcomes)) stop("Gamma_hat column names must be unique outcome names")
  if (!is.null(rownames(R)) && !is.null(colnames(R))) {
    if (anyDuplicated(rownames(R)) || anyDuplicated(colnames(R))) {
      stop("R row and column names must be unique outcome names")
    }
    if (!setequal(rownames(R), colnames(R))) {
      stop("R row names and column names must contain the same outcomes")
    }
    R <- R[rownames(R), rownames(R), drop = FALSE]
    if (!is.null(outcomes)) {
      if (!setequal(outcomes, rownames(R))) {
        stop("Outcome names in Gamma_hat and R must match exactly")
      }
      R <- R[outcomes, outcomes, drop = FALSE]
    }
  }
  if (ncol(Gamma_hat) != nrow(R) || ncol(Gamma_hat) != ncol(R)) {
    stop("R dimensions must equal the number of outcome columns in Gamma_hat")
  }
  if (!all(is.finite(gamma_hat)) || !all(is.finite(Gamma_hat))) stop("Summary-statistic betas must be finite")
  if (!all(is.finite(R))) stop("R matrix must be finite")
  if (length(n1) != 1L || length(n2) != 1L) stop("n1 and n2 must be scalar sample sizes")
  if (!is.finite(n1) || !is.finite(n2) || n1 <= 0 || n2 <= 0) stop("n1 and n2 must be positive")
  if (max(abs(R - t(R))) > 1e-8) stop("R matrix must be symmetric")
  if (max(abs(diag(R) - 1)) > 1e-8) stop("R matrix diagonal must be 1")
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 0) {
    stop("R matrix must be positive definite; use estimate_outcome_correlation(..., near_pd = TRUE) or provide a positive-definite working correlation matrix")
  }
  list(gamma_hat = gamma_hat, Gamma_hat = Gamma_hat, R = R, n1 = as.integer(n1), n2 = as.integer(n2),
       min_eigenvalue = min(ev), positive_definite = all(ev > 0))
}

#' Estimate a working outcome correlation matrix
#'
#' Estimates cross-outcome correlation from a matrix of putatively null
#' outcome-GWAS z-scores. Rows are variants and columns are outcomes.
#'
#' @param null_z_scores Numeric matrix or data frame of null-variant z-scores.
#' @param near_pd If `TRUE`, project a non-positive-definite estimate to the
#'   nearest positive-definite correlation matrix using [Matrix::nearPD()].
#' @return A numeric correlation matrix with attributes recording whether
#'   nearPD projection was used and the final minimum eigenvalue.
#' @export
estimate_outcome_correlation <- function(null_z_scores, near_pd = TRUE) {
  z <- as.matrix(null_z_scores)
  storage.mode(z) <- "numeric"
  if (nrow(z) < 2L) stop("At least two null-variant rows are required to estimate a correlation matrix")
  if (ncol(z) < 1L) stop("At least one outcome column is required")
  if (is.null(colnames(z))) colnames(z) <- paste0("outcome_", seq_len(ncol(z)))
  R <- if (ncol(z) == 1L) {
    matrix(1, nrow = 1L, ncol = 1L, dimnames = list(colnames(z), colnames(z)))
  } else {
    stats::cor(z, use = "pairwise.complete.obs")
  }
  if (!all(is.finite(R))) stop("Estimated correlation matrix contains non-finite values; check missing or constant outcome z-score columns")
  diag(R) <- 1
  R <- (R + t(R)) / 2
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  near_pd_used <- FALSE
  if (near_pd && any(ev <= 0)) {
    R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
    near_pd_used <- TRUE
  }
  attr(R, "nearPD_used") <- near_pd_used
  attr(R, "min_eigenvalue") <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
  R
}

default_theta0 <- function(m, beta0 = 0) {
  c(rep(beta0, m), 0.001, rep(0.001, m), 1, rep(0.8, m))
}

validate_optimizer_controls <- function(maxiter, miniter, rd, tol) {
  maxiter <- suppressWarnings(as.numeric(maxiter))
  miniter <- suppressWarnings(as.numeric(miniter))
  rd <- suppressWarnings(as.numeric(rd))
  tol <- suppressWarnings(as.numeric(tol))
  if (length(maxiter) != 1L || !is.finite(maxiter) || maxiter < 1L) stop("maxiter must be a positive scalar")
  if (length(miniter) != 1L || !is.finite(miniter) || miniter < 0L) stop("miniter must be a non-negative scalar")
  if (miniter > maxiter) stop("miniter must be less than or equal to maxiter")
  if (length(rd) != 1L || !is.finite(rd) || rd <= 0) stop("rd must be a positive scalar")
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) stop("tol must be a positive scalar")
  list(maxiter = as.integer(maxiter), miniter = as.integer(miniter),
       rd = as.numeric(rd), tol = as.numeric(tol))
}

#' Fit the MR-MOSS model
#'
#' Fits the likelihood-based multi-outcome MR-MOSS model to harmonized summary
#' statistics and a working outcome correlation matrix.
#'
#' @param summary_stats Either a path readable by [read_mrmoss_summary_stats()]
#'   or the list returned by that function. If supplied, `gamma_hat`,
#'   `Gamma_hat`, `n1`, and `n2` are filled from it unless explicitly supplied.
#' @param gamma_hat Numeric vector of SNP-exposure association estimates.
#' @param Gamma_hat Numeric matrix of SNP-outcome association estimates.
#' @param R Numeric correlation matrix or path readable by
#'   [read_outcome_correlation()].
#' @param n1 Exposure GWAS sample size.
#' @param n2 Outcome GWAS sample size.
#' @param theta0 Optional starting vector for the C++ optimizer. Must have
#'   length `3 * number_of_outcomes + 2` with positive variance and scale
#'   components.
#' @param test Optional fitted outcome names or one-based outcome indices for
#'   outcome-specific LRTs.
#' @param maxiter Maximum number of C++ core iterations.
#' @param rd Numeric parameter passed to the C++ core.
#' @param tol Convergence tolerance passed to the C++ core.
#' @param miniter Minimum number of C++ core iterations.
#' @return An object of class `mrmoss_fit`.
#' @export
fit_mrmoss <- function(summary_stats = NULL, gamma_hat = NULL, Gamma_hat = NULL, R,
                       n1 = NULL, n2 = NULL, theta0 = NULL, test = NULL,
                       maxiter = 1000L, rd = 1, tol = 1e-5, miniter = 20L) {
  if (!is.null(summary_stats)) {
    parsed <- if (is.character(summary_stats)) read_mrmoss_summary_stats(summary_stats) else summary_stats
    gamma_hat <- gamma_hat %||% parsed$gamma_hat
    Gamma_hat <- Gamma_hat %||% parsed$Gamma_hat
    n1 <- n1 %||% parsed$n1
    n2 <- n2 %||% parsed$n2
  }
  if (is.character(R)) R <- read_outcome_correlation(R)
  checked <- check_mrmoss_inputs(gamma_hat, Gamma_hat, R, n1, n2)
  m <- ncol(checked$Gamma_hat)
  if (is.null(theta0)) theta0 <- default_theta0(m)
  if (length(theta0) != 3L * m + 2L) stop("theta0 must have length 3 * number_of_outcomes + 2")
  theta0 <- as.numeric(theta0)
  if (!all(is.finite(theta0))) stop("theta0 must contain only finite values")
  variance_idx <- (m + 1L):(3L * m + 2L)
  if (any(theta0[variance_idx] <= 0)) stop("theta0 variance and scale components must be positive")
  outcomes <- colnames(checked$Gamma_hat)
  if (is.null(outcomes)) outcomes <- paste0("outcome_", seq_len(m))
  if (is.null(test)) {
    test <- seq_len(m)
  } else if (is.character(test)) {
    test <- match(test, outcomes)
  } else {
    test <- as.integer(test)
  }
  if (!length(test) || any(is.na(test)) || any(test < 1L | test > m)) {
    stop("test must contain fitted outcome names or one-based outcome indices")
  }
  if (anyDuplicated(test)) stop("test outcomes must be unique")
  controls <- validate_optimizer_controls(maxiter, miniter, rd, tol)
  res <- MRMOSS_PX_xprod_lrt_diag_cpp(
    checked$gamma_hat, checked$Gamma_hat, checked$R,
    n1 = checked$n1, n2 = checked$n2, theta0 = theta0,
    test = as.integer(test), maxiter = controls$maxiter, rd = controls$rd,
    tol = controls$tol, miniter = controls$miniter
  )
  fit <- list(
    beta = as.numeric(res$beta),
    theta = as.numeric(res$theta),
    lrt = pmax(as.numeric(res$LRT), 0),
    pvalue = pmin(pmax(as.numeric(res$pvalue), 0), 1),
    log_pvalue = if (!is.null(res$log_pvalue)) as.numeric(res$log_pvalue) else rep(NA_real_, m),
    lrt_overall = max(as.numeric(res$LRT_overall), 0),
    pvalue_overall = pmin(pmax(as.numeric(res$pvalue_overall), 0), 1),
    log_pvalue_overall = if (!is.null(res$log_pvalue_overall)) as.numeric(res$log_pvalue_overall) else NA_real_,
    iteration = as.integer(res$iteration),
    raw = res,
    inputs = checked,
    outcomes = outcomes
  )
  class(fit) <- "mrmoss_fit"
  fit
}

#' Outcome-specific likelihood-ratio tests
#'
#' Extracts outcome-specific MR-MOSS likelihood-ratio statistics, p-values, and
#' capped `-log10(P)` summaries from a fitted model.
#'
#' @param fit An object returned by [fit_mrmoss()].
#' @param cap Maximum value used when reporting capped `-log10(P)`.
#' @return A data frame with one row per outcome.
#' @export
outcome_lrt <- function(fit, cap = 300) {
  if (!inherits(fit, "mrmoss_fit")) stop("fit must be returned by fit_mrmoss()")
  lp <- log10_p_value(fit$pvalue, log_p = fit$log_pvalue, cap = cap)
  data.frame(
    outcome = fit$outcomes,
    beta_hat = fit$beta,
    lrt = fit$lrt,
    p_value = fit$pvalue,
    minus_log10_p = lp$minus_log10_p,
    finite_p = lp$finite_p,
    censored = lp$censored,
    cap_value = cap,
    stringsAsFactors = FALSE
  )
}

#' All-outcome global likelihood-ratio test
#'
#' Extracts the MR-MOSS all-outcome global likelihood-ratio statistic and
#' p-value from a fitted model.
#'
#' @param fit An object returned by [fit_mrmoss()].
#' @param cap Maximum value used when reporting capped `-log10(P)`.
#' @return A one-row data frame.
#' @export
global_lrt <- function(fit, cap = 300) {
  if (!inherits(fit, "mrmoss_fit")) stop("fit must be returned by fit_mrmoss()")
  lp <- log10_p_value(fit$pvalue_overall, log_p = fit$log_pvalue_overall, cap = cap)
  data.frame(
    test = "all_outcome_global",
    df = length(fit$outcomes),
    lrt = fit$lrt_overall,
    p_value = fit$pvalue_overall,
    minus_log10_p = lp$minus_log10_p,
    finite_p = lp$finite_p,
    censored = lp$censored,
    cap_value = cap,
    stringsAsFactors = FALSE
  )
}

#' Subset likelihood-ratio tests
#'
#' Runs MR-MOSS likelihood-ratio tests for user-defined subsets of outcomes.
#'
#' @param fit An object returned by [fit_mrmoss()].
#' @param subsets A list of outcome names or one-based outcome indices.
#' @param subset_names Optional names for `subsets`.
#' @param maxiter Maximum number of C++ core iterations.
#' @param rd Numeric parameter passed to the C++ core.
#' @param tol Convergence tolerance passed to the C++ core.
#' @param miniter Minimum number of C++ core iterations.
#' @param cap Maximum value used when reporting capped `-log10(P)`.
#' @return A data frame with one row per subset.
#' @export
subset_lrt <- function(fit, subsets, subset_names = names(subsets), maxiter = 1000L,
                       rd = 1, tol = 1e-5, miniter = 20L, cap = 300) {
  if (!inherits(fit, "mrmoss_fit")) stop("fit must be returned by fit_mrmoss()")
  if (!is.list(subsets) || !length(subsets)) stop("subsets must be a non-empty list")
  if (is.null(subset_names) || any(!nzchar(subset_names))) subset_names <- paste0("subset_", seq_along(subsets))
  if (length(subset_names) != length(subsets)) stop("subset_names must have the same length as subsets")
  idx <- lapply(subsets, function(x) {
    if (is.character(x)) match(x, fit$outcomes) else as.integer(x)
  })
  if (any(vapply(idx, function(x) !length(x) || any(is.na(x)), logical(1)))) {
    stop("All subset outcomes must match fitted outcome names")
  }
  if (any(vapply(idx, function(x) any(x < 1L | x > length(fit$outcomes)), logical(1)))) {
    stop("Subset outcome indices must be between 1 and the number of fitted outcomes")
  }
  if (any(vapply(idx, anyDuplicated, integer(1)) > 0L)) {
    stop("Subset outcomes must be unique within each subset")
  }
  controls <- validate_optimizer_controls(maxiter, miniter, rd, tol)
  res <- MRMOSS_PX_xprod_subset_lrt_cpp(
    fit$inputs$gamma_hat, fit$inputs$Gamma_hat, fit$inputs$R,
    n1 = fit$inputs$n1, n2 = fit$inputs$n2, theta0 = fit$theta,
    subsets = idx, subset_names = as.character(subset_names),
    maxiter = controls$maxiter, rd = controls$rd,
    tol = controls$tol, miniter = controls$miniter
  )
  p <- pmin(pmax(as.numeric(res$pvalue_subset), 0), 1)
  lp <- log10_p_value(p, log_p = as.numeric(res$log_pvalue_subset), cap = cap)
  data.frame(
    subset = as.character(res$subset_names),
    df = as.integer(res$df_subset),
    lrt = pmax(as.numeric(res$LRT_subset), 0),
    p_value = p,
    minus_log10_p = lp$minus_log10_p,
    finite_p = lp$finite_p,
    censored = lp$censored,
    cap_value = cap,
    stringsAsFactors = FALSE
  )
}

#' Domain likelihood-ratio tests
#'
#' Runs subset LRTs from a domain map with columns `outcome` and `domain`.
#'
#' @param fit An object returned by [fit_mrmoss()].
#' @param domain_map A data frame or CSV path with columns `outcome` and
#'   `domain`.
#' @param ... Additional arguments passed to [subset_lrt()].
#' @return A data frame with one row per domain.
#' @export
domain_lrt <- function(fit, domain_map, ...) {
  if (is.character(domain_map)) domain_map <- mrmoss_read_table(domain_map)
  required <- c("outcome", "domain")
  missing <- setdiff(required, names(domain_map))
  if (length(missing)) stop("domain_map missing columns: ", paste(missing, collapse = ", "))
  subsets <- split(domain_map$outcome, domain_map$domain)
  subset_lrt(fit, subsets, subset_names = names(subsets), ...)
}

#' Check minimal MR-MOSS input columns
#'
#' Checks the minimal wide-format columns expected by the package. This helper
#' does not align alleles, flip beta signs, resolve palindromic variants, or
#' replace allele harmonization against raw GWAS files.
#'
#' @param summary_stats A data frame or CSV path.
#' @return The input data frame after checking required columns.
#' @export
harmonize_inputs <- function(summary_stats) {
  dat <- if (is.character(summary_stats)) mrmoss_read_table(summary_stats) else summary_stats
  required <- c("SNP", "effect_allele", "other_allele", "beta_exposure", "se_exposure")
  missing <- setdiff(required, names(dat))
  if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))
  dat
}

#' Numerically stable `-log10(P)` summaries
#'
#' Converts p-values, and optionally same-length natural-log p-values, to raw
#' and capped `-log10(P)` summaries for plotting and reporting very small
#' p-values.
#'
#' @param p Numeric p-values.
#' @param log_p Optional natural-log p-values.
#' @param cap Maximum reported value for capped `-log10(P)`.
#' @return A data frame with raw/capped values and flags for finite/censored
#'   entries.
#' @export
log10_p_value <- function(p, log_p = NULL, cap = 300) {
  p <- suppressWarnings(as.numeric(p))
  if (is.null(log_p)) {
    log_p <- rep(NA_real_, length(p))
    ok <- is.finite(p) & p > 0
    log_p[ok] <- log(p[ok])
    log_p[p == 0] <- -Inf
  } else {
    log_p <- suppressWarnings(as.numeric(log_p))
    if (length(log_p) != length(p)) stop("log_p must have the same length as p")
  }
  minus <- -as.numeric(log_p) / log(10)
  raw_minus <- -log10(p)
  use_raw <- is.finite(raw_minus)
  minus[use_raw] <- raw_minus[use_raw]
  data.frame(
    minus_log10_p = pmin(minus, cap),
    raw_minus_log10_p = minus,
    finite_p = is.finite(p) & p >= 0 & p <= 1,
    censored = !is.finite(minus) | minus > cap,
    stringsAsFactors = FALSE
  )
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x

mrmoss_read_table <- function(path) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("path must be a non-empty scalar character path")
  }
  if (!file.exists(path)) stop("File does not exist: ", path)
  lower <- tolower(path)
  is_tab <- grepl("\\.(tsv|tab|txt|sumstats)(\\.gz)?$", lower)
  sep <- if (is_tab) "\t" else ","
  con <- if (grepl("\\.gz$", lower)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  utils::read.table(
    con,
    header = TRUE,
    sep = sep,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
