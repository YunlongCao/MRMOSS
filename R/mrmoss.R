#' MRMOSS: Joint Mendelian Randomization of Correlated Outcomes
#'
#' Fits the MR-MOSS model to GWAS summary statistics for one exposure and
#' several correlated outcomes, and tests the causal effects with likelihood
#' ratio tests.
#'
#' @keywords internal
#' @useDynLib MRMOSS, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"

#' Read summary statistics
#'
#' Reads a tab-delimited file, which may be gzip-compressed, with one row per
#' instrument and the columns `SNP`, `effect_allele`, `other_allele`,
#' `exposure_sample_size`, `outcome_sample_size`, `beta_exposure` and one
#' column `beta_outcome_<name>` for each outcome. The associations are aligned
#' to the effect allele and standardized as the z-score divided by the square
#' root of the GWAS sample size.
#'
#' @param path Path to the file.
#' @return A list with the SNP identifiers `snp`, the SNP-exposure
#'   associations `gamma_hat`, the matrix of SNP-outcome associations
#'   `Gamma_hat` (instruments in rows, outcomes in columns) and the exposure
#'   and outcome sample sizes `n1` and `n2`.
#' @export
read_mrmoss_summary_stats <- function(path) {
  dat <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("SNP", "effect_allele", "other_allele", "exposure_sample_size",
                "outcome_sample_size", "beta_exposure")
  missing <- setdiff(required, names(dat))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
  outcome_cols <- grep("^beta_outcome_", names(dat), value = TRUE)
  if (!length(outcome_cols)) stop("No beta_outcome_<name> columns found")
  Gamma_hat <- as.matrix(dat[, outcome_cols, drop = FALSE])
  dimnames(Gamma_hat) <- list(dat$SNP, sub("^beta_outcome_", "", outcome_cols))
  list(snp = dat$SNP, gamma_hat = dat$beta_exposure, Gamma_hat = Gamma_hat,
       n1 = sample_size(dat$exposure_sample_size, "exposure_sample_size"),
       n2 = sample_size(dat$outcome_sample_size, "outcome_sample_size"))
}

sample_size <- function(x, name) {
  x <- unique(x)
  if (length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(name, " must be one positive value for all instruments")
  }
  x
}

#' Read an outcome correlation matrix
#'
#' Reads a tab-delimited file whose first column holds the outcome names and
#' whose remaining columns hold the correlation matrix.
#'
#' @param path Path to the file.
#' @return A numeric matrix with the outcome names as row and column names.
#' @export
read_outcome_correlation <- function(path) {
  dat <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  R <- as.matrix(dat[, -1, drop = FALSE])
  rownames(R) <- dat[[1]]
  R
}

#' Estimate the outcome correlation matrix
#'
#' Computes the correlation between outcomes from the z-scores of variants
#' that are not associated with any of the outcomes.
#'
#' @param z Matrix of z-scores with variants in rows and outcomes in columns.
#' @return The correlation matrix, replaced by the nearest positive definite
#'   correlation matrix if it is not positive definite.
#' @export
estimate_outcome_correlation <- function(z) {
  R <- stats::cor(as.matrix(z))
  if (min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
  }
  R
}

#' Fit MR-MOSS
#'
#' Estimates the causal effects of the exposure on all outcomes jointly by
#' maximum likelihood, using the PX-EM algorithm, and computes the likelihood
#' ratio tests for each outcome and for all outcomes together.
#'
#' @param data List returned by [read_mrmoss_summary_stats()], or a list with
#'   the elements `gamma_hat`, `Gamma_hat`, `n1` and `n2`.
#' @param R Outcome correlation matrix with the outcome names as row and
#'   column names.
#' @param maxiter Maximum number of PX-EM iterations.
#' @param tol Tolerance for the relative change in the log-likelihood.
#' @return An object of class `mrmoss_fit`.
#' @export
fit_mrmoss <- function(data, R, maxiter = 1e6, tol = 1e-12) {
  if (maxiter < 1 || tol <= 0) stop("maxiter and tol must be positive")
  Gamma_hat <- as.matrix(data$Gamma_hat)
  outcomes <- colnames(Gamma_hat)
  m <- ncol(Gamma_hat)
  if (length(data$gamma_hat) != nrow(Gamma_hat)) {
    stop("gamma_hat and Gamma_hat must have the same number of instruments")
  }
  if (!all(is.finite(data$gamma_hat)) || !all(is.finite(Gamma_hat))) {
    stop("Associations must be finite")
  }
  if (!setequal(rownames(R), outcomes) || !setequal(colnames(R), outcomes)) {
    stop("R must have the outcome names of Gamma_hat as row and column names")
  }
  R <- as.matrix(R)[outcomes, outcomes, drop = FALSE]
  if (max(abs(R - t(R))) > 1e-8 || max(abs(diag(R) - 1)) > 1e-8) {
    stop("R must be a correlation matrix")
  }
  if (min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop("R must be positive definite")
  }
  model <- list(gamma_hat = as.numeric(data$gamma_hat), Gamma_hat = Gamma_hat,
                R = R, n1 = data$n1, n2 = data$n2,
                theta0 = c(rep(0, m), 0.001, rep(0.001, m), 1, rep(0.8, m)),
                maxiter = as.integer(maxiter), tol = tol)
  full <- pxem(model, integer())
  null_loglik <- vapply(seq_len(m), function(i) pxem(model, i)$loglik, 0)
  global_loglik <- pxem(model, seq_len(m))$loglik
  lrt <- pmax(2 * (full$loglik - null_loglik), 0)
  lrt_all <- max(2 * (full$loglik - global_loglik), 0)
  fit <- list(outcomes = outcomes, beta = full$theta[seq_len(m)], theta = full$theta,
              loglik = full$loglik, iterations = full$iteration,
              lrt = lrt, pvalue = stats::pchisq(lrt, 1, lower.tail = FALSE),
              global_lrt = lrt_all,
              global_pvalue = stats::pchisq(lrt_all, m, lower.tail = FALSE),
              model = model)
  class(fit) <- "mrmoss_fit"
  fit
}

# Fits the model with the causal effects on the outcomes in `fixed` set to zero.
pxem <- function(model, fixed) {
  fit <- mrmoss_pxem(model$gamma_hat, model$Gamma_hat, model$R, model$n1, model$n2,
                     model$theta0, as.integer(fixed) - 1L, model$maxiter, model$tol)
  fit$theta <- as.numeric(fit$theta)
  fit
}

#' Outcome-specific tests
#'
#' @param fit An object returned by [fit_mrmoss()].
#' @return A data frame with the estimated causal effect, the likelihood ratio
#'   statistic and the P value for each outcome.
#' @export
outcome_lrt <- function(fit) {
  data.frame(outcome = fit$outcomes, beta = fit$beta, lrt = fit$lrt,
             p_value = fit$pvalue)
}

#' Global test
#'
#' Tests whether the causal effects on all outcomes are zero.
#'
#' @param fit An object returned by [fit_mrmoss()].
#' @return A data frame with the degrees of freedom, the likelihood ratio
#'   statistic and the P value.
#' @export
global_lrt <- function(fit) {
  data.frame(df = length(fit$outcomes), lrt = fit$global_lrt,
             p_value = fit$global_pvalue)
}

#' Tests for sets of outcomes
#'
#' Tests whether the causal effects on a set of outcomes are all zero, with
#' the effects on the other outcomes unrestricted.
#'
#' @param fit An object returned by [fit_mrmoss()].
#' @param subsets A named list of character vectors of outcome names.
#' @return A data frame with the degrees of freedom, the likelihood ratio
#'   statistic and the P value for each set.
#' @export
subset_lrt <- function(fit, subsets) {
  if (!is.list(subsets) || is.null(names(subsets))) stop("subsets must be a named list")
  index <- lapply(subsets, match, table = fit$outcomes)
  if (anyNA(unlist(index))) stop("Unknown outcome names in subsets")
  null_loglik <- vapply(index, function(i) pxem(fit$model, unique(i))$loglik, 0)
  df <- lengths(lapply(index, unique))
  lrt <- pmax(2 * (fit$loglik - null_loglik), 0)
  data.frame(subset = names(subsets), df = df, lrt = lrt,
             p_value = stats::pchisq(lrt, df, lower.tail = FALSE), row.names = NULL)
}
