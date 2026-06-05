test_that("toy fit completes and dimensions match", {
  toy <- read_mrmoss_summary_stats(system.file("extdata", "toy_summary_stats.csv", package = "MRMOSS"))
  R <- read_outcome_correlation(system.file("extdata", "toy_outcome_correlation.csv", package = "MRMOSS"))
  fit <- fit_mrmoss(summary_stats = toy, R = R, maxiter = 500)
  expect_s3_class(fit, "mrmoss_fit")
  expect_equal(length(fit$beta), ncol(toy$Gamma_hat))
  expect_equal(nrow(outcome_lrt(fit)), ncol(toy$Gamma_hat))
})

test_that("toy summary statistic file has expected public columns", {
  toy_path <- system.file("extdata", "toy_summary_stats.csv", package = "MRMOSS")
  dat <- utils::read.csv(toy_path, check.names = FALSE)
  expect_gte(nrow(dat), 50)
  expect_true(all(c("SNP", "effect_allele", "other_allele", "beta_exposure", "se_exposure") %in% names(dat)))
  expect_equal(length(grep("^beta_outcome_", names(dat))), 3)
  expect_equal(length(grep("^se_outcome_", names(dat))), 3)
  toy <- read_mrmoss_summary_stats(toy_path)
  expect_equal(dim(toy$se_Gamma), dim(toy$Gamma_hat))
  expect_true(all(toy$se_gamma > 0))
  expect_true(all(toy$se_Gamma > 0))
})

test_that("outcome names are parsed from documented beta column patterns", {
  tmp <- tempfile(fileext = ".csv")
  dat <- data.frame(
    SNP = c("rs1", "rs2"),
    effect_allele = c("A", "C"),
    other_allele = c("G", "T"),
    beta_exposure = c(0.1, 0.2),
    se_exposure = c(0.01, 0.02),
    exposure_sample_size = c(1000, 1000),
    outcome_sample_size = c(900, 900),
    beta_outcome_CAD = c(0.03, 0.04),
    se_outcome_CAD = c(0.01, 0.01),
    beta_outcome_1 = c(0.05, 0.06),
    se_outcome_1 = c(0.02, 0.02),
    check.names = FALSE
  )
  utils::write.csv(dat, tmp, row.names = FALSE)
  parsed <- read_mrmoss_summary_stats(tmp)
  expect_equal(parsed$outcomes, c("CAD", "outcome_1"))
})

test_that("fit_mrmoss validates user-controlled optimizer inputs", {
  toy <- read_mrmoss_summary_stats(system.file("extdata", "toy_summary_stats.csv", package = "MRMOSS"))
  R <- read_outcome_correlation(system.file("extdata", "toy_outcome_correlation.csv", package = "MRMOSS"))
  expect_error(fit_mrmoss(summary_stats = toy, R = R, theta0 = 1), "theta0")
  expect_error(fit_mrmoss(summary_stats = toy, R = R, test = 0), "test")
  expect_error(fit_mrmoss(summary_stats = toy, R = R, maxiter = 0), "maxiter")
  expect_error(fit_mrmoss(summary_stats = toy, R = R, miniter = 10, maxiter = 5), "miniter")
  expect_error(fit_mrmoss(summary_stats = toy, R = R, rd = 0), "rd")
  expect_error(fit_mrmoss(summary_stats = toy, R = R, tol = 0), "tol")
  expect_error(fit_mrmoss(summary_stats = toy, gamma_hat = toy$gamma_hat[-1], R = R),
               "gamma_hat length")
})
