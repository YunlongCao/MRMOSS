test_that("R matrix checks and nearPD handling work", {
  R <- read_outcome_correlation(system.file("extdata", "toy_outcome_correlation.csv", package = "MRMOSS"))
  expect_equal(max(abs(R - t(R))), 0)
  expect_equal(unname(diag(R)), rep(1, nrow(R)))
  expect_true(min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) > 0)
  expect_equal(rownames(R), colnames(R))
  bad <- matrix(c(1, 1.2, 1.2, 1), 2)
  fixed <- estimate_outcome_correlation(bad, near_pd = TRUE)
  expect_true(attr(fixed, "min_eigenvalue") > 0)
})

test_that("R matrix names are aligned to outcome columns", {
  Gamma <- matrix(0, nrow = 2, ncol = 3)
  colnames(Gamma) <- c("outcome_1", "outcome_2", "outcome_3")
  R <- matrix(
    c(1, 0.1, 0.2, 0.1, 1, 0.3, 0.2, 0.3, 1),
    nrow = 3,
    dimnames = list(c("outcome_3", "outcome_1", "outcome_2"),
                    c("outcome_3", "outcome_1", "outcome_2"))
  )
  checked <- check_mrmoss_inputs(rep(0.1, 2), Gamma, R, 100, 100)
  expect_equal(rownames(checked$R), colnames(Gamma))
  expect_error(check_mrmoss_inputs(rep(0.1, 2), Gamma, R[1:2, 1:2], 100, 100),
               "must match")
  non_pd <- matrix(c(1, 1.2, 1.2, 1), 2)
  expect_error(check_mrmoss_inputs(rep(0.1, 2), Gamma[, 1:2], non_pd, 100, 100),
               "positive definite")
})

test_that("toy config is valid YAML for examples", {
  skip_if_not_installed("yaml")
  cfg <- yaml::read_yaml(system.file("extdata", "toy_config.yml", package = "MRMOSS"))
  expect_true(all(c("summary_stats", "outcome_correlation", "domain_map", "maxiter", "rd") %in% names(cfg)))
  expect_equal(cfg$seed, 20260604)
  expect_equal(cfg$n_snps, 80)
})
