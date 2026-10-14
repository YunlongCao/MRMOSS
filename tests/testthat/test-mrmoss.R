path <- system.file("extdata", "positive_control", package = "MRMOSS")
R <- read_outcome_correlation(file.path(path, "outcome_correlation.tsv"))

test_that("input files are read correctly", {
  data <- read_mrmoss_summary_stats(file.path(path, "LPA_5e-08.tsv.gz"))
  expect_equal(dim(data$Gamma_hat), c(42L, 9L))
  expect_equal(colnames(data$Gamma_hat), rownames(R))
  expect_equal(data$n1, 343821)
})

test_that("the Lp(a) panel reproduces the published global P value", {
  data <- read_mrmoss_summary_stats(file.path(path, "LPA_5e-08.tsv.gz"))
  published <- utils::read.delim(file.path(path, "published_global_p.tsv"),
                                 colClasses = c(threshold = "character"))
  p <- published$published_p[published$exposure == "LPA" & published$threshold == "5e-08"]
  fit <- fit_mrmoss(data, R)
  expect_equal(global_lrt(fit)$p_value, p, tolerance = 0.01)
  expect_equal(nrow(outcome_lrt(fit)), 9L)
  expect_equal(subset_lrt(fit, list(all = rownames(R)))$lrt, global_lrt(fit)$lrt,
               tolerance = 1e-6)
})

test_that("the correlation matrix is checked", {
  data <- read_mrmoss_summary_stats(file.path(path, "LPA_5e-08.tsv.gz"))
  bad <- R
  bad[1, 2] <- bad[2, 1] <- 2
  expect_error(fit_mrmoss(data, bad), "positive definite")
  expect_error(fit_mrmoss(data, R[-1, -1]), "outcome names")
})
