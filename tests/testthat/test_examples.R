test_that("example registry is well formed even when no processed example is installed", {
  examples <- list_mrmoss_examples()
  expect_true(all(c("example_id", "title", "description") %in% names(examples)))
  expect_true(all(c("analysis", "exposure", "n_instruments", "n_outcomes") %in% names(examples)))
})

test_that("processed CVD example can be inspected and run", {
  examples <- list_mrmoss_examples()
  id <- "cvd__Apolipoprotein_B_levels__p5e_08"
  expect_true(id %in% examples$example_id)
  loaded <- load_mrmoss_example(id)
  expect_s3_class(loaded, "mrmoss_example_data")
  expect_equal(ncol(loaded$summary_stats$Gamma_hat), 5)
  info <- show_mrmoss_example(id)
  expect_s3_class(info, "mrmoss_example_info")
  expect_equal(info$n_instruments, 380)
  expect_equal(length(info$outcomes), 5)
  console_result <- run_mrmoss_example(id, open_report = FALSE, verbose = FALSE,
                                       maxiter = 500)
  expect_s3_class(console_result, "mrmoss_example_result")
  expect_equal(length(console_result$paths), 0)
  tmp <- tempfile("mrmoss-example-")
  result <- run_mrmoss_example(id, out_dir = tmp,
                               open_report = FALSE, verbose = FALSE,
                               maxiter = 500)
  expect_s3_class(result, "mrmoss_example_result")
  expect_true(file.exists(file.path(tmp, "mrmoss_report.html")))
  expect_equal(nrow(result$outcome_lrt), 5)
})

test_that("delimited readers support tab-delimited files", {
  tmp_summary <- tempfile(fileext = ".tsv")
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
    check.names = FALSE
  )
  utils::write.table(dat, tmp_summary, sep = "\t", row.names = FALSE, quote = FALSE)
  parsed <- read_mrmoss_summary_stats(tmp_summary)
  expect_equal(parsed$outcomes, "CAD")

  tmp_R <- tempfile(fileext = ".tsv")
  Rdat <- data.frame(outcome = "CAD", CAD = 1, check.names = FALSE)
  utils::write.table(Rdat, tmp_R, sep = "\t", row.names = FALSE, quote = FALSE)
  R <- read_outcome_correlation(tmp_R)
  expect_equal(rownames(R), "CAD")
  expect_equal(colnames(R), "CAD")
})
