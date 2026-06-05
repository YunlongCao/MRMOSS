tempfile_input <- function(dat, ext) {
  path <- tempfile(fileext = ext)
  utils::write.table(dat, path, sep = "\t", row.names = FALSE, quote = FALSE)
  path
}

test_that("built-in null SNP panel is available and excludes instruments", {
  panel <- mrmoss_null_snp_panel(n = 10)
  expect_equal(length(panel), 10)
  expect_true(all(grepl("^rs[0-9]+$", panel)))
  reduced <- mrmoss_null_snp_panel(n = 10, exclude = panel[1:3])
  expect_equal(length(reduced), 10)
  expect_false(any(panel[1:3] %in% reduced))
  all_panel <- mrmoss_null_snp_panel(n = Inf)
  expect_true(length(all_panel) >= 10000)
})

test_that("TwoSampleMR harmonised output converts to MR-MOSS wide input", {
  dat <- data.frame(
    SNP = rep(c("rs1", "rs2", "rs3"), times = 2),
    id.outcome = rep(c("ieu-a-7", "ieu-a-26"), each = 3),
    outcome = rep(c("Coronary heart disease", "Type 2 diabetes"), each = 3),
    beta.exposure = rep(c(0.10, 0.20, 0.30), times = 2),
    se.exposure = rep(c(0.01, 0.02, 0.03), times = 2),
    beta.outcome = c(0.03, 0.04, 0.05, 0.01, 0.02, 0.03),
    se.outcome = rep(0.01, 6),
    effect_allele.exposure = rep(c("A", "C", "G"), times = 2),
    other_allele.exposure = rep(c("G", "T", "A"), times = 2),
    samplesize.exposure = 100000,
    samplesize.outcome = rep(c(90000, 80000), each = 3),
    mr_keep = TRUE,
    stringsAsFactors = FALSE
  )
  wide <- twosamplemr_to_mrmoss_input(
    dat,
    outcome_names = c("ieu-a-7" = "CHD", "ieu-a-26" = "T2D")
  )
  expect_equal(nrow(wide), 3)
  expect_true(all(c("beta_outcome_CHD", "se_outcome_CHD",
                    "beta_outcome_T2D", "se_outcome_T2D") %in% names(wide)))
  expect_equal(attr(wide, "outcome_map")$mrmoss_outcome, c("CHD", "T2D"))
  parsed <- read_mrmoss_summary_stats(tempfile_input(wide, ".tsv"))
  expect_equal(parsed$outcomes, c("CHD", "T2D"))
})

test_that("TwoSampleMR null outcome data converts to a correlation matrix", {
  dat <- data.frame(
    SNP = rep(paste0("rs", 1:5), times = 2),
    id.outcome = rep(c("ieu-a-7", "ieu-a-26"), each = 5),
    beta.outcome = c(0.01, 0.02, -0.01, 0.03, 0.00,
                     0.02, 0.01, -0.02, 0.01, 0.01),
    se.outcome = rep(0.01, 10),
    pval.outcome = rep(0.5, 10),
    stringsAsFactors = FALSE
  )
  R <- twosamplemr_to_outcome_correlation(
    dat,
    outcome_names = c("ieu-a-7" = "CHD", "ieu-a-26" = "T2D"),
    pval_threshold = 1e-5
  )
  expect_equal(rownames(R), c("CHD", "T2D"))
  expect_equal(colnames(R), c("CHD", "T2D"))
  expect_equal(diag(R), c(1, 1))
})
