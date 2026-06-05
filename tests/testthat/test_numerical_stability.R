test_that("-log10 p-value capping tracks censored values", {
  out <- log10_p_value(c(1, 1e-10, 0), cap = 20)
  expect_equal(out$minus_log10_p[1], 0)
  expect_true(out$minus_log10_p[3] == 20)
  expect_true(out$censored[3])
  from_log <- log10_p_value(0, log_p = -1000, cap = 500)
  expect_gt(from_log$raw_minus_log10_p, 400)
  expect_false(from_log$censored)
  expect_error(log10_p_value(c(0.1, 0.2), log_p = -1), "same length")
})

test_that("citation metadata is available", {
  cit <- utils::citation("MRMOSS")
  expect_true(length(cit) >= 1)
})
