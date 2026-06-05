library(MRMOSS)

summary_path <- system.file("extdata", "toy_summary_stats.csv", package = "MRMOSS")
cor_path <- system.file("extdata", "toy_outcome_correlation.csv", package = "MRMOSS")

toy <- read_mrmoss_summary_stats(summary_path)
R <- read_outcome_correlation(cor_path)
checked <- check_mrmoss_inputs(toy$gamma_hat, toy$Gamma_hat, R, toy$n1, toy$n2)
print(checked[c("positive_definite", "min_eigenvalue")])

fit <- fit_mrmoss(summary_stats = toy, R = R, maxiter = 500)
print(data.frame(outcome = fit$outcomes, beta_hat = fit$beta))
