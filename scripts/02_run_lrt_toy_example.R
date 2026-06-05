library(MRMOSS)

summary_path <- system.file("extdata", "toy_summary_stats.csv", package = "MRMOSS")
cor_path <- system.file("extdata", "toy_outcome_correlation.csv", package = "MRMOSS")
domain_path <- system.file("extdata", "toy_domain_map.csv", package = "MRMOSS")

toy <- read_mrmoss_summary_stats(summary_path)
R <- read_outcome_correlation(cor_path)
fit <- fit_mrmoss(summary_stats = toy, R = R, maxiter = 500)

print(outcome_lrt(fit))
print(global_lrt(fit))
print(domain_lrt(fit, domain_path))
print(subset_lrt(fit, list(first_two = c("outcome_1", "outcome_2"))))
