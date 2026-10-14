# Reproduces the MR-MOSS global P values of the positive control analysis
# (Supplementary Table 7).

library(MRMOSS)

path <- system.file("extdata", "positive_control", package = "MRMOSS")
R <- read_outcome_correlation(file.path(path, "outcome_correlation.tsv"))
results <- utils::read.delim(file.path(path, "published_global_p.tsv"),
                             colClasses = c(threshold = "character"))

results$reproduced_p <- NA_real_
for (i in seq_len(nrow(results))) {
  file <- sprintf("%s_%s.tsv.gz", results$exposure[i], results$threshold[i])
  message("Fitting ", file)
  fit <- fit_mrmoss(read_mrmoss_summary_stats(file.path(path, file)), R)
  results$reproduced_p[i] <- global_lrt(fit)$p_value
}

print(results, digits = 3, row.names = FALSE)
