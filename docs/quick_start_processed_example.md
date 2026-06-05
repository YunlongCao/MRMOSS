# Built-in processed examples

The primary quick-start is now in the repository `README.md`. MR-MOSS ships a
set of manuscript-derived processed examples rather than a single toy dataset.
Each example includes:

- a complete-case harmonized IV set;
- a wide SNP-by-outcome summary-statistic file;
- a working outcome-correlation matrix;
- an optional outcome-domain map.

List available examples from R:

```r
library(MRMOSS)

list_mrmoss_examples()
list_mrmoss_examples(analysis = "cvd_positive_control")
list_mrmoss_examples(analysis = "mvp72", max_rows = 10)
```

Load and inspect one example:

```r
example_id <- "cvd__Apolipoprotein_B_levels__p5e_08"
dat <- load_mrmoss_example(example_id)
dat

dim(dat$summary_stats$Gamma_hat)
dat$R
dat$domain_map
```

Run MR-MOSS and print concise results:

```r
res <- run_mrmoss_example(example_id)
res$global_lrt
head(res$outcome_lrt)
head(res$domain_lrt)
```
