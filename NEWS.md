# MRMOSS 1.0.0

- Manuscript submission release.
- Provides a compact public R package for fitting the likelihood-based MR-MOSS
  model to harmonized GWAS summary statistics and a working outcome correlation
  matrix.
- Adds user-facing helpers for reading toy input files, checking dimensions,
  fitting the joint model, reporting outcome-specific likelihood-ratio tests,
  reporting all-outcome global tests, and running subset/domain-level tests.
- Strengthens R-side validation for outcome-correlation matrices, required
  standard-error columns, scalar sample sizes, user-supplied optimizer controls,
  and subset/test outcome indices.
- Reports log p-values from the outcome-specific and global C++ LRT paths so
  `-log10(P)` summaries remain informative when p-values underflow.
- Adds a GitHub-ready quick start, vignette, toy data, example scripts,
  testthat tests, and citation metadata.
- Adds a processed manuscript-style example, a built-in candidate SNP panel for
  outcome-correlation estimation, and TwoSampleMR/OpenGWAS bridge helpers for
  converting harmonized online GWAS lookups into MR-MOSS inputs.
- Expands real-GWAS onboarding with copy-and-edit analysis and input-preparation
  templates, explicit allele harmonization guidance, scalar sample-size caveats,
  and a data-availability table for public release.
