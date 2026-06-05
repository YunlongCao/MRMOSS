\name{list_mrmoss_examples}
\alias{list_mrmoss_examples}
\title{List installed MR-MOSS example datasets}
\usage{
list_mrmoss_examples(analysis = NULL, max_rows = NULL)
}
\arguments{
\item{analysis}{Optional analysis filter. Common values are
\code{"cvd_positive_control"}, \code{"mvp72"} and \code{"negative_control"}.}
\item{max_rows}{Optional maximum number of rows to return after filtering.}
}
\value{
A data frame with one row per installed example.
}
\description{
MR-MOSS example datasets are analysis-ready, processed summary-statistic
examples stored under \code{inst/extdata/examples/<example_id>/}. Each example
directory should contain \code{summary_stats.tsv}, \code{outcome_correlation.tsv},
optional \code{domain_map.tsv}, and optional \code{metadata.tsv}.
}
