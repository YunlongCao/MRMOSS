\name{twosamplemr_to_outcome_correlation}
\alias{twosamplemr_to_outcome_correlation}
\title{Estimate an MR-MOSS outcome-correlation matrix from TwoSampleMR outcome data}
\usage{
twosamplemr_to_outcome_correlation(outcome_dat, outcome_names = NULL,
  pval_threshold = NULL, near_pd = TRUE, output_path = NULL)
}
\arguments{
\item{outcome_dat}{Data frame containing \code{SNP}, outcome identifier,
\code{beta.outcome} and \code{se.outcome} columns.}

\item{outcome_names}{Optional labels passed to the MR-MOSS correlation matrix.
The same mapping should be used in \code{\link{twosamplemr_to_mrmoss_input}}.}

\item{pval_threshold}{Optional threshold. If \code{pval.outcome} is present,
rows with \code{pval.outcome <= pval_threshold} are excluded before estimating
the matrix.}

\item{near_pd}{If \code{TRUE}, allow \code{\link{estimate_outcome_correlation}}
to project a non-positive-definite estimate to a nearest positive-definite
correlation matrix.}

\item{output_path}{Optional path to write the correlation matrix.}
}
\value{
A numeric outcome-correlation matrix.
}
\description{
Converts a data frame returned by \code{TwoSampleMR::extract_outcome_data()} for
a candidate SNP set into outcome z-score columns, filters on outcome association
strength when requested, then estimates the working outcome-correlation matrix
with \code{\link{estimate_outcome_correlation}}.
}
