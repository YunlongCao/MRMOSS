\name{mrmoss_null_snp_panel}
\alias{mrmoss_null_snp_panel}
\title{Return the built-in MR-MOSS candidate SNP panel for estimating outcome correlation}
\usage{
mrmoss_null_snp_panel(n = 5000L, exclude = NULL)
}
\arguments{
\item{n}{Number of SNPs to return. Use \code{Inf} to return all packaged SNPs.}

\item{exclude}{Optional vector of SNP IDs to remove, for example exposure
instruments.}
}
\value{
A character vector of rsIDs.
}
\description{
Returns a built-in panel of common rsID variants that can be queried with
\code{TwoSampleMR::extract_outcome_data()} and then filtered as putatively null
variants by \code{\link{twosamplemr_to_outcome_correlation}}. The panel is
intended as a convenient default for workflow setup and small analyses; users can
replace it with an analysis-specific null-variant set for final studies.
}
