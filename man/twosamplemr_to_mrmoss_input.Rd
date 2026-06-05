\name{twosamplemr_to_mrmoss_input}
\alias{twosamplemr_to_mrmoss_input}
\title{Convert harmonized TwoSampleMR output to MR-MOSS wide input}
\usage{
twosamplemr_to_mrmoss_input(harmonised_dat, outcome_names = NULL,
  output_path = NULL)
}
\arguments{
\item{harmonised_dat}{A data frame returned by \code{TwoSampleMR::harmonise_data()}
for one exposure and multiple outcomes.}

\item{outcome_names}{Optional character vector giving MR-MOSS outcome labels. If
named, names should match \code{id.outcome} values. If unnamed, values are used
in the order of outcomes in \code{harmonised_dat}.}

\item{output_path}{Optional path to write the wide input file. \code{.csv},
\code{.tsv}, \code{.txt} and \code{.gz} paths are supported.}
}
\value{
A data frame in MR-MOSS wide input format. The attribute \code{outcome_map}
records the mapping from TwoSampleMR outcome IDs to MR-MOSS outcome labels.
}
\description{
Converts the data frame returned by \code{TwoSampleMR::harmonise_data()} into
the wide summary-statistic format expected by \code{\link{read_mrmoss_summary_stats}}.
This helper does not call TwoSampleMR or OpenGWAS directly; it only reshapes an
already harmonized data frame.
}
