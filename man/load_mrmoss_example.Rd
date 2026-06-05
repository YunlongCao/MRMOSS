\name{load_mrmoss_example}
\alias{load_mrmoss_example}
\title{Load an installed MR-MOSS example dataset}
\usage{
load_mrmoss_example(example_id = NULL)
}
\arguments{
\item{example_id}{Installed example identifier. Use
\code{list_mrmoss_examples()} to see available choices.}
}
\value{
A list of class \code{mrmoss_example_data} containing parsed inputs, metadata
and source paths.
}
\description{
Loads one processed example dataset and returns the analysis-ready MR-MOSS input
object. This is the clearest way to inspect the input contract: one exposure
vector, one SNP-by-outcome association matrix, a working outcome-correlation
matrix and an optional domain map.
}
