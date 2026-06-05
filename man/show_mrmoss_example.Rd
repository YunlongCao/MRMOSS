\name{show_mrmoss_example}
\alias{show_mrmoss_example}
\title{Show details for an installed MR-MOSS example dataset}
\usage{
show_mrmoss_example(example_id = NULL)
}
\arguments{
\item{example_id}{Installed example identifier. If omitted and exactly one
example is installed, that example is used.}
}
\value{
Invisibly returns an object of class \code{mrmoss_example_info}.
}
\description{
Prints the exposure, outcomes, input files, dimensions and core MR-MOSS
functions used by an installed processed example. This is a read-only helper
intended to orient new users before fitting.
}
