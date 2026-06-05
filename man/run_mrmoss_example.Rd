\name{run_mrmoss_example}
\alias{run_mrmoss_example}
\title{Run an installed MR-MOSS example dataset}
\usage{
run_mrmoss_example(example_id = NULL, out_dir = NULL, maxiter = 1000L,
  open_report = FALSE, verbose = TRUE, output_format = c("tsv", "csv"))
}
\arguments{
\item{example_id}{Installed example identifier. If omitted and exactly one
example is installed, that example is used. Use \code{list_mrmoss_examples()} to
see available examples.}
\item{out_dir}{Optional directory where result tables and the HTML report are
written. The default \code{NULL} keeps the first-run workflow in the R console
and writes no files.}
\item{maxiter}{Maximum number of C++ core iterations passed to \code{fit_mrmoss()}.}
\item{open_report}{If \code{TRUE}, open the generated HTML report. The default is
\code{FALSE} so that the first-run workflow stays in the console.}
\item{verbose}{If \code{TRUE}, print a concise console summary of the input data,
core functions used and key results.}
\item{output_format}{Result-table format, either \code{"tsv"} or \code{"csv"}.}
}
\value{
Invisibly returns an object of class \code{mrmoss_example_result}.
}
\description{
Runs a processed example dataset shipped with the package. This is the
recommended first-run workflow once a manuscript-style example has been added
under \code{inst/extdata/examples/}.
}
