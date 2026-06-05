\name{write_mrmoss_report}
\alias{write_mrmoss_report}
\title{Write a lightweight MR-MOSS HTML report}
\usage{
write_mrmoss_report(report_path, input_qc, global_results, outcome_results,
  domain_results = NULL, title = "MR-MOSS report", note = NULL)
}
\arguments{
\item{report_path}{Output HTML path.}
\item{input_qc}{Data frame with input QC summaries.}
\item{global_results}{Data frame returned by \code{global_lrt()}.}
\item{outcome_results}{Data frame returned by \code{outcome_lrt()}.}
\item{domain_results}{Optional data frame returned by \code{domain_lrt()}.}
\item{title}{Report title.}
\item{note}{Optional short note shown below the title.}
}
\value{
Invisibly returns \code{report_path}.
}
\description{
Writes a static HTML report from MR-MOSS result tables. This helper is used by
\code{run_mrmoss_example()} and can also be used by analysis scripts that want a
simple local report without adding an R Markdown dependency.
}
