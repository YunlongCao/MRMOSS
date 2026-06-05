#' List installed MR-MOSS example datasets
#'
#' MR-MOSS example datasets are analysis-ready, processed summary-statistic
#' examples stored under `inst/extdata/examples/<example_id>/`. Each example
#' directory should contain `summary_stats.tsv`, `outcome_correlation.tsv`,
#' optional `domain_map.tsv`, and optional `metadata.tsv`.
#'
#' @param analysis Optional analysis filter. Common values are
#'   `"cvd_positive_control"`, `"mvp72"` and `"negative_control"`.
#' @param max_rows Optional maximum number of rows to return after filtering.
#' @return A data frame with one row per installed example.
#' @export
list_mrmoss_examples <- function(analysis = NULL, max_rows = NULL) {
  root <- system.file("extdata", "examples", package = "MRMOSS")
  if (!nzchar(root) || !dir.exists(root)) {
    return(data.frame(
      example_id = character(),
      analysis = character(),
      panel_id = character(),
      exposure_id = character(),
      exposure = character(),
      iv_threshold = character(),
      n_instruments = integer(),
      n_outcomes = integer(),
      title = character(),
      description = character(),
      stringsAsFactors = FALSE
    ))
  }
  manifest_path <- file.path(root, "manifest.tsv")
  if (file.exists(manifest_path)) {
    out <- mrmoss_read_table(manifest_path)
    if (!is.null(analysis)) out <- out[out$analysis %in% analysis, , drop = FALSE]
    if (!is.null(max_rows)) out <- utils::head(out, max_rows)
    return(out)
  }
  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[dirs != root]
  if (!length(dirs)) {
    return(data.frame(
      example_id = character(),
      title = character(),
      description = character(),
      stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(dirs, function(d) {
    meta_path <- file.path(d, "metadata.tsv")
    if (file.exists(meta_path)) {
      meta <- mrmoss_read_table(meta_path)
      if (all(c("field", "value") %in% names(meta))) {
        vals <- stats::setNames(as.character(meta$value), as.character(meta$field))
        return(data.frame(
          example_id = basename(d),
          title = vals[["title"]] %||% basename(d),
          description = vals[["description"]] %||% "",
          stringsAsFactors = FALSE
        ))
      }
    }
    data.frame(
      example_id = basename(d),
      title = basename(d),
      description = "",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (!is.null(max_rows)) out <- utils::head(out, max_rows)
  out
}

#' Load an installed MR-MOSS example dataset
#'
#' Loads one processed example dataset and returns the analysis-ready MR-MOSS
#' input object. This is the clearest way to inspect the input contract: one
#' exposure vector, one SNP-by-outcome association matrix, a working
#' outcome-correlation matrix and an optional domain map.
#'
#' @param example_id Installed example identifier. Use
#'   [list_mrmoss_examples()] to see available choices.
#' @return A list of class `mrmoss_example_data` containing parsed inputs,
#'   metadata and source paths.
#' @export
load_mrmoss_example <- function(example_id = NULL) {
  resolved <- mrmoss_load_example(example_id)
  summary_stats <- read_mrmoss_summary_stats(resolved$summary_path)
  R <- read_outcome_correlation(resolved$cor_path)
  domain_map <- if (!is.na(resolved$domain_path)) mrmoss_read_table(resolved$domain_path) else NULL
  checked <- check_mrmoss_inputs(
    gamma_hat = summary_stats$gamma_hat,
    Gamma_hat = summary_stats$Gamma_hat,
    R = R,
    n1 = summary_stats$n1,
    n2 = summary_stats$n2
  )
  out <- list(
    example_id = resolved$example_id,
    metadata = resolved$metadata,
    summary_stats = summary_stats,
    R = checked$R,
    domain_map = domain_map,
    paths = list(
      summary_stats = resolved$summary_path,
      outcome_correlation = resolved$cor_path,
      domain_map = resolved$domain_path
    ),
    input_qc = data.frame(
      metric = c(
        "example_id",
        "n_instruments",
        "n_outcomes",
        "n1",
        "n2",
        "correlation_positive_definite",
        "correlation_min_eigenvalue"
      ),
      value = c(
        resolved$example_id,
        length(summary_stats$gamma_hat),
        ncol(summary_stats$Gamma_hat),
        summary_stats$n1,
        summary_stats$n2,
        checked$positive_definite,
        signif(checked$min_eigenvalue, 6)
      ),
      stringsAsFactors = FALSE
    )
  )
  class(out) <- "mrmoss_example_data"
  out
}

#' Show details for an installed MR-MOSS example dataset
#'
#' Prints the exposure, outcomes, input files, dimensions and core MR-MOSS
#' functions used by an installed processed example. This is a read-only helper
#' intended to orient new users before fitting.
#'
#' @param example_id Installed example identifier. If omitted and exactly one
#'   example is installed, that example is used.
#' @return Invisibly returns an object of class `mrmoss_example_info`.
#' @export
show_mrmoss_example <- function(example_id = NULL) {
  loaded <- load_mrmoss_example(example_id)
  info <- list(
    example_id = loaded$example_id,
    metadata = loaded$metadata,
    paths = loaded$paths,
    n_instruments = length(loaded$summary_stats$gamma_hat),
    outcomes = loaded$summary_stats$outcomes,
    n1 = loaded$summary_stats$n1,
    n2 = loaded$summary_stats$n2,
    min_eigenvalue = as.numeric(loaded$input_qc$value[loaded$input_qc$metric == "correlation_min_eigenvalue"]),
    positive_definite = loaded$input_qc$value[loaded$input_qc$metric == "correlation_positive_definite"],
    summary_columns = names(loaded$summary_stats$raw),
    domain_map = loaded$domain_map
  )
  class(info) <- "mrmoss_example_info"
  print(info)
  invisible(info)
}

#' Run an installed MR-MOSS example dataset
#'
#' Runs a processed example dataset shipped with the package. This is the
#' recommended first-run workflow once a manuscript-style example has been added
#' under `inst/extdata/examples/`.
#'
#' @param example_id Installed example identifier. If omitted and exactly one
#'   example is installed, that example is used. Use [list_mrmoss_examples()] to
#'   see available examples.
#' @param out_dir Optional directory where result tables and the HTML report are
#'   written. The default `NULL` keeps the first-run workflow in the R console
#'   and writes no files.
#' @param maxiter Maximum number of C++ core iterations passed to [fit_mrmoss()].
#' @param open_report If `TRUE`, open the generated HTML report. The default is
#'   `FALSE` so that the first-run workflow stays in the console.
#' @param verbose If `TRUE`, print a concise console summary of the input data,
#'   core functions used and key results.
#' @param output_format Result-table format, either `"tsv"` or `"csv"`.
#' @return Invisibly returns an object of class `mrmoss_example_result`.
#' @export
run_mrmoss_example <- function(example_id = NULL,
                               out_dir = NULL,
                               maxiter = 1000L,
                               open_report = FALSE,
                               verbose = TRUE,
                               output_format = c("tsv", "csv")) {
  output_format <- match.arg(output_format)
  loaded <- load_mrmoss_example(example_id)
  example_id <- loaded$example_id
  summary_stats <- loaded$summary_stats
  R <- loaded$R
  domain_map <- loaded$domain_map
  metadata <- loaded$metadata

  fit <- fit_mrmoss(summary_stats = summary_stats, R = R, maxiter = maxiter)
  outcome_results <- outcome_lrt(fit)
  global_results <- global_lrt(fit)
  domain_results <- if (!is.null(domain_map)) domain_lrt(fit, domain_map) else NULL

  title <- metadata[["title"]] %||% paste("MR-MOSS example:", example_id)
  paths <- list()
  if (isTRUE(open_report) && is.null(out_dir)) out_dir <- file.path("mrmoss_examples", example_id)
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    ext <- output_format
    paths <- list(
      input_qc = file.path(out_dir, paste0("input_qc.", ext)),
      outcome_lrt = file.path(out_dir, paste0("outcome_lrt.", ext)),
      global_lrt = file.path(out_dir, paste0("global_lrt.", ext)),
      fit = file.path(out_dir, "mrmoss_fit.rds"),
      report = file.path(out_dir, "mrmoss_report.html")
    )
    if (!is.null(domain_results)) paths$domain_lrt <- file.path(out_dir, paste0("domain_lrt.", ext))

    mrmoss_write_table(loaded$input_qc, paths$input_qc, output_format)
    mrmoss_write_table(outcome_results, paths$outcome_lrt, output_format)
    mrmoss_write_table(global_results, paths$global_lrt, output_format)
    if (!is.null(domain_results)) mrmoss_write_table(domain_results, paths$domain_lrt, output_format)
    saveRDS(fit, paths$fit)
    write_mrmoss_report(
      report_path = paths$report,
      input_qc = loaded$input_qc,
      global_results = global_results,
      outcome_results = outcome_results,
      domain_results = domain_results,
      title = title,
      note = metadata[["description"]] %||% NULL
    )
  }

  if (isTRUE(open_report) && !is.null(paths$report)) utils::browseURL(normalizePath(paths$report))
  result <- list(
    example_id = example_id,
    title = title,
    metadata = metadata,
    fit = fit,
    input_qc = loaded$input_qc,
    outcome_lrt = outcome_results,
    domain_lrt = domain_results,
    global_lrt = global_results,
    paths = paths,
    source_paths = loaded$paths,
    outcomes = summary_stats$outcomes
  )
  class(result) <- "mrmoss_example_result"
  if (isTRUE(verbose)) print(result)
  invisible(result)
}

#' Write a lightweight MR-MOSS HTML report
#'
#' Writes a static HTML report from MR-MOSS result tables. This helper is used
#' by [run_mrmoss_example()] and can also be used by analysis scripts that want
#' a simple local report without adding an R Markdown dependency.
#'
#' @param report_path Output HTML path.
#' @param input_qc Data frame with input QC summaries.
#' @param global_results Data frame returned by [global_lrt()].
#' @param outcome_results Data frame returned by [outcome_lrt()].
#' @param domain_results Optional data frame returned by [domain_lrt()].
#' @param title Report title.
#' @param note Optional short note shown below the title.
#' @return Invisibly returns `report_path`.
#' @export
write_mrmoss_report <- function(report_path, input_qc, global_results,
                                outcome_results, domain_results = NULL,
                                title = "MR-MOSS report", note = NULL) {
  if (!is.character(report_path) || length(report_path) != 1L || !nzchar(report_path)) {
    stop("report_path must be a non-empty scalar character path")
  }
  report_dir <- dirname(report_path)
  if (!dir.exists(report_dir)) dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

  sections <- c(
    "<h2>Input QC</h2>",
    mrmoss_table_to_html(input_qc),
    "<h2>All-outcome global LRT</h2>",
    mrmoss_table_to_html(global_results),
    "<h2>Outcome-specific LRTs</h2>",
    mrmoss_table_to_html(outcome_results)
  )
  if (!is.null(domain_results)) {
    sections <- c(sections, "<h2>Domain or subset LRTs</h2>", mrmoss_table_to_html(domain_results))
  }
  note_html <- if (is.null(note) || !nzchar(note)) {
    character(0)
  } else {
    paste0("<p class=\"note\">", mrmoss_escape_html(note), "</p>")
  }
  html <- c(
    "<!doctype html>",
    "<html>",
    "<head>",
    "<meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    paste0("<title>", mrmoss_escape_html(title), "</title>"),
    "<style>",
    "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;max-width:980px;margin:40px auto;line-height:1.45;color:#20242a;padding:0 18px}",
    "h1,h2{line-height:1.15}",
    "table{border-collapse:collapse;width:100%;margin:12px 0 28px}",
    "th{background:#17365d;color:white;text-align:left}",
    "th,td{border:1px solid #d6dbe1;padding:6px 8px;font-size:13px;vertical-align:top}",
    "code{background:#f4f6f8;padding:1px 4px;border-radius:3px}",
    ".note{background:#f7f9fc;border-left:4px solid #2f6fba;padding:10px 12px}",
    "</style>",
    "</head>",
    "<body>",
    paste0("<h1>", mrmoss_escape_html(title), "</h1>"),
    note_html,
    sections,
    "</body>",
    "</html>"
  )
  writeLines(html, report_path)
  invisible(report_path)
}

mrmoss_find_example_file <- function(example_dir, stem, required = TRUE) {
  candidates <- file.path(
    example_dir,
    paste0(stem, c(".tsv", ".tsv.gz", ".csv", ".csv.gz", ".txt", ".txt.gz"))
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    if (required) stop("Example is missing required file: ", stem, ".tsv or .csv")
    return(NA_character_)
  }
  hit[[1L]]
}

mrmoss_resolve_example_id <- function(example_id = NULL) {
  available <- list_mrmoss_examples()
  if (is.null(example_id)) {
    if (nrow(available) == 1L) return(available$example_id[[1L]])
    if (nrow(available) == 0L) {
      stop("No processed MR-MOSS examples are installed")
    }
    stop("Multiple examples are installed; choose one of: ",
         paste(available$example_id, collapse = ", "))
  }
  if (!is.character(example_id) || length(example_id) != 1L || !nzchar(example_id)) {
    stop("example_id must be a non-empty scalar character value")
  }
  if (!example_id %in% available$example_id) {
    msg <- if (nrow(available)) paste(available$example_id, collapse = ", ") else "none installed"
    stop("Example '", example_id, "' is not installed. Available examples: ", msg)
  }
  example_id
}

mrmoss_load_example <- function(example_id = NULL) {
  example_id <- mrmoss_resolve_example_id(example_id)
  root <- system.file("extdata", "examples", package = "MRMOSS")
  example_dir <- file.path(root, example_id)
  list(
    example_id = example_id,
    example_dir = example_dir,
    summary_path = mrmoss_find_example_file(example_dir, "summary_stats"),
    cor_path = mrmoss_find_example_file(example_dir, "outcome_correlation"),
    domain_path = mrmoss_find_example_file(example_dir, "domain_map", required = FALSE),
    metadata = mrmoss_read_example_metadata(example_dir)
  )
}

mrmoss_read_example_metadata <- function(example_dir) {
  path <- file.path(example_dir, "metadata.tsv")
  if (!file.exists(path)) return(character())
  meta <- mrmoss_read_table(path)
  if (!all(c("field", "value") %in% names(meta))) return(character())
  stats::setNames(as.character(meta$value), as.character(meta$field))
}

print.mrmoss_example_info <- function(x, ...) {
  cat("MR-MOSS processed example\n")
  cat("  example_id: ", x$example_id, "\n", sep = "")
  cat("  title:      ", x$metadata[["title"]] %||% x$example_id, "\n", sep = "")
  if (!is.null(x$metadata[["description"]])) {
    cat("  about:      ", x$metadata[["description"]], "\n", sep = "")
  }
  if (!is.null(x$metadata[["exposure"]])) {
    cat("  exposure:   ", x$metadata[["exposure"]], "\n", sep = "")
  }
  cat("  outcomes:   ", paste(x$outcomes, collapse = ", "), "\n", sep = "")
  cat("  instruments:", x$n_instruments, "\n")
  cat("  n1/n2:      ", x$n1, " / ", x$n2, "\n", sep = "")
  cat("  R min eigen:", signif(x$min_eigenvalue, 4), "\n")
  cat("\nInput files\n")
  cat("  summary stats:       ", normalizePath(x$paths$summary_stats), "\n", sep = "")
  cat("  outcome correlation: ", normalizePath(x$paths$outcome_correlation), "\n", sep = "")
  if (!is.na(x$paths$domain_map)) {
    cat("  domain map:          ", normalizePath(x$paths$domain_map), "\n", sep = "")
  }
  cat("\nSummary-statistic columns\n")
  cat("  ", paste(x$summary_columns, collapse = ", "), "\n", sep = "")
  cat("\nCore MR-MOSS calls used by this example\n")
  cat("  summary_stats <- read_mrmoss_summary_stats(summary_path)\n")
  cat("  R <- read_outcome_correlation(correlation_path)\n")
  cat("  check_mrmoss_inputs(summary_stats$gamma_hat, summary_stats$Gamma_hat, R, summary_stats$n1, summary_stats$n2)\n")
  cat("  fit <- fit_mrmoss(summary_stats = summary_stats, R = R)\n")
  cat("  outcome_lrt(fit); global_lrt(fit); domain_lrt(fit, domain_path)\n")
  invisible(x)
}

print.mrmoss_example_data <- function(x, ...) {
  cat("MR-MOSS example input\n")
  cat("  example_id: ", x$example_id, "\n", sep = "")
  cat("  exposure:   ", x$metadata[["exposure"]] %||% "", "\n", sep = "")
  cat("  analysis:   ", x$metadata[["analysis"]] %||% "", "\n", sep = "")
  cat("  instruments:", length(x$summary_stats$gamma_hat), "\n")
  cat("  outcomes:   ", length(x$summary_stats$outcomes), "\n", sep = "")
  shown <- utils::head(x$summary_stats$outcomes, 8)
  cat("  outcome labels: ", paste(shown, collapse = ", "),
      if (length(x$summary_stats$outcomes) > length(shown)) ", ..." else "",
      "\n", sep = "")
  cat("\nMatrices available\n")
  cat("  gamma_hat: ", length(x$summary_stats$gamma_hat), "\n", sep = "")
  cat("  Gamma_hat: ", nrow(x$summary_stats$Gamma_hat), " x ", ncol(x$summary_stats$Gamma_hat), "\n", sep = "")
  cat("  R:         ", nrow(x$R), " x ", ncol(x$R), "\n", sep = "")
  if (!is.null(x$domain_map)) {
    cat("  domain_map:", nrow(x$domain_map), " rows, ", length(unique(x$domain_map$domain)), " domains\n", sep = "")
  }
  cat("\nRun this example with:\n")
  cat("  res <- run_mrmoss_example(\"", x$example_id, "\")\n", sep = "")
  invisible(x)
}

print.mrmoss_example_result <- function(x, ...) {
  cat("MR-MOSS example result\n")
  cat("  example_id: ", x$example_id, "\n", sep = "")
  cat("  title:      ", x$title, "\n", sep = "")
  if (!is.null(x$metadata[["exposure"]])) {
    cat("  exposure:   ", x$metadata[["exposure"]], "\n", sep = "")
  }
  cat("  outcomes:   ", paste(x$outcomes, collapse = ", "), "\n", sep = "")
  cat("  instruments:", x$input_qc$value[x$input_qc$metric == "n_instruments"], "\n")
  cat("\nKey result\n")
  g <- x$global_lrt[1, , drop = FALSE]
  cat("  all-outcome global LRT: df=", g$df,
      ", -log10(P)=", signif(g$minus_log10_p, 4),
      ", P=", format(g$p_value, scientific = TRUE, digits = 4), "\n", sep = "")
  if (nrow(x$outcome_lrt)) {
    top <- x$outcome_lrt[order(-x$outcome_lrt$minus_log10_p), , drop = FALSE]
    top <- utils::head(top[, c("outcome", "beta_hat", "minus_log10_p", "p_value")], 5)
    cat("\nTop outcome-specific LRTs\n")
    print(top, row.names = FALSE)
  }
  if (!is.null(x$domain_lrt) && nrow(x$domain_lrt)) {
    topd <- x$domain_lrt[order(-x$domain_lrt$minus_log10_p), , drop = FALSE]
    topd <- utils::head(topd[, c("subset", "df", "minus_log10_p", "p_value")], 5)
    cat("\nTop domain/subset LRTs\n")
    print(topd, row.names = FALSE)
  }
  if (length(x$paths)) {
    cat("\nOutputs written\n")
    out_dir <- dirname(x$paths$report)
    cat("  directory: ", normalizePath(out_dir), "\n", sep = "")
    for (nm in names(x$paths)) {
      cat("  ", nm, ": ", normalizePath(x$paths[[nm]]), "\n", sep = "")
    }
    cat("\nOpen the report manually with:\n")
    cat("  browseURL(\"", normalizePath(x$paths$report), "\")\n", sep = "")
  } else {
    cat("\nNo files were written. Supply out_dir = \"...\" to export tables and an HTML report.\n")
  }
  invisible(x)
}

mrmoss_write_table <- function(x, path, format = c("tsv", "csv")) {
  format <- match.arg(format)
  if (format == "tsv") {
    utils::write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  } else {
    utils::write.csv(x, path, row.names = FALSE, quote = FALSE, na = "")
  }
}

mrmoss_escape_html <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

mrmoss_table_to_html <- function(x) {
  x <- as.data.frame(x)
  header <- paste(sprintf("<th>%s</th>", mrmoss_escape_html(names(x))), collapse = "")
  if (!nrow(x)) {
    body <- paste0("<tr><td colspan=\"", ncol(x), "\">No rows</td></tr>")
  } else {
    body <- paste(apply(x, 1, function(row) {
      paste0(
        "<tr>",
        paste(sprintf("<td>%s</td>", mrmoss_escape_html(row)), collapse = ""),
        "</tr>"
      )
    }), collapse = "\n")
  }
  paste0(
    "<table><thead><tr>", header, "</tr></thead><tbody>",
    body,
    "</tbody></table>"
  )
}
