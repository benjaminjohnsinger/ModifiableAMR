#!/usr/bin/env Rscript
# promote_metrics.R — promote passing pending metrics to canonical
#
# Usage:
#   Rscript scripts/promote_metrics.R [scenario]
#
# Reads manifests/expected_metrics.csv, evaluates all pending rows against
# current output files, and promotes rows that pass tolerance to canonical.
# Rows with missing files or failing values are left as pending with a note.
#
# Run this after completing make manuscript SCENARIO=main (and/or supplementary
# scenario runs) to lock numeric reproducibility values.

source("config.R")

args <- commandArgs(trailingOnly = TRUE)
scenario_filter <- if (length(args) >= 1 && args[[1]] != "") args[[1]] else NULL

metrics_path <- file.path(AMR_CONFIG$paths$manifests, "expected_metrics.csv")
metrics <- read.csv(metrics_path, stringsAsFactors = FALSE, check.names = FALSE)

pending <- metrics[metrics$status == "pending", , drop = FALSE]

if (!is.null(scenario_filter)) {
  pending <- pending[pending$scenario == scenario_filter, , drop = FALSE]
}

if (nrow(pending) == 0) {
  message("[promote-metrics] No pending metrics",
          if (!is.null(scenario_filter)) paste0(" for scenario: ", scenario_filter) else "",
          ". Nothing to promote.")
  quit(status = 0)
}

compute_metric <- function(df, agg, column) {
  if (agg == "nrow") return(nrow(df))
  if (!column %in% names(df)) stop("Column not found: ", column)
  x <- df[[column]]
  if (agg == "sum") return(sum(x, na.rm = TRUE))
  if (agg == "mean") return(mean(x, na.rm = TRUE))
  if (agg == "median") return(median(x, na.rm = TRUE))
  if (agg == "unique_n") return(length(unique(x)))
  stop("Unsupported aggregate: ", agg)
}

today <- format(Sys.Date(), "%Y-%m-%d")
promoted <- character(0)
skipped  <- character(0)

for (i in seq_len(nrow(pending))) {
  row <- pending[i, ]
  idx <- which(metrics$metric_id == row$metric_id)

  if (!file.exists(row$file_path)) {
    skipped <- c(skipped, sprintf("%s: file not found — %s", row$metric_id, row$file_path))
    next
  }

  df <- read.csv(row$file_path, stringsAsFactors = FALSE)

  if (!is.na(row$filter_col) && row$filter_col != "" && row$filter_col %in% names(df)) {
    df <- df[df[[row$filter_col]] == row$filter_value, , drop = FALSE]
  }

  observed <- tryCatch(compute_metric(df, row$aggregate, row$column),
                       error = function(e) NA)

  if (is.na(observed) || !is.numeric(observed)) {
    skipped <- c(skipped, sprintf("%s: observed value is NA or non-numeric", row$metric_id))
    next
  }

  expected <- as.numeric(row$expected_value)
  tol_abs  <- as.numeric(row$tol_abs)
  tol_rel  <- as.numeric(row$tol_rel)
  abs_err  <- abs(observed - expected)
  rel_err  <- ifelse(expected == 0, abs_err, abs_err / abs(expected))

  if (abs_err <= tol_abs || rel_err <= tol_rel) {
    metrics$status[idx] <- "canonical"
    metrics$notes[idx]  <- paste0("Promoted to canonical on ", today,
                                   " (observed=", signif(observed, 7), ")")
    promoted <- c(promoted, row$metric_id)
  } else {
    # Update expected value to observed and leave as pending for manual review
    skipped <- c(skipped, sprintf(
      "%s: MISMATCH — observed=%s expected=%s abs_err=%s rel_err=%s",
      row$metric_id, signif(observed, 7), expected, signif(abs_err, 4), signif(rel_err, 4)
    ))
  }
}

write.csv(metrics, metrics_path, row.names = FALSE, quote = TRUE)

if (length(promoted) > 0) {
  message("[promote-metrics] Promoted ", length(promoted), " metric(s) to canonical:")
  for (m in promoted) message("  + ", m)
}
if (length(skipped) > 0) {
  message("[promote-metrics] Skipped ", length(skipped), " metric(s):")
  for (m in skipped) message("  ! ", m)
}
if (length(promoted) == 0 && length(skipped) == 0) {
  message("[promote-metrics] No pending metrics to evaluate.")
}
