#!/usr/bin/env Rscript

source("config.R")

args <- commandArgs(trailingOnly = TRUE)
scenario <- if (length(args) >= 1) args[[1]] else get_amr_scenario()

metrics_path <- file.path(AMR_CONFIG$paths$manifests, "expected_metrics.csv")
if (!file.exists(metrics_path)) {
  message("[validate-metrics] Metrics manifest not found. Skipping.")
  quit(status = 0)
}

metrics <- read.csv(metrics_path, stringsAsFactors = FALSE)
if (nrow(metrics) == 0) {
  message("[validate-metrics] No metrics defined. Skipping.")
  quit(status = 0)
}

if (!"scenario" %in% names(metrics)) {
  metrics$scenario <- "main"
}

active <- metrics[
  metrics$status == "canonical" & (metrics$scenario == scenario | metrics$scenario == "all"),
  ,
  drop = FALSE
]
if (nrow(active) == 0) {
  message("[validate-metrics] No canonical metrics defined for scenario: ", scenario, ". Skipping.")
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

failures <- character(0)

for (i in seq_len(nrow(active))) {
  row <- active[i, ]
  if (!file.exists(row$file_path)) {
    failures <- c(failures, sprintf("%s: file not found %s", row$metric_id, row$file_path))
    next
  }

  df <- read.csv(row$file_path, stringsAsFactors = FALSE)

  if (!is.na(row$filter_col) && row$filter_col != "" && row$filter_col %in% names(df)) {
    df <- df[df[[row$filter_col]] == row$filter_value, , drop = FALSE]
  }

  observed <- compute_metric(df, row$aggregate, row$column)
  expected <- as.numeric(row$expected_value)
  tol_abs <- ifelse(is.na(row$tol_abs), AMR_CONFIG$parameters$tolerance_abs, as.numeric(row$tol_abs))
  tol_rel <- ifelse(is.na(row$tol_rel), AMR_CONFIG$parameters$tolerance_rel, as.numeric(row$tol_rel))

  if (is.na(observed) || !is.numeric(observed)) {
    failures <- c(failures, sprintf(
      "%s: observed value is NA or non-numeric (got: %s)", row$metric_id, observed
    ))
    next
  }

  abs_err <- abs(observed - expected)
  rel_err <- ifelse(expected == 0, abs_err, abs_err / abs(expected))

  if (!(abs_err <= tol_abs || rel_err <= tol_rel)) {
    failures <- c(
      failures,
      sprintf(
        "%s: observed=%s expected=%s abs_err=%s rel_err=%s",
        row$metric_id, observed, expected, abs_err, rel_err
      )
    )
  }
}

if (length(failures) > 0) {
  message("[validate-metrics] Metric mismatches:")
  for (f in failures) message("  - ", f)
  quit(status = 1)
}

message("[validate-metrics] All canonical metrics within tolerance for scenario: ", scenario)
