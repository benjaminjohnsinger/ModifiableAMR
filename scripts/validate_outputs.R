#!/usr/bin/env Rscript

source("config.R")

args <- commandArgs(trailingOnly = TRUE)
scenario <- if (length(args) >= 1) args[[1]] else get_amr_scenario()
producer_filter <- if (length(args) >= 2) args[[2]] else ""

fig_manifest <- read.csv(file.path(AMR_CONFIG$paths$manifests, "figure_manifest.csv"))
expected <- fig_manifest[fig_manifest$scenario == scenario, , drop = FALSE]

expected <- expected[expected$status == "canonical", , drop = FALSE]

if (!identical(producer_filter, "")) {
  expected <- expected[expected$producer == producer_filter, , drop = FALSE]
}

if (nrow(expected) == 0) {
  message("[validate] No manifest outputs defined for scenario: ", scenario)
  quit(status = 0)
}

missing_inputs <- expected$input_hint[!is.na(expected$input_hint) & expected$input_hint != "" & !file.exists(expected$input_hint)]
if (length(missing_inputs) > 0) {
  message("[validate] Missing expected input files declared in manifest:")
  for (path in unique(missing_inputs)) message("  - ", path)
  quit(status = 1)
}

missing <- expected$output_path[!file.exists(expected$output_path)]
if (length(missing) > 0) {
  message("[validate] Missing outputs:")
  for (path in missing) message("  - ", path)
  quit(status = 1)
}

message("[validate] All expected outputs found for scenario: ", scenario)
