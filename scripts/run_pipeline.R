#!/usr/bin/env Rscript

source("config.R")

args <- commandArgs(trailingOnly = TRUE)
scenario <- if (length(args) >= 1) args[[1]] else get_amr_scenario()
Sys.setenv(ANALYSIS_SCENARIO = scenario)
Sys.setenv(AMR_SCENARIO = scenario)

message("[pipeline] scenario: ", scenario)

source("scripts/run_prepare_data.R")
source("scripts/run_fit_models.R")
source("scripts/run_estimate_burden.R")
source("scripts/run_generate_figures.R")
source("scripts/run_tables.R")

message("[pipeline] complete")
