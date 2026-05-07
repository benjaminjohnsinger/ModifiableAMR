# Central configuration for reproducible AMR analysis workflow

options(stringsAsFactors = FALSE)

AMR_CONFIG <- list(
  paths = list(
    root = ".",
    outputs = "Outputs",
    manifests = "manifests",
    logs = "logs"
  ),
  parameters = list(
    year_cutoff = 2018,
    min_entry_count = 10,
    min_isolates = 100,
    n_bootstraps = 1000,
    global_normalization_year = 2018,
    tolerance_abs = 1e-6,
    tolerance_rel = 1e-4
  ),
  scenarios = list(
    default = "main",
    all = c(
      "main",
      "hic",
      "lmic",
      "raw_iqvia",
      "hospital_nagorsen",
      "burden_optimistic",
      "burden_pessimistic",
      "burden_lower_region",
      "burden_upper_region",
      "burden_drug_region",
      "burden_pathogen_region",
      "permutation"
    ),
    manuscript = c("main")
  ),
  canonical = list(
    figure_module = "scripts/plotting_canonical.R",
    model_script = "regression_models.r",
    burden_script = "avertable_burden.R",
    data_script = "data_processing.r",
    table_script = "generate_supplementary_table3.R"
  ),
  # Declared input paths for burden figures — single source of truth
  burden_inputs = list(
    figure3_pathogen = "Outputs/10pc_avertable_burden_by_pathogen_canonical_weighted_lower_region_v2.csv",
    figure3_optimistic = "Outputs/10pc_avertable_burden_by_pathogen_canonical_weighted_upper_region_optimistic_overall.csv",
    figure3_pessimistic = "Outputs/10pc_avertable_burden_by_pathogen_canonical_weighted_upper_region_pessimistic_overall.csv",
    figure4_region = "Outputs/10pc_avertable_burden_by_region_canonical_weighted_lower_region_v2.csv",
    figure4_gdp = "Outputs/gdp_by_lower_ihme_region_2018_test.csv",
    figure4_use = "Outputs/use_by_lower_ihme_region_2018_test.csv",
    lower_burden_region = "Outputs/total_bacterial_disease_burden_by_lower_ihme_region_v2.csv"
  ),
  # Output directories for figure artifacts
  output_dirs = list(
    manuscript = "Outputs/manuscript",
    slides = "Outputs/slides",
    data = "Outputs"
  )
)

get_amr_scenario <- function(default = AMR_CONFIG$scenarios$default) {
  scenario <- Sys.getenv("ANALYSIS_SCENARIO", unset = "")
  if (identical(scenario, "")) {
    scenario <- Sys.getenv("SCENARIO", unset = "")
  }
  if (identical(scenario, "")) {
    scenario <- Sys.getenv("AMR_SCENARIO", unset = default)
  }
  if (!scenario %in% AMR_CONFIG$scenarios$all) {
    stop(
      sprintf(
        "Unknown scenario '%s'. Allowed scenarios: %s",
        scenario,
        paste(AMR_CONFIG$scenarios$all, collapse = ", ")
      )
    )
  }
  scenario
}

is_manuscript_target <- function() {
  identical(Sys.getenv("AMR_TARGET", unset = ""), "manuscript")
}

# Stage input validation ---------------------------------------------------

#' Fail with an actionable message if a required input file is missing.
require_input <- function(path, stage = "stage") {
  if (!file.exists(path)) {
    stop(sprintf(
      "[%s] Required input not found: %s\n  Ensure upstream stages have been run and the file exists at the expected path.",
      stage, path
    ), call. = FALSE)
  }
  invisible(path)
}

#' Check multiple required inputs, failing on the first missing file.
require_inputs <- function(paths, stage = "stage") {
  invisible(lapply(paths, require_input, stage = stage))
}

#' Check inputs and warn (not fail) if any are missing.
warn_missing_inputs <- function(paths, stage = "stage") {
  for (p in paths) {
    if (!file.exists(p)) {
      warning(sprintf("[%s] Optional input not found (skipping): %s", stage, p))
    }
  }
  invisible(all(file.exists(paths)))
}

# Figure metadata ----------------------------------------------------------

#' Write a small sidecar JSON-like text file alongside a figure recording
#' the producing scenario, input files, and generation timestamp.
write_figure_metadata <- function(figure_path, inputs = character(0),
                                  scenario = Sys.getenv("ANALYSIS_SCENARIO", "unknown")) {
  meta_path <- sub("\\.[^.]+$", ".meta.txt", figure_path)
  lines <- c(
    paste0("figure: ", basename(figure_path)),
    paste0("scenario: ", scenario),
    paste0("generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    if (length(inputs) > 0) paste0("inputs:\n", paste0("  - ", inputs, collapse = "\n")) else "inputs: []"
  )
  writeLines(lines, meta_path)
  invisible(meta_path)
}
