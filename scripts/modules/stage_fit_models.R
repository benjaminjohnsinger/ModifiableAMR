source("config.R")

run_fit_models_stage <- function(scenario = get_amr_scenario()) {
  message("[fit_models] scenario: ", scenario)

  if (scenario %in% c("main", "hic", "lmic")) {
    require_input("merged_data_N_PC3_GDP.csv", stage = "fit_models")
  } else if (scenario == "raw_iqvia") {
    require_input("merged_data_N_PC3_GDP_IQVIA.csv", stage = "fit_models")
  } else if (scenario == "hospital_nagorsen") {
    require_input("merged_data_Nagorsen_hospital_to_all_filtered.csv", stage = "fit_models")
  } else if (scenario == "permutation") {
    require_input("merged_data_N_PC3_GDP.csv", stage = "fit_models")
    perm_class <- Sys.getenv("AMR_PERMUTATION_CLASS", unset = "")
    if (identical(perm_class, "")) {
      stop("[fit_models] permutation scenario requires AMR_PERMUTATION_CLASS to be set (e.g. J01A)",
           call. = FALSE)
    }
    message("[fit_models] Permutation class: ", perm_class)
  }

  if (scenario %in% c("main", "hic", "lmic", "raw_iqvia", "hospital_nagorsen",
                       "exploratory_lagged", "permutation")) {
    message("[fit_models] Running linear model script...")
    source(AMR_CONFIG$canonical$model_script)
  } else {
    message("[fit_models] No model-fitting action defined for scenario: ", scenario)
  }

  message("[fit_models] done")
}
