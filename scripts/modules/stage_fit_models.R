source("config.R")

run_fit_models_stage <- function(scenario = get_amr_scenario()) {
  message("[fit_models] scenario: ", scenario)

  if (scenario %in% c("main", "continuity", "hic", "lmic", "main_binomial", "main_ppml", "hic_ppml", "lmic_ppml", "exploratory_lagged", "consumption_lagged", "consumption_lagged_ppml", "extra_pcs", "extra_pcs_ppml", "mi", "main_2000")) {
    require_input("merged_data_N_PC3_GDP.csv", stage = "fit_models")
  } else if (scenario == "main_finer") {
    require_input("finer_data_new.csv", stage = "fit_models")
  } else if (scenario == "raw_iqvia") {
    require_input("merged_data_N_PC3_GDP_IQVIA.csv", stage = "fit_models")
  } else if (scenario == "hospital_nagorsen") {
    require_input("merged_data_Nagorsen_hospital_to_all_filtered.csv", stage = "fit_models")
  } else if (scenario %in% c("permutation", "permutations_ppml")) {
    require_input("merged_data_N_PC3_GDP.csv", stage = "fit_models")
    perm_class <- Sys.getenv("AMR_PERMUTATION_CLASS", unset = "")
    if (identical(perm_class, "")) {
      stop("[fit_models] permutation scenario requires AMR_PERMUTATION_CLASS to be set (e.g. J01A)",
           call. = FALSE)
    }
    message("[fit_models] Permutation class: ", perm_class)
  }

  if (scenario %in% c("main", "continuity", "main_finer", "main_binomial", "main_ppml", "hic", "hic_ppml", "lmic", "lmic_ppml", "raw_iqvia", "hospital_nagorsen",
                       "exploratory_lagged", "consumption_lagged", "consumption_lagged_ppml", "extra_pcs", "extra_pcs_ppml", "mi", "permutation", "permutations_ppml", "main_2000")) {
    message("[fit_models] Running linear model script...")
    # Keep sourced model script scenario resolution aligned with this stage argument.
    Sys.setenv(ANALYSIS_SCENARIO = scenario, SCENARIO = scenario, AMR_SCENARIO = scenario)
    source(AMR_CONFIG$canonical$model_script)
  } else {
    message("[fit_models] No model-fitting action defined for scenario: ", scenario)
  }

  message("[fit_models] done")
}
