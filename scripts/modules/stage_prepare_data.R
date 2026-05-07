source("config.R")

run_prepare_data_stage <- function(scenario = get_amr_scenario()) {
  message("[prepare_data] scenario: ", scenario)

  if (scenario %in% c("main", "hic", "lmic", "raw_iqvia")) {
    require_inputs(c(
      "pathogen_abx_analysis_all_variables_(class-specific).csv",
      "ATLAS_data/ATLAS_data_renamed.csv",
      "ATLAS_more/ATLAS_more_renamed.csv",
      "ATLAS_Enterococcus/ATLAS_Enterococcus_renamed.csv",
      "GASP_N_renamed.csv",
      "DDD_country_year_class.csv",
      "Chungman/Chungman_pca_renamed.csv"
    ), stage = "prepare_data")
  } else if (scenario == "hospital_nagorsen") {
    require_inputs(c(
      "Nagorsen_clean.csv",
      "Chungman/Chungman_pca_renamed.csv"
    ), stage = "prepare_data")
  }

  if (scenario %in% c("main", "hic", "lmic", "raw_iqvia")) {
    message("[prepare_data] Running main regression data preparation...")
    source("data_processing.r")
    prepare_main_regression_data()
  } else if (scenario == "hospital_nagorsen") {
    message("[prepare_data] Running hospital regression data preparation...")
    source("data_processing.r")
    prepare_nagorsen_hospital_regression_data()
  } else {
    message("[prepare_data] No data-preparation action defined for scenario: ", scenario)
  }

  message("[prepare_data] done")
}
