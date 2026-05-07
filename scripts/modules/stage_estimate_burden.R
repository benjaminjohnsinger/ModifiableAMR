source("config.R")

run_estimate_burden_stage <- function(scenario = get_amr_scenario()) {
  message("[estimate_burden] scenario: ", scenario)

  # Resolve output tags matching those used by regression_models.r
  tags <- list(
    pathogen_tag = "main",
    class_tag    = "all"
  )
  if (scenario == "hic")  tags <- list(pathogen_tag = "HIC",  class_tag = "HIC")
  if (scenario == "lmic") tags <- list(pathogen_tag = "LMIC", class_tag = "LMIC")

  gradients_path           <- paste0("Outputs/database_gradients_ATC3_PCA_canonical_weighted_",           tags$class_tag,    ".csv")
  gradients_bootstrap_path <- paste0("Outputs/database_gradients_bootstraps_ATC3_PCA_canonical_weighted_", tags$class_tag,    ".csv")
  results_path             <- paste0("Outputs/database_gradients_pathogen_ATC3_PCA_canonical_weighted_",  tags$pathogen_tag, ".csv")
  results_bootstrap_path   <- paste0("Outputs/database_gradients_bootstraps_pathogen_ATC3_PCA_canonical_weighted_", tags$pathogen_tag, ".csv")

  if (grepl("^burden_", scenario) || scenario %in% c("main", "hic", "lmic")) {
    require_inputs(c(
      "IHME_AMR/IHME_AMR_fitted_gammas_v2.csv",
      "IHME_AMR/IHME_AMR_PATHOGEN_2019_DATA_COUNTED_AB.CSV",
      gradients_path,
      results_path,
      gradients_bootstrap_path
    ), stage = "estimate_burden")
  }

  if (grepl("^burden_", scenario) || scenario %in% c("main", "hic", "lmic")) {
    message("[estimate_burden] Running burden estimation script...")
    # Pass the resolved paths and scenario as global options before sourcing
    options(
      amr_scenario                        = scenario,
      amr_smoke_mode                      = identical(Sys.getenv("AMR_DEV_SMOKE", "0"), "1"),
      amr_burden_gradients_path           = gradients_path,
      amr_burden_gradients_bootstrap_path = gradients_bootstrap_path,
      amr_burden_results_path             = results_path,
      amr_burden_results_bootstrap_path   = results_bootstrap_path
    )
    source(AMR_CONFIG$canonical$burden_script)
  } else {
    message("[estimate_burden] No burden-estimation action defined for scenario: ", scenario)
  }

  message("[estimate_burden] done")
}
