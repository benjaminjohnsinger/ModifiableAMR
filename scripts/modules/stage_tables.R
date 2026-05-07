source("config.R")

run_tables_stage <- function(scenario = NULL) {
  if (is.null(scenario)) {
    scenario <- get("get_amr_scenario", envir = .GlobalEnv)()
  }
  cfg <- get("AMR_CONFIG", envir = .GlobalEnv)

  message("[tables] scenario: ", scenario)
  message("[tables] Note: supplementary tables are not part of the default manuscript target.")
  message("[tables]   Tables 2 and 3 are generated here; Table 2 requires permutation model")
  message("[tables]   outputs (see: make permutation-models).")

  if (scenario == "main") {
    message("[tables] Running Wilcoxon permutation test (Table 2 prerequisite)...")
    tryCatch(
      source("scripts/run_wilcoxon_test.R"),
      error = function(e) {
        if (conditionMessage(e) == ".wilcoxon_skip") {
          # graceful skip — message already printed
        } else {
          stop(e)
        }
      }
    )

    message("[tables] Generating Supplementary Table 2...")
    source("scripts/generate_supplementary_table2.R")

    message("[tables] Generating Supplementary Table 3...")
    # Check if the gradient file has the R_squared columns required by Table 3
    .tab3_path <- getOption("amr_table3_gradients_path",
      "Outputs/database_gradients_pathogen_ATC3_PCA_canonical_weighted_main.csv")
    if (!file.exists(.tab3_path)) {
      message("[tables] Table 3 gradient file not found — skipping: ", .tab3_path)
    } else {
      .tab3_cols <- names(data.table::fread(.tab3_path, nrows = 0))
      if (!"R_squared" %in% .tab3_cols) {
        message("[tables] Table 3: R_squared column not present in gradient file — skipping.")
        message("[tables]   Re-run 'make models SCENARIO=main' to regenerate the gradient file.")
      } else {
        source(cfg$canonical$table_script)
      }
    }
  } else {
    message("[tables] No table action defined for scenario: ", scenario)
  }

  message("[tables] done")
}
