source("config.R")

resolve_burden_model_tags <- function(scenario) {
  tags <- list(
    pathogen_tag = "main",
    class_tag = "all",
    use_nagorsen_prefix = FALSE
  )

  if (scenario == "continuity") {
    custom_correction <- as.integer(Sys.getenv("CONTINUITY_EXP", unset = "2"))
    if (is.na(custom_correction)) custom_correction <- 2L
    tag_suffix <- if (custom_correction == 2L) "ccorrected" else paste0("ccorrected", custom_correction)
    tags <- list(
      pathogen_tag = paste0("main", tag_suffix),
      class_tag = paste0("all", tag_suffix),
      use_nagorsen_prefix = FALSE
    )
  } else if (scenario == "main_2000") {
    tags <- list(pathogen_tag = "main_2000", class_tag = "all_2000", use_nagorsen_prefix = FALSE)
  } else if (scenario == "main_finer") {
    tags <- list(pathogen_tag = "main_finer", class_tag = "all", use_nagorsen_prefix = FALSE)
  } else if (scenario == "main_binomial") {
    tags <- list(pathogen_tag = "main_binomial", class_tag = "all_binomial", use_nagorsen_prefix = FALSE)
  } else if (scenario == "main_ppml") {
    tags <- list(pathogen_tag = "main_ppml", class_tag = "all_ppml", use_nagorsen_prefix = FALSE)
  } else if (scenario == "hic") {
    tags <- list(pathogen_tag = "HIC", class_tag = "HIC", use_nagorsen_prefix = FALSE)
  } else if (scenario == "hic_ppml") {
    tags <- list(pathogen_tag = "HIC_ppml", class_tag = "HIC_ppml", use_nagorsen_prefix = FALSE)
  } else if (scenario == "lmic") {
    tags <- list(pathogen_tag = "LMIC", class_tag = "LMIC", use_nagorsen_prefix = FALSE)
  } else if (scenario == "lmic_ppml") {
    tags <- list(pathogen_tag = "LMIC_ppml", class_tag = "LMIC_ppml", use_nagorsen_prefix = FALSE)
  } else if (scenario == "raw_iqvia") {
    tags <- list(pathogen_tag = "main_IQVIA", class_tag = "all_IQVIA", use_nagorsen_prefix = FALSE)
  } else if (scenario == "exploratory_lagged") {
    custom_lag <- as.integer(Sys.getenv("AMR_LAG_N", unset = "1"))
    if (is.na(custom_lag)) custom_lag <- 1L
    tag_suffix <- if (custom_lag == 1L) "lagged" else paste0("lagged_", custom_lag, "y")
    tags <- list(
      pathogen_tag = tag_suffix,
      class_tag = paste0("all_", tag_suffix),
      use_nagorsen_prefix = FALSE
    )
  } else if (scenario == "consumption_lagged") {
    custom_lag <- as.integer(Sys.getenv("AMR_LAG_N", unset = "1"))
    if (is.na(custom_lag)) custom_lag <- 1L
    tag_suffix <- if (custom_lag == 1L) "clagged" else paste0("clagged_", custom_lag, "y")
    tags <- list(
      pathogen_tag = tag_suffix,
      class_tag = paste0("all_", tag_suffix),
      use_nagorsen_prefix = FALSE
    )
  } else if (scenario == "consumption_lagged_ppml") {
    custom_lag <- as.integer(Sys.getenv("AMR_LAG_N", unset = "1"))
    if (is.na(custom_lag)) custom_lag <- 1L
    tag_suffix <- if (custom_lag == 1L) "clagged_ppml" else paste0("clagged_ppml_", custom_lag, "y")
    tags <- list(
      pathogen_tag = tag_suffix,
      class_tag = paste0("all_", tag_suffix),
      use_nagorsen_prefix = FALSE
    )
  } else if (scenario == "extra_pcs") {
    tags <- list(pathogen_tag = "extra_pcs", class_tag = "all_extra_pcs", use_nagorsen_prefix = FALSE)
  } else if (scenario == "extra_pcs_ppml") {
    tags <- list(pathogen_tag = "extra_pcs_ppml", class_tag = "all_extra_pcs_ppml", use_nagorsen_prefix = FALSE)
  } else if (scenario == "mi") {
    tags <- list(pathogen_tag = "mi", class_tag = "all_mi", use_nagorsen_prefix = FALSE)
  } else if (scenario == "hospital_nagorsen") {
    tags <- list(
      pathogen_tag = "hospital_to_all_filtered",
      class_tag = "hospital_to_all_filtered",
      use_nagorsen_prefix = TRUE
    )
  } else if (scenario == "permutation") {
    perm_class <- Sys.getenv("AMR_PERMUTATION_CLASS", unset = "")
    if (identical(perm_class, "")) {
      stop("[estimate_burden] permutation scenario requires AMR_PERMUTATION_CLASS to be set (e.g. J01A)", call. = FALSE)
    }
    tags <- list(pathogen_tag = paste0("permutation", perm_class), class_tag = NA_character_, use_nagorsen_prefix = FALSE)
  } else if (scenario == "permutations_ppml") {
    perm_class <- Sys.getenv("AMR_PERMUTATION_CLASS", unset = "")
    if (identical(perm_class, "")) {
      stop("[estimate_burden] permutations_ppml scenario requires AMR_PERMUTATION_CLASS to be set (e.g. J01A)", call. = FALSE)
    }
    tags <- list(pathogen_tag = paste0("permutations_ppml", perm_class), class_tag = NA_character_, use_nagorsen_prefix = FALSE)
  }

  tags
}

run_estimate_burden_stage <- function(scenario = get_amr_scenario()) {
  message("[estimate_burden] scenario: ", scenario)

  filter_implausible_env <- Sys.getenv("AMR_BURDEN_FILTER_IMPLAUSIBLE_GRADIENTS", unset = "1")
  filter_implausible <- !tolower(filter_implausible_env) %in% c("0", "false", "no")
  max_abs_gradient_env <- Sys.getenv("AMR_BURDEN_MAX_ABS_GRADIENT", unset = "10")
  max_abs_gradient <- suppressWarnings(as.numeric(max_abs_gradient_env))
  if (is.na(max_abs_gradient) || max_abs_gradient <= 0) {
    stop(
      "[estimate_burden] AMR_BURDEN_MAX_ABS_GRADIENT must be a positive number. Got: ",
      max_abs_gradient_env,
      call. = FALSE
    )
  }

  elasticity_scale <- 1
  if (scenario == "scale_up") {
    elasticity_scale <- 2
  } else if (scenario == "scale_down") {
    elasticity_scale <- 0.5
  }

  tags <- resolve_burden_model_tags(scenario)

  if (isTRUE(tags$use_nagorsen_prefix)) {
    gradients_path <- paste0("Outputs/Nagorsen_gradients_ATC3_PCA_canonical_", tags$class_tag, ".csv")
    gradients_bootstrap_path <- paste0("Outputs/Nagorsen_gradients_bootstraps_ATC3_PCA_canonical_", tags$class_tag, ".csv")
    results_path <- paste0("Outputs/Nagorsen_gradients_pathogen_ATC3_PCA_canonical_", tags$pathogen_tag, ".csv")
    results_bootstrap_path <- paste0("Outputs/Nagorsen_gradients_bootstraps_pathogen_ATC3_PCA_canonical_", tags$pathogen_tag, ".csv")
  } else {
    gradients_path <- if (!is.na(tags$class_tag)) {
      paste0("Outputs/database_gradients_ATC3_PCA_canonical_weighted_", tags$class_tag, ".csv")
    } else {
      NA_character_
    }
    gradients_bootstrap_path <- if (!is.na(tags$class_tag)) {
      paste0("Outputs/database_gradients_bootstraps_ATC3_PCA_canonical_weighted_", tags$class_tag, ".csv")
    } else {
      NA_character_
    }
    results_path <- paste0("Outputs/database_gradients_pathogen_ATC3_PCA_canonical_weighted_", tags$pathogen_tag, ".csv")
    results_bootstrap_path <- paste0("Outputs/database_gradients_bootstraps_pathogen_ATC3_PCA_canonical_weighted_", tags$pathogen_tag, ".csv")
  }

  model_like_scenarios <- c(
    "main", "continuity", "main_finer", "main_binomial", "main_ppml", "main_2000",
    "hic", "hic_ppml", "lmic", "lmic_ppml", "raw_iqvia", "hospital_nagorsen",
    "exploratory_lagged", "consumption_lagged", "consumption_lagged_ppml",
    "extra_pcs", "extra_pcs_ppml", "mi", "scale_up", "scale_down",
    "permutation", "permutations_ppml"
  )

  should_run_burden <- grepl("^burden_", scenario) || scenario %in% model_like_scenarios

  if (should_run_burden) {
    if (is.na(tags$class_tag)) {
      message(
        "[estimate_burden] Scenario '", scenario,
        "' does not emit class-level gradients required by burden estimation; skipping burden run."
      )
      message("[estimate_burden] done")
      return(invisible(NULL))
    }

    required_inputs <- c(
      "IHME_AMR/IHME_AMR_fitted_gammas_v2.csv",
      "IHME_AMR/IHME_AMR_PATHOGEN_2019_DATA_COUNTED_AB.CSV",
      results_path,
      results_bootstrap_path
    )
    if (!is.na(gradients_path) && !is.na(gradients_bootstrap_path)) {
      required_inputs <- c(required_inputs, gradients_path, gradients_bootstrap_path)
    }
    require_inputs(required_inputs, stage = "estimate_burden")
  }

  if (should_run_burden) {
    message("[estimate_burden] Running burden estimation script...")
    if (scenario != "main") {
      message("[estimate_burden] Non-main scenario detected; canonical burden outputs will not be overwritten.")
    }
    # Pass the resolved paths and scenario as global options before sourcing
    options(
      amr_scenario                        = scenario,
      amr_smoke_mode                      = identical(Sys.getenv("AMR_DEV_SMOKE", "0"), "1"),
      amr_burden_gradients_path           = gradients_path,
      amr_burden_gradients_bootstrap_path = gradients_bootstrap_path,
      amr_burden_results_path             = results_path,
      amr_burden_results_bootstrap_path   = results_bootstrap_path,
      amr_burden_filter_implausible_gradients = filter_implausible,
      amr_burden_max_abs_gradient         = max_abs_gradient,
      amr_burden_elasticity_scale         = elasticity_scale
    )
    message(
      "[estimate_burden] Gradient filter: enabled=", filter_implausible,
      ", max_abs_gradient=", max_abs_gradient
    )
    message("[estimate_burden] Elasticity scale multiplier: ", elasticity_scale)
    source(AMR_CONFIG$canonical$burden_script)
  } else {
    message("[estimate_burden] No burden-estimation action defined for scenario: ", scenario)
  }

  message("[estimate_burden] done")
}
