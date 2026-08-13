source("config.R")

run_generate_figures_stage <- function(scenario = get_amr_scenario()) {
  message("[generate_figures] scenario: ", scenario)

  # Resolve model-output tags for figure inputs.
  # Defaults preserve canonical main behavior.
  plot_tags <- list(
    pathogen_tag = "main",
    class_tag = "all",
    random_pathogen_tag = "all"
  )
  if (scenario == "main_binomial") {
    plot_tags <- list(
      pathogen_tag = "main_binomial",
      class_tag = "all_binomial",
      random_pathogen_tag = "all_binomial"
    )
  } else if (scenario == "continuity") {
    plot_tags <- list(
      pathogen_tag = "mainccorrected",
      class_tag = "allccorrected",
      random_pathogen_tag = "allccorrected"
    )
  } else if (scenario == "continuity6") {
    plot_tags <- list(
      pathogen_tag = "mainccorrected6",
      class_tag = "allccorrected6",
      random_pathogen_tag = "allccorrected6"
    )
  } else if (scenario == "continuity100") {
    plot_tags <- list(
      pathogen_tag = "mainccorrected100",
      class_tag = "allccorrected100",
      random_pathogen_tag = "allccorrected100"
    )
  } else if (scenario == "raw_iqvia") {
    plot_tags <- list(
      pathogen_tag = "main_iqvia",
      class_tag = "all_iqvia",
      random_pathogen_tag = "all_iqvia"
    )
  } else if (scenario == "main_ppml") {
    plot_tags <- list(
      pathogen_tag = "main_ppml",
      class_tag = "all_ppml",
      random_pathogen_tag = "all_ppml"
    )
  } else if (scenario == "main_2000") {
    plot_tags <- list(
      pathogen_tag = "main_2000",
      class_tag = "all_2000",
      random_pathogen_tag = "all_2000"
    )
  } else if (scenario == "main_2000_plot") {
    plot_tags <- list(
      pathogen_tag = "main_2000",
      class_tag = "all_2000",
      random_pathogen_tag = "all_2000"
    )
  } else if (scenario == "main_finer") {
    plot_tags <- list(
      pathogen_tag = "main_finer",
      class_tag = "all_finer",
      random_pathogen_tag = "all_finer"
    )
  } else if (scenario == "hic_ppml") {
    plot_tags <- list(
      pathogen_tag = "hic_ppml",
      class_tag = "hic_ppml",
      random_pathogen_tag = "hic_ppml"
    )
  } else if (scenario == "lmic_ppml") {
    plot_tags <- list(
      pathogen_tag = "lmic_ppml",
      class_tag = "lmic_ppml",
      random_pathogen_tag = "lmic_ppml"
    )
  } else if (scenario == "hic") {
    plot_tags <- list(
      pathogen_tag = "HIC",
      class_tag = "HIC",
      random_pathogen_tag = "HIC"
    )
  } else if (scenario == "lmic") {
    plot_tags <- list(
      pathogen_tag = "LMIC",
      class_tag = "LMIC",
      random_pathogen_tag = "LMIC"
    )
  } else if (scenario == "exploratory_lagged") {
    # Match the custom lag parsing from the model script
    custom_lag <- as.integer(Sys.getenv("AMR_LAG_N", unset = "1"))
    if (is.na(custom_lag)) custom_lag <- 1
    
    tag_suffix <- if (custom_lag == 1) "lagged" else paste0("lagged_", custom_lag, "y")
    
    plot_tags <- list(
      pathogen_tag = tag_suffix,
      class_tag = paste0("all_", tag_suffix),
      random_pathogen_tag = paste0("all_", tag_suffix)
    )
  } else if (scenario == "consumption_lagged") {
    # Match the custom lag parsing from the model script
    custom_lag <- as.integer(Sys.getenv("AMR_LAG_N", unset = "1"))
    if (is.na(custom_lag)) custom_lag <- 1
    
    tag_suffix <- if (custom_lag == 1) "clagged" else paste0("clagged_", custom_lag, "y")
    
    plot_tags <- list(
      pathogen_tag = tag_suffix,
      class_tag = paste0("all_", tag_suffix),
      random_pathogen_tag = paste0("all_", tag_suffix)
    )
  } else if (scenario == "consumption_lagged_ppml") {
    custom_lag <- as.integer(Sys.getenv("AMR_LAG_N", unset = "1"))
    if (is.na(custom_lag)) custom_lag <- 1

    tag_suffix <- if (custom_lag == 1) "clagged_ppml" else paste0("clagged_ppml_", custom_lag, "y")

    plot_tags <- list(
      pathogen_tag = tag_suffix,
      class_tag = paste0("all_", tag_suffix),
      random_pathogen_tag = paste0("all_", tag_suffix)
    )
  } else if (scenario == "extra_pcs") {
    plot_tags <- list(
      pathogen_tag = "extra_pcs",
      class_tag = "all_extra_pcs",
      random_pathogen_tag = "all_extra_pcs"
    )
  } else if (scenario == "extra_pcs_ppml") {
    plot_tags <- list(
      pathogen_tag = "extra_pcs_ppml",
      class_tag = "all_extra_pcs_ppml",
      random_pathogen_tag = "all_extra_pcs_ppml"
    )
  } else if (scenario == "mi") {
    plot_tags <- list(
      pathogen_tag = "mi",
      class_tag = "all_mi",
      random_pathogen_tag = "all_mi"
    )
  }

  fig1_pathogen_input <- paste0(
    "Outputs/database_gradients_pathogen_ATC3_PCA_canonical_weighted_",
    plot_tags$pathogen_tag,
    ".csv"
  )
  fig1_class_gradients_input <- paste0(
    "Outputs/database_gradients_ATC3_PCA_canonical_weighted_",
    plot_tags$class_tag,
    ".csv"
  )
  fig1_class_bootstrap_input <- paste0(
    "Outputs/database_gradients_bootstraps_ATC3_PCA_canonical_weighted_",
    plot_tags$class_tag,
    ".csv"
  )
  fig2_pathogen_input <- paste0(
    "Outputs/database_gradients_pathogen_PCA_canonical_weighted_",
    plot_tags$random_pathogen_tag,
    ".csv"
  )
  fig2_pathogen_bootstrap_input <- paste0(
    "Outputs/database_gradients_bootstraps_pathogen_PCA_canonical_weighted_",
    plot_tags$random_pathogen_tag,
    ".csv"
  )

  message(
    "[generate_figures] figure tags: pathogen_tag=", plot_tags$pathogen_tag,
    ", class_tag=", plot_tags$class_tag,
    ", random_pathogen_tag=", plot_tags$random_pathogen_tag
  )
  message("[generate_figures] Figure1 pathogen input: ", fig1_pathogen_input)
  message("[generate_figures] Figure1 class bootstrap input: ", fig1_class_bootstrap_input)
  message("[generate_figures] Figure2 pathogen input: ", fig2_pathogen_input)
  message("[generate_figures] Figure2 pathogen bootstrap input: ", fig2_pathogen_bootstrap_input)

  # Added exploratory_lagged to require the correct inputs
  if (scenario %in% c("main", "continuity", "continuity6", "continuity100",  "main_2000", "raw_iqvia", "main_2000_plot", "main_finer", "main_binomial", "main_ppml", "hic", "hic_ppml", "lmic", "lmic_ppml", "exploratory_lagged", "consumption_lagged", "consumption_lagged_ppml", "extra_pcs", "extra_pcs_ppml", "mi")) {
    require_inputs(c(
      fig1_pathogen_input,
      fig1_class_gradients_input,
      fig1_class_bootstrap_input,
      fig2_pathogen_input,
      fig2_pathogen_bootstrap_input
    ), stage = "generate_figures")
  } else if (scenario == "hospital_nagorsen") {
    require_input(
      "Outputs/Nagorsen_gradients_pathogen_ATC3_PCA_canonical_hospital_to_all_filtered.csv",
      stage = "generate_figures"
    )
  }

  if (scenario %in% c("main", "continuity", "continuity6", "continuity100", "main_2000", "raw_iqvia", "main_2000_plot", "main_finer", "hic", "lmic", "main_ppml", "hic_ppml", "lmic_ppml")) {
    burden_files_present <- warn_missing_inputs(c(
      AMR_CONFIG$burden_inputs$figure3_pathogen,
      AMR_CONFIG$burden_inputs$figure3_optimistic,
      AMR_CONFIG$burden_inputs$figure3_pessimistic
    ), stage = "generate_figures")
  } else {
    burden_files_present <- FALSE
  }

  dir.create(AMR_CONFIG$output_dirs$manuscript, recursive = TRUE, showWarnings = FALSE)
  dir.create(AMR_CONFIG$output_dirs$slides, recursive = TRUE, showWarnings = FALSE)

  # Added exploratory_lagged to actually run the figure module
  if (scenario %in% c("main", "continuity", "continuity6", "continuity100", "main_2000", "raw_iqvia", "main_2000_plot", "main_finer", "main_binomial", "main_ppml", "hic", "hic_ppml", "lmic", "lmic_ppml", "hospital_nagorsen", "extra_pcs", "extra_pcs_ppml", "exploratory_lagged", "consumption_lagged", "consumption_lagged_ppml", "mi")) {
    message("[generate_figures] Running canonical figure module...")
    old_options <- options(
      amr_plot_pathogen_tag = plot_tags$pathogen_tag,
      amr_plot_class_tag = plot_tags$class_tag,
      amr_plot_random_pathogen_tag = plot_tags$random_pathogen_tag
    )
    on.exit(options(old_options), add = TRUE)
    source(AMR_CONFIG$canonical$figure_module)
    if (scenario == "main_finer") {
      write_figure_metadata("Figure1_finer_part1.pdf",
        inputs = c(fig1_pathogen_input,
                   fig1_class_bootstrap_input),
        scenario = scenario)
      write_figure_metadata("Figure1_finer_part2.pdf",
        inputs = c(fig1_pathogen_input,
                   fig1_class_bootstrap_input),
        scenario = scenario)
    } else {
      write_figure_metadata("Figure1.pdf",
        inputs = c(fig1_pathogen_input,
                   fig1_class_bootstrap_input),
        scenario = scenario)
    }
    write_figure_metadata("Figure2.pdf",
      inputs = c(fig2_pathogen_input,
                 fig2_pathogen_bootstrap_input),
      scenario = scenario)
    if (burden_files_present) {
      write_figure_metadata("Figure3.pdf",
        inputs = c(AMR_CONFIG$burden_inputs$figure3_pathogen,
                   AMR_CONFIG$burden_inputs$figure3_optimistic,
                   AMR_CONFIG$burden_inputs$figure3_pessimistic),
        scenario = scenario)
    }
    # Figure 4 is generated by plotting.R (generate_figure4()) when its burden
    # input CSVs are present (written by make burden).
    fig4_inputs_present <- warn_missing_inputs(c(
      AMR_CONFIG$burden_inputs$figure4_region,
      AMR_CONFIG$burden_inputs$figure4_gdp,
      AMR_CONFIG$burden_inputs$figure4_use,
      AMR_CONFIG$burden_inputs$lower_burden_region
    ), stage = "generate_figures")
    if (fig4_inputs_present) {
      write_figure_metadata("Figure4.pdf",
        inputs = c(AMR_CONFIG$burden_inputs$figure4_region,
                   AMR_CONFIG$burden_inputs$figure4_gdp,
                   AMR_CONFIG$burden_inputs$figure4_use),
        scenario = scenario)
    }
    # Supplementary Figure S1: class-level elasticities HIC vs LMIC.
    # For main_ppml, use ppml-tagged subgroup inputs.
    suppfig_s1_hic_suffix <- "HIC"
    suppfig_s1_lmic_suffix <- "LMIC"
    if (scenario == "main_ppml") {
      suppfig_s1_hic_suffix <- "HIC_ppml"
      suppfig_s1_lmic_suffix <- "LMIC_ppml"
    }
    suppfig_s1_hic_input <- paste0(
      "Outputs/database_gradients_ATC3_PCA_canonical_weighted_",
      suppfig_s1_hic_suffix,
      ".csv"
    )
    suppfig_s1_lmic_input <- paste0(
      "Outputs/database_gradients_ATC3_PCA_canonical_weighted_",
      suppfig_s1_lmic_suffix,
      ".csv"
    )

    suppfig_s1_inputs_present <- warn_missing_inputs(c(
      suppfig_s1_hic_input,
      suppfig_s1_lmic_input
    ), stage = "generate_figures")
    if (suppfig_s1_inputs_present) {
      write_figure_metadata(
        file.path(AMR_CONFIG$output_dirs$slides, "Supplementary_Figure_S1_Slide_narrow.pdf"),
        inputs = c(
          suppfig_s1_hic_input,
          suppfig_s1_lmic_input
        ),
        scenario = scenario)
    }
  }

  message("[generate_figures] done")
}