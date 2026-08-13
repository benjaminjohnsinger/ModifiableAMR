#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
  "Usage:",
  "Rscript scripts/export_model_fit_tables.R [--output-dir DIR] model_fit.csv [more_model_fit.csv ...]",
  sep = "\n"
)

if (length(args) == 0) {
  stop(usage, call. = FALSE)
}

output_dir <- NULL
if (length(args) >= 2 && args[[1]] == "--output-dir") {
  output_dir <- args[[2]]
  args <- args[-c(1, 2)]
}

if (length(args) == 0) {
  stop(usage, call. = FALSE)
}

input_paths <- args

script_file <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))])
script_dir <- if (length(script_file) > 0) dirname(normalizePath(script_file[[1]])) else getwd()
source(file.path(script_dir, "..", "utils.r"))

detect_model_scenario <- function(path) {
  name <- tolower(tools::file_path_sans_ext(basename(path)))

  if (grepl("clagged_3y", name, fixed = TRUE)) {
    return("clagged_3y")
  }
  if (grepl("clagged_5y", name, fixed = TRUE)) {
    return("clagged_5y")
  }
  if (grepl("extra_pcs", name, fixed = TRUE)) {
    return("extra_pcs")
  }
  if (grepl("clagged", name, fixed = TRUE) || grepl("consumption_lagged", name, fixed = TRUE)) {
    return("clagged")
  }
  if (grepl("main", name, fixed = TRUE)) {
    return("main")
  }

  NA_character_
}

format_response_ci <- function(response, lower_ci, upper_ci) {
  if (is.na(response) || is.na(lower_ci) || is.na(upper_ci)) {
    return("")
  }

  sprintf("%.2f (%.2f to %.2f)", response, lower_ci, upper_ci)
}

format_burden_ci <- function(avertable, lower_ci, upper_ci) {
  if (is.na(avertable) || is.na(lower_ci) || is.na(upper_ci)) {
    return("")
  }

  mean_txt <- formatC(round(avertable), format = "f", digits = 0)
  lower_txt <- formatC(round(lower_ci), format = "f", digits = 0)
  upper_txt <- formatC(round(upper_ci), format = "f", digits = 0)
  paste0(mean_txt, " (", lower_txt, " to ", upper_txt, ")")
}

build_response_table_from_file <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    return(NULL)
  }

  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  has_response_cols <- all(c("Antibiotic", "Pathogen", "Response", "Lower_CI", "Upper_CI") %in% names(df))
  if (!has_response_cols) {
    return(NULL)
  }

  response_df <- df[c("Antibiotic", "Pathogen", "Response", "Lower_CI", "Upper_CI")]
  response_df <- rename_atc_classes(response_df)
  response_df$Response_fmt <- vapply(
    seq_len(nrow(response_df)),
    function(i) format_response_ci(
      response_df$Response[[i]],
      response_df$Lower_CI[[i]],
      response_df$Upper_CI[[i]]
    ),
    character(1)
  )

  response_df[c("Antibiotic", "Pathogen", "Response_fmt")]
}

build_response_sensitivity_table <- function(response_tables, scenario_order, scenario_labels) {
  available <- Filter(Negate(is.null), response_tables)
  if (length(available) == 0) {
    return(NULL)
  }

  key_pairs <- unique(do.call(rbind, lapply(available, function(df) {
    df[c("Antibiotic", "Pathogen")]
  })))

  out <- key_pairs
  for (scenario_id in scenario_order) {
    column_name <- scenario_labels[[scenario_id]]
    scenario_df <- response_tables[[scenario_id]]
    if (is.null(scenario_df)) {
      out[[column_name]] <- ""
      next
    }

    merged <- merge(
      out[c("Antibiotic", "Pathogen")],
      scenario_df,
      by = c("Antibiotic", "Pathogen"),
      all.x = TRUE,
      sort = FALSE
    )

    vals <- merged$Response_fmt
    vals[is.na(vals)] <- ""
    out[[column_name]] <- vals
  }

  out[order(out$Antibiotic, out$Pathogen), c("Antibiotic", "Pathogen", unname(unlist(scenario_labels[scenario_order])))]
}

find_existing_file <- function(candidate_names, search_dirs) {
  for (dir_path in search_dirs) {
    for (candidate_name in candidate_names) {
      candidate_path <- file.path(dir_path, candidate_name)
      if (file.exists(candidate_path)) {
        return(candidate_path)
      }
    }
  }

  NULL
}

build_avertable_region_sensitivity_table <- function(search_dirs) {
  scenario_order <- c("main", "clagged", "clagged_3y", "clagged_5y", "extra_pcs", "scale_down", "scale_up")
  scenario_labels <- c(
    main = "Main",
    clagged = "clagged",
    clagged_3y = "clagged_3y",
    clagged_5y = "clagged_5y",
    extra_pcs = "extra_pcs",
    scale_down = "scale_down",
    scale_up = "scale_up"
  )

  scenario_candidates <- list(
    main = c(
      "10pc_avertable_burden_by_region_canonical_weighted_lower_region_v2.csv",
      "10pc_avertable_burden_by_region_canonical_weighted_upper_region_main_overall.csv"
    ),
    clagged = c(
      "10pc_avertable_burden_by_region_canonical_weighted_lower_region_v2_clagged.csv",
      "10pc_avertable_burden_by_region_canonical_weighted_lower_region_v2_consumption_lagged.csv"
    ),
    clagged_3y = c(
      "10pc_avertable_burden_by_region_canonical_weighted_lower_region_v2_clagged_3y.csv",
      "10pc_avertable_burden_by_region_canonical_weighted_lower_region_v2_consumption_lagged_3y.csv"
    ),
    clagged_5y = c(
      "10pc_avertable_burden_by_region_canonical_weighted_lower_region_v2_clagged_5y.csv",
      "10pc_avertable_burden_by_region_canonical_weighted_lower_region_v2_consumption_lagged_5y.csv"
    ),
    extra_pcs = c("10pc_avertable_burden_by_region_canonical_weighted_lower_region_v2_extra_pcs.csv"),
    scale_down = c("10pc_avertable_burden_by_region_canonical_weighted_lower_region_v2_scale_down.csv"),
    scale_up = c("10pc_avertable_burden_by_region_canonical_weighted_lower_region_v2_scale_up.csv")
  )

  per_scenario <- list()
  per_scenario_totals <- list()
  missing_scenarios <- character(0)

  for (scenario_id in scenario_order) {
    scenario_path <- find_existing_file(scenario_candidates[[scenario_id]], search_dirs)
    if (is.null(scenario_path)) {
      missing_scenarios <- c(missing_scenarios, scenario_id)
      per_scenario[[scenario_id]] <- NULL
      next
    }

    scenario_df <- read.csv(scenario_path, stringsAsFactors = FALSE, check.names = FALSE)
    required_cols <- c("region", "avertable_burden", "lower_bound", "upper_bound")
    if (!all(required_cols %in% names(scenario_df))) {
      warning(
        "Skipping avertable scenario '", scenario_id,
        "' because required columns are missing in ",
        scenario_path,
        call. = FALSE
      )
      per_scenario[[scenario_id]] <- NULL
      next
    }

    scenario_df$Avertable_10pc <- vapply(
      seq_len(nrow(scenario_df)),
      function(i) format_burden_ci(
        scenario_df$avertable_burden[[i]],
        scenario_df$lower_bound[[i]],
        scenario_df$upper_bound[[i]]
      ),
      character(1)
    )

    total_row <- data.frame(
      region = "Total",
      Avertable_10pc = format_burden_ci(
        sum(scenario_df$avertable_burden, na.rm = TRUE),
        sum(scenario_df$lower_bound, na.rm = TRUE),
        sum(scenario_df$upper_bound, na.rm = TRUE)
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    per_scenario[[scenario_id]] <- scenario_df[c("region", "Avertable_10pc")]
    per_scenario_totals[[scenario_id]] <- total_row
  }

  available <- Filter(Negate(is.null), per_scenario)
  if (length(available) == 0) {
    warning(
      "No 10pc avertable burden by region files found in search directories: ",
      paste(search_dirs, collapse = ", "),
      call. = FALSE
    )
    return(NULL)
  }

  regions <- unique(do.call(c, lapply(available, function(df) df$region)))
  out <- data.frame(region = regions, stringsAsFactors = FALSE, check.names = FALSE)

  for (scenario_id in scenario_order) {
    column_name <- scenario_labels[[scenario_id]]
    scenario_df <- per_scenario[[scenario_id]]
    if (is.null(scenario_df)) {
      out[[column_name]] <- ""
      next
    }

    merged <- merge(out["region"], scenario_df, by = "region", all.x = TRUE, sort = FALSE)
    vals <- merged$Avertable_10pc
    vals[is.na(vals)] <- ""
    out[[column_name]] <- vals
  }

  total_out <- data.frame(region = "Total", stringsAsFactors = FALSE, check.names = FALSE)
  for (scenario_id in scenario_order) {
    column_name <- scenario_labels[[scenario_id]]
    total_df <- per_scenario_totals[[scenario_id]]
    if (is.null(total_df)) {
      total_out[[column_name]] <- ""
    } else {
      total_out[[column_name]] <- total_df$Avertable_10pc[[1]]
    }
  }

  if (length(missing_scenarios) > 0) {
    warning(
      "Missing avertable burden scenarios for sensitivity table: ",
      paste(missing_scenarios, collapse = ", "),
      call. = FALSE
    )
  }

  out <- out[order(out$region), c("region", unname(unlist(scenario_labels[scenario_order])))]
  rbind(out, total_out)
}

bind_rows_fill <- function(dfs) {
  dfs <- Filter(Negate(is.null), dfs)
  if (length(dfs) == 0) {
    return(NULL)
  }

  all_names <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  aligned <- lapply(dfs, function(df) {
    missing_names <- setdiff(all_names, names(df))
    if (length(missing_names) > 0) {
      for (name in missing_names) {
        df[[name]] <- NA
      }
    }
    df[all_names]
  })

  do.call(rbind, aligned)
}

format_numeric_columns <- function(df, sig_digits = 4) {
  if (is.null(df)) {
    return(df)
  }

  for (col_name in names(df)) {
    if (is.numeric(df[[col_name]])) {
      rounded <- signif(df[[col_name]], sig_digits)
      df[[col_name]] <- vapply(rounded, function(value) {
        if (is.na(value)) {
          return("")
        }

        format(
          value,
          digits = sig_digits,
          trim = TRUE,
          scientific = NA,
          drop0trailing = TRUE
        )
      }, character(1))
    }
  }

  df
}

make_long_table <- function(df, value_cols, id_cols, prefix_pattern, value_name) {
  if (length(value_cols) == 0) {
    return(NULL)
  }

  rows <- lapply(value_cols, function(col_name) {
    variable_name <- sub(prefix_pattern, "", col_name)
    out <- df[id_cols]
    out$Variable <- variable_name
    out[[value_name]] <- df[[col_name]]
    out
  })

  bind_rows_fill(rows)
}

make_wide_metric_table <- function(df, value_cols, id_cols, prefix_pattern) {
  if (length(value_cols) == 0) {
    return(NULL)
  }

  out <- df[id_cols]
  for (col_name in value_cols) {
    variable_name <- sub(prefix_pattern, "", col_name)
    out[[variable_name]] <- df[[col_name]]
  }

  out
}

drop_columns <- function(df, columns) {
  if (is.null(df)) {
    return(df)
  }

  keep_names <- setdiff(names(df), columns)
  df[keep_names]
}

rename_atc_classes <- function(df) {
  if (is.null(df) || !("Antibiotic" %in% names(df))) {
    return(df)
  }

  mapped <- df$Antibiotic
  known_codes <- mapped %in% names(atc_names)
  mapped[known_codes] <- unname(unlist(atc_names[mapped[known_codes]], use.names = FALSE))
  df$Antibiotic <- mapped
  df
}

add_response_fmt <- function(df) {
  if (is.null(df)) {
    return(df)
  }

  has_response_cols <- all(c("Response", "Lower_CI", "Upper_CI") %in% names(df))
  if (!has_response_cols) {
    return(df)
  }

  response_fmt <- sprintf("%.2f (%.2f to %.2f)", df$Response, df$Lower_CI, df$Upper_CI)
  insert_at <- match("Upper_CI", names(df))

  df$Response <- NULL
  df$Lower_CI <- NULL
  df$Upper_CI <- NULL

  if ("Pathogen" %in% names(df)) {
    pathogen_index <- match("Pathogen", names(df))
    left_names <- names(df)[seq_len(pathogen_index)]
    right_names <- setdiff(names(df), left_names)
    out <- c(df[left_names], list(Response_fmt = response_fmt), df[right_names])
    return(as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE))
  }

  if (is.na(insert_at) || insert_at <= 0) {
    df$Response_fmt <- response_fmt
    return(df)
  }

  left_names <- names(df)[seq_len(min(insert_at - 3, ncol(df)))]
  right_names <- setdiff(names(df), left_names)
  out <- c(df[left_names], list(Response_fmt = response_fmt), df[right_names])
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

diagnostic_name_map <- c(
  R_squared = "R_squared",
  RMSE = "RMSE",
  AIC = "AIC",
  BIC = "BIC",
  White_test_p_value = "White p-value",
  `White p-value` = "White p-value"
)

diagnostic_tables <- list()
variation_tables <- list()
vif_tables <- list()
response_sensitivity_inputs <- list()
response_scenarios <- vapply(input_paths, detect_model_scenario, character(1))
if (length(input_paths) == 1 && is.na(response_scenarios[[1]])) {
  response_scenarios[[1]] <- "main"
}

for (input_idx in seq_along(input_paths)) {
  input_path <- input_paths[[input_idx]]
  if (!file.exists(input_path)) {
    stop("Input file not found: ", input_path, call. = FALSE)
  }

  df <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)
  source_name <- tools::file_path_sans_ext(basename(input_path))

  scenario_id <- response_scenarios[[input_idx]]

  diagnostic_cols <- intersect(names(diagnostic_name_map), names(df))
  variation_cols <- grep("^Vari(?:ance|ation)_Explained\\.", names(df), value = TRUE)
  vif_cols <- grep("^VIF\\.", names(df), value = TRUE)
  singularity_cols <- intersect(
    c("Model_Singular", "Bootstrap_Singular_Proportion"),
    names(df)
  )
  metric_cols <- unique(c(diagnostic_cols, variation_cols, vif_cols, singularity_cols))
  id_cols <- setdiff(names(df), metric_cols)

  diagnostics_df <- df[c(id_cols, diagnostic_cols)]
  if (length(diagnostic_cols) > 0) {
    renamed <- names(diagnostics_df)
    renamed[match(diagnostic_cols, names(diagnostics_df))] <- unname(diagnostic_name_map[diagnostic_cols])
    names(diagnostics_df) <- renamed
  }

  variation_df <- make_wide_metric_table(
    df,
    value_cols = variation_cols,
    id_cols = id_cols,
    prefix_pattern = "^Vari(?:ance|ation)_Explained\\."
  )

  vif_df <- make_wide_metric_table(
    df,
    value_cols = vif_cols,
    id_cols = id_cols,
    prefix_pattern = "^VIF\\."
  )

  diagnostic_tables[[length(diagnostic_tables) + 1]] <- diagnostics_df
  variation_tables[[length(variation_tables) + 1]] <- variation_df
  vif_tables[[length(vif_tables) + 1]] <- vif_df

  has_response_cols <- all(c("Antibiotic", "Pathogen", "Response", "Lower_CI", "Upper_CI") %in% names(df))
  if (!is.na(scenario_id) && has_response_cols && !(scenario_id %in% names(response_sensitivity_inputs))) {
    response_df <- df[c("Antibiotic", "Pathogen", "Response", "Lower_CI", "Upper_CI")]
    response_df <- rename_atc_classes(response_df)
    response_df$Response_fmt <- vapply(
      seq_len(nrow(response_df)),
      function(i) format_response_ci(
        response_df$Response[[i]],
        response_df$Lower_CI[[i]],
        response_df$Upper_CI[[i]]
      ),
      character(1)
    )
    response_sensitivity_inputs[[scenario_id]] <- response_df[c("Antibiotic", "Pathogen", "Response_fmt")]
  }
}

diagnostics_out <- bind_rows_fill(diagnostic_tables)
variation_out <- bind_rows_fill(variation_tables)
vif_out <- bind_rows_fill(vif_tables)

diagnostics_out <- rename_atc_classes(diagnostics_out)
variation_out <- rename_atc_classes(variation_out)
vif_out <- rename_atc_classes(vif_out)

diagnostics_out <- add_response_fmt(diagnostics_out)
variation_out <- drop_columns(variation_out, c("Response", "Lower_CI", "Upper_CI"))
vif_out <- drop_columns(vif_out, c("Response", "Lower_CI", "Upper_CI"))

diagnostics_out <- format_numeric_columns(diagnostics_out)
variation_out <- format_numeric_columns(variation_out)
vif_out <- format_numeric_columns(vif_out)

if (is.null(output_dir) || identical(output_dir, "")) {
  output_dir <- dirname(input_paths[[1]])
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

diagnostics_path <- file.path(output_dir, "model_fit_diagnostics.csv")
variation_path <- file.path(output_dir, "model_fit_variance_explained.csv")
vif_path <- file.path(output_dir, "model_fit_vif.csv")
response_sensitivity_path <- file.path(output_dir, "model_fit_response_sensitivity.csv")
avertable_region_sensitivity_path <- file.path(output_dir, "10pc_avertable_burden_by_region_sensitivity.csv")

search_dirs <- unique(c(
  normalizePath(output_dir, mustWork = FALSE),
  normalizePath(dirname(input_paths[[1]]), mustWork = FALSE),
  normalizePath(file.path(script_dir, "..", "Outputs"), mustWork = FALSE)
))

response_scenario_order <- c("main", "clagged", "clagged_3y", "clagged_5y", "extra_pcs")
response_scenario_labels <- c(
  main = "Main",
  clagged = "clagged",
  clagged_3y = "clagged_3y",
  clagged_5y = "clagged_5y",
  extra_pcs = "extra_pcs"
)

response_scenario_candidates <- list(
  main = c(
    "database_gradients_pathogen_ATC3_PCA_canonical_weighted_main.csv",
    "model_fit.csv"
  ),
  clagged = c(
    "database_gradients_pathogen_ATC3_PCA_canonical_weighted_clagged.csv",
    "database_gradients_pathogen_ATC3_PCA_canonical_weighted_consumption_lagged.csv"
  ),
  clagged_3y = c(
    "database_gradients_pathogen_ATC3_PCA_canonical_weighted_clagged_3y.csv",
    "database_gradients_pathogen_ATC3_PCA_canonical_weighted_consumption_lagged_3y.csv"
  ),
  clagged_5y = c(
    "database_gradients_pathogen_ATC3_PCA_canonical_weighted_clagged_5y.csv",
    "database_gradients_pathogen_ATC3_PCA_canonical_weighted_consumption_lagged_5y.csv"
  ),
  extra_pcs = c("database_gradients_pathogen_ATC3_PCA_canonical_weighted_extra_pcs.csv")
)

for (scenario_id in response_scenario_order) {
  if (scenario_id %in% names(response_sensitivity_inputs)) {
    next
  }

  scenario_path <- find_existing_file(response_scenario_candidates[[scenario_id]], search_dirs)
  response_from_file <- build_response_table_from_file(scenario_path)
  if (!is.null(response_from_file)) {
    response_sensitivity_inputs[[scenario_id]] <- response_from_file
  }
}

response_sensitivity_out <- build_response_sensitivity_table(
  response_sensitivity_inputs,
  response_scenario_order,
  response_scenario_labels
)

avertable_region_sensitivity_out <- build_avertable_region_sensitivity_table(search_dirs)

write.table(diagnostics_out, diagnostics_path, row.names = FALSE, quote = FALSE, sep = ",", na = "")
write.table(variation_out, variation_path, row.names = FALSE, quote = FALSE, sep = ",", na = "")
write.table(vif_out, vif_path, row.names = FALSE, quote = FALSE, sep = ",", na = "")
if (!is.null(response_sensitivity_out)) {
  write.table(response_sensitivity_out, response_sensitivity_path, row.names = FALSE, quote = FALSE, sep = ",", na = "")
}
if (!is.null(avertable_region_sensitivity_out)) {
  write.table(avertable_region_sensitivity_out, avertable_region_sensitivity_path, row.names = FALSE, quote = FALSE, sep = ",", na = "")
}

message("Wrote diagnostics table: ", diagnostics_path)
message("Wrote variance explained table: ", variation_path)
message("Wrote VIF table: ", vif_path)
if (!is.null(response_sensitivity_out)) {
  message("Wrote response sensitivity table: ", response_sensitivity_path)
} else {
  message("Skipped response sensitivity table (no eligible scenario-tagged inputs).")
}
if (!is.null(avertable_region_sensitivity_out)) {
  message("Wrote 10pc avertable burden by region sensitivity table: ", avertable_region_sensitivity_path)
} else {
  message("Skipped 10pc avertable burden by region sensitivity table (no eligible burden files found).")
}