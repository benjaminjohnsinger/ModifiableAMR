#!/usr/bin/env Rscript

source("config.R")

scenario <- get_amr_scenario()

expected_merged_path <- switch(
  scenario,
  main = "merged_data_N_PC3_GDP.csv",
  hic = "merged_data_N_PC3_GDP.csv",
  lmic = "merged_data_N_PC3_GDP.csv",
  raw_iqvia = "merged_data_N_PC3_GDP_IQVIA.csv",
  hospital_nagorsen = "merged_data_Nagorsen_hospital_to_all_filtered.csv",
  NA_character_
)

if (is.na(expected_merged_path)) {
  message("[validate-stage-data] No merged-data checks defined for scenario: ", scenario)
  quit(status = 0)
}

if (!file.exists(expected_merged_path)) {
  message("[validate-stage-data] Missing merged dataset: ", expected_merged_path)
  quit(status = 1)
}

df <- read.csv(expected_merged_path, stringsAsFactors = FALSE)

required_cols <- c(
  "ISO3", "Year", "Pathogen", "ATC.Class", "Total.Isolates",
  "Percent.Resistant.Isolates", "Antibiotic.Consumption"
)

if (scenario %in% c("main", "hic", "lmic", "raw_iqvia")) {
  required_cols <- c(required_cols, "PC3", "GDP")
}

missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  message("[validate-stage-data] Missing required columns in ", expected_merged_path, ":")
  for (col in missing_cols) message("  - ", col)
  quit(status = 1)
}

key_cols <- c("ISO3", "Year", "Pathogen", "ATC.Class")
if (all(key_cols %in% names(df))) {
  key_counts <- as.data.frame(table(interaction(df[key_cols], drop = TRUE)), stringsAsFactors = FALSE)
  repeated_keys <- sum(key_counts$Freq > 1)
  if (repeated_keys > 0) {
    message("[validate-stage-data] Note: repeated analytic keys detected for ISO3-Year-Pathogen-ATC.Class: ", repeated_keys)
    message("[validate-stage-data] This is allowed at this stage when multiple observations exist within a key.")
  }
}

exact_dup_rows <- duplicated(df)
if (any(exact_dup_rows)) {
  message("[validate-stage-data] Note: exact duplicate rows detected: ", sum(exact_dup_rows))
  message("[validate-stage-data] Review upstream harmonization if this count grows unexpectedly.")
}

join_critical_cols <- intersect(c("Antibiotic.Consumption", "PC3", "GDP"), names(df))
join_missing <- vapply(join_critical_cols, function(col) sum(is.na(df[[col]])), integer(1))
if (any(join_missing > 0)) {
  message("[validate-stage-data] Missing values found in join-critical columns:")
  for (col in names(join_missing)) {
    if (join_missing[[col]] > 0) {
      message(sprintf("  - %s: %d missing", col, join_missing[[col]]))
    }
  }
  quit(status = 1)
}

numeric_checks <- intersect(c("Antibiotic.Consumption", "Percent.Resistant.Isolates", "PC3", "GDP"), names(df))
zero_var <- character(0)
for (col in numeric_checks) {
  x <- suppressWarnings(as.numeric(df[[col]]))
  x <- x[is.finite(x)]
  if (length(x) == 0 || sd(x) == 0) {
    zero_var <- c(zero_var, col)
  }
}

if (length(zero_var) > 0) {
  message("[validate-stage-data] Zero-variance numeric columns detected:")
  for (col in zero_var) message("  - ", col)
  quit(status = 1)
}

message("[validate-stage-data] Checks passed for scenario: ", scenario)
