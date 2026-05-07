#!/usr/bin/env Rscript

library(dplyr)

input_path <- "Outputs/wilcoxon_bootstrapped_vs_permuted_gradients.csv"
output_path <- "Outputs/supplementary_table_2_wilcoxon.csv"

if (!file.exists(input_path)) {
  message(sprintf("[table2] Skipping: required input not found: %s", input_path))
  message("[table2]   Run the Wilcoxon test stage (or provide permuted bootstrap CSVs) first.")
  if (sys.nframe() == 0L) quit(save = "no", status = 0)
  return(invisible(NULL))
}

df <- read.csv(input_path, stringsAsFactors = FALSE)

required_cols <- c("Antibiotic", "Pathogen", "p_value")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop(sprintf("[table2] Missing required columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
}

if (!"wilcoxon_statistic" %in% names(df)) {
  df$wilcoxon_statistic <- NA_real_
}

out <- df %>%
  transmute(
    Antibiotic,
    Pathogen,
    p_value = as.numeric(p_value),
    wilcoxon_statistic = as.numeric(wilcoxon_statistic),
    significant_p_lt_0_001 = ifelse(as.numeric(p_value) < 0.001, "yes", "no")
  ) %>%
  arrange(Antibiotic, Pathogen)

write.csv(out, output_path, row.names = FALSE)
message("[table2] Wrote: ", output_path)
