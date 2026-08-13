library(dplyr)
library(readr)

apply_fdr_to_bootstraps <- function(results_path, bootstraps_path, output_path) {
  
  # 1. Load the data
  results <- read_csv(results_path, show_col_types = FALSE)
  bootstraps <- read_csv(bootstraps_path, show_col_types = FALSE)
  
  print(head(results))
  print(head(bootstraps))

  # 2. Calculate empirical p-values from bootstraps
  p_vals_df <- bootstraps %>%
    group_by(Antibiotic, Pathogen) %>%
    summarise(
      n_bootstraps = n(),
      
      # Detect fallback: if variance is 0, the script copied the point estimate
      is_fallback = var(Gradient, na.rm = TRUE) == 0,
      
      # Calculate empirical p-value
      prop_below = mean(Gradient <= 0, na.rm = TRUE),
      prop_above = mean(Gradient >= 0, na.rm = TRUE),
      
      # 2 * min(proportion_below, proportion_above)
      empirical_p = if_else(
        is_fallback, 
        NA_real_, # Exclude fallback models from this calculation
        2 * min(prop_below, prop_above)
      ),
      .groups = "drop"
    )
  
  # 3. Merge and apply Benjamini-Hochberg correction
  final_results <- results %>%
    left_join(p_vals_df %>% dplyr::select(Antibiotic, Pathogen, empirical_p, is_fallback), 
              by = c("Antibiotic", "Pathogen")) %>%
    mutate(
      # p.adjust ignores NAs by default and calculates FDR on valid tests
      BH_FDR = p.adjust(empirical_p, method = "BH"),
      is_significant = !is.na(BH_FDR) & BH_FDR < 0.05
    )
  
  # 4. Save and summarize
  write_csv(final_results, output_path)
  
  cat("Processed", nrow(final_results), "drug-bug combinations.\n")
  cat("Found", sum(final_results$is_fallback, na.rm = TRUE), "fallback models (excluded from FDR).\n")
  cat("Significant results (FDR < 0.05):", sum(final_results$is_significant, na.rm = TRUE), "\n")
  
  # print drug-bugs for which FDR is significant
  significant_results <- final_results %>%
    filter(is_significant) %>%
    dplyr::select(Antibiotic, Pathogen, Response, empirical_p, BH_FDR)
  print(significant_results)

  return(final_results)
}

# Example execution:
updated_results <- apply_fdr_to_bootstraps(
  results_path = "Outputs/database_gradients_pathogen_ATC3_PCA_canonical_weighted_main.csv",
  bootstraps_path = "Outputs/database_gradients_bootstraps_pathogen_ATC3_PCA_canonical_weighted_main.csv",
  output_path = "Outputs/database_gradients_FDR_corrected_main.csv"
)