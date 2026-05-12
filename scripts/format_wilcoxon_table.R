#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

input_path <- getOption(
  "amr_wilcoxon_table_input_path",
  "Outputs/wilcoxon_bootstrapped_vs_permuted_gradients.csv"
)
output_path <- getOption(
  "amr_wilcoxon_table_output_path",
  "Outputs/wilcoxon_bootstrapped_vs_permuted_gradients_table.txt"
)

antibiotic_map <- c(
  "J01A" = "Tetracyclines",
  "J01C" = "Penicillins",
  "J01D" = "Non-Penicillin Beta-Lactams",
  "J01E" = "Sulfonamides and Trimethoprim",
  "J01F" = "Macrolides",
  "J01G" = "Aminoglycosides",
  "J01M" = "Quinolones"
)

class_order <- c(
  "Tetracyclines",
  "Penicillins",
  "Non-Penicillin Beta-Lactams",
  "Sulfonamides and Trimethoprim",
  "Macrolides",
  "Aminoglycosides",
  "Quinolones"
)

pathogen_order <- c(
  "Acinetobacter spp.",
  "E. coli",
  "E. faecalis",
  "E. faecium",
  "Enterococcus spp.",
  "H. influenzae",
  "K. pneumoniae",
  "Morganella spp.",
  "N. gonorrhoeae",
  "P. aeruginosa",
  "S. agalactiae",
  "S. aureus",
  "S. pneumoniae",
  "S. pyogenes",
  "Salmonella spp."
)

format_p_value <- function(x) {
  out <- rep("—", length(x))
  non_missing <- !is.na(x)
  out[non_missing & x == 0] <- "0"
  out[non_missing & x == 1] <- "1"

  needs_scientific <- non_missing & x != 0 & x != 1
  if (any(needs_scientific)) {
    sci <- formatC(x[needs_scientific], format = "e", digits = 1)
    sci_parts <- strsplit(sci, "e", fixed = TRUE)
    sci <- vapply(sci_parts, function(parts) {
      exponent <- as.integer(parts[[2]])
      if (identical(parts[[1]], "1.0") && identical(exponent, 0L)) {
        return("1")
      }
      paste0(parts[[1]], "×10", exponent)
    }, character(1))
    out[needs_scientific] <- sci
  }

  out
}

if (!file.exists(input_path)) {
  message("[wilcoxon-table] Skipping: input not found: ", input_path)
  if (sys.nframe() == 0L) quit(save = "no", status = 0)
  stop("[wilcoxon-table] Input file missing", call. = FALSE)
}

df <- fread(input_path)

required_cols <- c("Antibiotic", "Pathogen", "p_value")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop(
    "[wilcoxon-table] Missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

pathogens_in_data <- unique(df$Pathogen)
pathogen_order <- c(pathogen_order, setdiff(sort(pathogens_in_data), pathogen_order))

table_df <- df %>%
  mutate(
    Antibiotic = unname(antibiotic_map[Antibiotic]),
    Pathogen = as.character(Pathogen),
    p_value = as.numeric(p_value)
  ) %>%
  filter(!is.na(Antibiotic)) %>%
  mutate(
    Antibiotic = factor(Antibiotic, levels = class_order),
    Pathogen = factor(Pathogen, levels = pathogen_order)
  ) %>%
  select(Pathogen, Antibiotic, p_value) %>%
  complete(Pathogen, Antibiotic = class_order) %>%
  mutate(p_value = format_p_value(p_value)) %>%
  pivot_wider(names_from = Antibiotic, values_from = p_value) %>%
  mutate(Pathogen = as.character(Pathogen)) %>%
  arrange(match(Pathogen, pathogen_order)) %>%
  select(Pathogen, all_of(class_order))

fwrite(table_df, output_path, sep = "\t", quote = FALSE, na = "—")
message("[wilcoxon-table] Wrote: ", output_path)