rm(list = ls())
gc()

library(data.table)
library(forestplot)
library(magrittr)
library(dplyr)
library(tidyr)
library(tidyverse)

gradients <- fread(getOption("amr_table3_gradients_path",
    "Outputs/database_gradients_pathogen_ATC3_PCA_canonical_weighted_main.csv"))

# Calculate interactions as the difference between R_squared and sum of individual variations
df <- gradients %>%
    mutate(Variation_Explained.Interactions = R_squared - (Variation_Explained.Year + Variation_Explained.GDP + 
                                                          Variation_Explained.PC3 + Variation_Explained.PC2 + 
                                                          Variation_Explained.PC1 + Variation_Explained.Consumption)) %>%
    select(Pathogen, Antibiotic, R_squared,Variation_Explained.Interactions,
           Variation_Explained.Year, Variation_Explained.GDP, Variation_Explained.PC3, 
           Variation_Explained.PC2, Variation_Explained.PC1, Variation_Explained.Consumption,
           ) %>%
    pivot_longer(cols = c(Variation_Explained.Interactions,Variation_Explained.Year, Variation_Explained.GDP, Variation_Explained.PC3, 
                          Variation_Explained.PC2, Variation_Explained.PC1, Variation_Explained.Consumption,
                          ), 
                 names_to = "Variable", values_to = "Variation_Explained") %>%
    mutate(Variable = factor(Variable, levels = c("Variation_Explained.Interactions", "Variation_Explained.Year", "Variation_Explained.GDP", "Variation_Explained.PC3", 
                                                  "Variation_Explained.PC2", "Variation_Explained.PC1", "Variation_Explained.Consumption"
                                                  )))

data <- read.csv("merged_data_N_PC3_GDP.csv")
data <- data %>%
    rename(
        Antibiotic = ATC.Class,
        Consumption = Antibiotic.Consumption,
        Resistance = Percent.Resistant.Isolates,
        Pathogen = Pathogen,
        Weight = Total.Isolates
        )
# to find the total variation explained by Consumption, weight each pathogen-antibiotic combination by the number of isolates (sum of Total.Isolates in data)
data_counts <- data %>%
    group_by(Pathogen, Antibiotic) %>%
    summarise(Weight = sum(Weight, na.rm = TRUE))
df <- df %>%
    left_join(data_counts, by = c("Pathogen", "Antibiotic")) %>%
    mutate(Weighted_Variation_Explained = Variation_Explained * Weight)
# calculate total weighted variation explained by Consumption across all pathogen-antibiotic combinations
total_weighted_variation_explained_consumption <- sum(df %>% filter(Variable == "Variation_Explained.Consumption") %>% pull(Weighted_Variation_Explained), na.rm = TRUE)
total_isolates <- sum(df$Weight, na.rm = TRUE)
overall_variation_explained_consumption <- total_weighted_variation_explained_consumption / total_isolates
print(paste0("Overall variation explained by Consumption: ", round(overall_variation_explained_consumption, 4)))
# and for each other variable
for (var in c("Variation_Explained.Year", "Variation_Explained.GDP", "Variation_Explained.PC3", 
              "Variation_Explained.PC2", "Variation_Explained.PC1", "Variation_Explained.Interactions")) {
    total_weighted_variation_explained_var <- sum(df %>% filter(Variable == var) %>% pull(Weighted_Variation_Explained), na.rm = TRUE)
    overall_variation_explained_var <- total_weighted_variation_explained_var / total_isolates
    print(paste0("Overall variation explained by ", var, ": ", round(overall_variation_explained_var, 4)))
}

# 2. Define the mapping from Antibiotic codes to Class names
antibiotic_map <- c(
  "J01M" = "Quinolones",
  "J01G" = "Aminoglycosides",
  "J01D" = "Non-Penicillin Beta-Lactams",
  "J01C" = "Penicillins",
  "J01F" = "Macrolides",
  "J01E" = "Sulfonamides and Trimethoprim",
  "J01A" = "Tetracyclines"
)

# 3. Define the specific sort order for the Antibiotic Classes
class_order <- c(
  "Quinolones",
  "Aminoglycosides",
  "Non-Penicillin Beta-Lactams",
  "Penicillins",
  "Macrolides",
  "Sulfonamides and Trimethoprim",
  "Tetracyclines"
)

# 4. Process the data
final_table <- df %>%
  # Create the Antibiotic Class column using the map
  mutate(Antibiotic_Class = antibiotic_map[Antibiotic]) %>%
  # Clean up the Variable names (remove "Variation_Explained.")
  mutate(Variable = str_replace(Variable, "Variation_Explained\\.", "")) %>%
  # Select only the columns needed for the pivot
  select(Antibiotic_Class, Pathogen, Variable, Variation_Explained) %>%
  # Pivot from long to wide format
  pivot_wider(names_from = Variable, values_from = Variation_Explained) %>%
  # Reorder the columns to match your desired format
  select(Antibiotic_Class, Pathogen, Consumption, PC1, PC2, PC3, GDP, Year, Interactions) %>%
  # Convert Antibiotic_Class to a factor to enforce the specific sort order
  mutate(Antibiotic_Class = factor(Antibiotic_Class, levels = class_order)) %>%
  # Sort the rows
  arrange(Antibiotic_Class, Pathogen) %>%
  # Round all numeric columns to 3 decimal places
  mutate(across(where(is.numeric), ~ round(., 3)))

# 5. Output the result
# Print to console
print(final_table)

# Export to a Tab-Separated file (best for pasting into Word/Google Docs)
write_tsv(final_table, "formatted_variable_explained_table_for_word.txt")
