rm(list = ls())
gc()

library(data.table)
library(forestplot)
library(magrittr)
library(dplyr)
library(tidyr)
library(tidyverse)

gradients <- fread('Outputs/database_gradients_pathogen_ATC3_PCA_joelike_weighted_test.csv')

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
# View in RStudio
View(final_table)

# Print to console
print(final_table)

# Export to a Tab-Separated file (best for pasting into Word/Google Docs)
write_tsv(final_table, "formatted_variable_explained_table_for_word.txt")

# library(ggplot2)

# p <- ggplot(df, aes(x = Pathogen, y = Variation_Explained, fill = Variable)) +
#     geom_bar(stat = "identity", position = "stack") +
#     scale_fill_manual(values = c(
#                                  "Variation_Explained.Interactions" = "#C0C0C0",
#                                  "Variation_Explained.Year" = "#648FFF",
#                                  "Variation_Explained.GDP" = "#785EF0", 
#                                  "Variation_Explained.PC3" = "#DC267F",
#                                  "Variation_Explained.PC2" = "#FF832B",
#                                  "Variation_Explained.PC1" = "#FFB000",
#                                  "Variation_Explained.Consumption" = "#000000"),
#                       labels = c("Variation_Explained.Year" = "Year",
#                                  "Variation_Explained.GDP" = "GDP", 
#                                  "Variation_Explained.PC3" = "PC3",
#                                  "Variation_Explained.PC2" = "PC2",
#                                  "Variation_Explained.PC1" = "PC1",
#                                  "Variation_Explained.Consumption" = "Consumption",
#                                  "Variation_Explained.Interactions" = "Interactions"),
#                       guide = guide_legend(reverse = TRUE)) +
#     scale_x_discrete(limits = rev) +
#     facet_wrap(~ factor(Antibiotic, levels = c("J01M", "J01G", "J01D", "J01C", "J01F", "J01E", "J01A"),
#                        labels = c("J01M" = "Quinolones", "J01G" = "Aminoglycosides", "J01D" = "Non-Penicillin Beta-Lactams", 
#                                  "J01C" = "Penicillins", "J01F" = "Macrolides", "J01E" = "Sulfonamides and Trimethoprim", 
#                                  "J01A" = "Tetracyclines")), ncol = 1, scales = "free_y") +
#     labs(title = "Variation Explained by Each Variable by Pathogen for Each Antibiotic Class",
#          x = "Pathogen",
#          y = "Variation Explained") +
#     theme_minimal() +
#     coord_flip()

# png('variation_explained_by_variable_stacked.png', width = 800, height = 1200)
# print(p)
# dev.off()







# gradients <- fread('Outputs/database_gradients_pathogen_ATC3_PCA_joelike_weighted.csv')
# bootstraps <- fread('Outputs/database_gradients_bootstraps_pathogen_ATC3_PCA_joelike_weighted.csv')

# bootstraps_permuted <- data.frame()
# for (ab in c("J01A", "J01C", "J01D", "J01E", "J01F", "J01G", "J01M")) {
#     bootstraps_permuted_ab <- fread(paste0('Outputs/database_gradients_bootstraps_pathogen_ATC3_PCA_joelike_weighted_permutation', ab, '.csv'))
#     bootstraps_permuted_ab$Antibiotic_permuted <- ab
#     bootstraps_permuted <- rbind(bootstraps_permuted, bootstraps_permuted_ab)
# }
# gradients_permuted <- data.frame()
# for (ab in c("J01A", "J01C", "J01D", "J01E", "J01F", "J01G", "J01M")) {
#     gradients_permuted_ab <- fread(paste0('Outputs/database_gradients_pathogen_ATC3_PCA_joelike_weighted_permutation', ab, '.csv'))
#     gradients_permuted_ab$Antibiotic_permuted <- ab
#     gradients_permuted <- rbind(gradients_permuted, gradients_permuted_ab)
# }

# # Wilcoxon signed rank test to compare distributions of bootstrapped gradients vs permuted gradients, per antibiotic class and pathogen
# results <- data.frame(Antibiotic = character(),
#                       Pathogen = character(),
#                       p_value = numeric(),
#                       stringsAsFactors = FALSE)

# pathogens <- unique(bootstraps$Pathogen)
# for (ab in c("J01A", "J01C", "J01D", "J01E", "J01F", "J01G", "J01M")) {
#     for (pathogen in pathogens) {
#         bootstrapped_values <- bootstraps %>%
#             filter(Antibiotic == ab, Pathogen == pathogen) %>%
#             pull(Gradient.Consumption)
#         permuted_values <- bootstraps_permuted %>%
#             filter(Antibiotic_permuted == ab, Pathogen == pathogen) %>%
#             pull(Gradient.Consumption)
        
#         if (length(bootstrapped_values) > 0 && length(permuted_values) > 0) {
#             test_result <- wilcox.test(bootstrapped_values, permuted_values, alternative = "greater")
#             results <- rbind(results, data.frame(Antibiotic = ab,
#                                                  Pathogen = pathogen,
#                                                  p_value = test_result$p.value,
#                                                  wilcoxon_statistic = test_result$statistic))
#         }
#     }
# }
# # save wilcoxon results
# fwrite(results, 'Outputs/wilcoxon_bootstrapped_vs_permuted_gradients.csv')

# # Add significance levels
# results <- results %>%
#     mutate(significance = case_when(
#         p_value < 0.001 ~ "X",
#         TRUE ~ ""
#     ))

# # # plot wilcoxon statistics
# # library(ggplot2)

# # p_wilcox <- ggplot(results, aes(x = Pathogen, y = wilcoxon_statistic)) +
# #     geom_bar(stat = "identity", fill = "steelblue") +
# #     geom_text(aes(label = significance, y = wilcoxon_statistic + max(wilcoxon_statistic) * 0.02), 
# #               size = 3, hjust = -0.1) +
# #     theme_minimal() +
# #     labs(x = "Pathogen",
# #          y = "Wilcoxon Statistic") +
# #     scale_x_discrete(limits = rev) +
# #     facet_wrap(~ factor(Antibiotic, levels = c("J01M", "J01G", "J01D", "J01C", "J01F", "J01E", "J01A"),
# #                        labels = c("J01M" = "Quinolones", "J01G" = "Aminoglycosides", "J01D" = "Non-Penicillin Beta-Lactams", 
# #                                   "J01C" = "Penicillins", "J01F" = "Macrolides", "J01E" = "Sulfonamides and Trimethoprim", 
# #                                   "J01A" = "Tetracyclines")), ncol = 1, scales = "free_y") +
# #     coord_flip()
# # png('wilcoxon_statistics_bootstrapped_vs_permuted_gradients_all_antibiotics.png', width = 1950, height = 2800, res = 300)
# # print(p_wilcox)
# # dev.off()

# # # plot distributions of bootstrapped vs permuted gradients for all antibiotic classes in one figure
# # library(ggplot2)
# # pathogens <- unique(bootstraps$Pathogen)

# # # Create mapping for antibiotic class names
# ab_labels <- c("J01M" = "Quinolones", "J01G" = "Aminoglycosides", "J01D" = "Non-Penicillin Beta-Lactams", 
#                "J01C" = "Penicillins", "J01F" = "Macrolides", "J01E" = "Sulfonamides and Trimethoprim", 
#                "J01A" = "Tetracyclines")

# # bootstraps$Type <- "Bootstrapped"
# # bootstraps_permuted$Type <- "Permuted"
# # # remove antibiotic_permuted==antibiotic
# # bootstraps_permuted <- bootstraps_permuted %>%
# #     filter(Antibiotic != Antibiotic_permuted)

# # combined_data <- rbind(
# #     bootstraps %>% filter(Gradient.Consumption >= -1.2, Gradient.Consumption <= 3.1) %>% select(Antibiotic, Gradient.Consumption, Type, Pathogen),
# #     bootstraps_permuted %>% filter(Gradient.Consumption >= -1.2, Gradient.Consumption <= 3.1) %>% select(Antibiotic = Antibiotic_permuted, Gradient.Consumption, Type, Pathogen)
# # )

# # # only include pathogens present in both datasets for each antibiotic
# # pathogens_in_both <- combined_data %>%
# #     group_by(Antibiotic, Pathogen) %>%
# #     summarise(has_both = length(unique(Type)) == 2, .groups = 'drop') %>%
# #     filter(has_both) %>%
# #     select(Antibiotic, Pathogen)

# # combined_data <- combined_data %>%
# #     inner_join(pathogens_in_both, by = c("Antibiotic", "Pathogen"))

# # # Create the combined plot

# # p <- ggplot(combined_data, aes(x = Pathogen, y = Gradient.Consumption, fill = Type)) +
# #     geom_violin(alpha = 0.5, position = position_dodge(0.9), width = 4) +
# #     geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
# #     theme_minimal() +
# #     labs(title = "Distribution of Bootstrapped vs Permuted Gradients by Antibiotic Class",
# #         x = "Pathogen",
# #         y = "Gradient Consumption") +
# #     scale_fill_manual(values = c("Bootstrapped" = "blue", "Permuted" = "red")) +
# #     scale_x_discrete(limits = rev) +
# #     facet_wrap(~ factor(Antibiotic, levels = c("J01M", "J01G", "J01D", "J01C", "J01F", "J01E", "J01A"),
# #                        labels = ab_labels[c("J01M", "J01G", "J01D", "J01C", "J01F", "J01E", "J01A")]), 
# #                ncol = 1, scales = "free_y") +
# #     coord_flip()

# # png('bootstrapped_vs_permuted_gradients_all_antibiotics_violin.png', width = 1200, height = 1600)
# # print(p)
# # dev.off()


# # calculate the distribution of the difference between bootstrapped and permuted gradients
# differences <- data.frame()
# for (ab in c("J01A", "J01C", "J01D", "J01E", "J01F", "J01G", "J01M")) {
#     print(paste("Processing antibiotic class:", ab))
#     for (pathogen in pathogens) {
#         bootstrapped_values <- bootstraps %>%
#             filter(Antibiotic == ab, Pathogen == pathogen) %>%
#             pull(Gradient.Consumption)
#         permuted_values <- bootstraps_permuted %>%
#             filter(Antibiotic_permuted == ab, Pathogen == pathogen) %>%
#             pull(Gradient.Consumption)
        
#         if (length(bootstrapped_values) > 0 && length(permuted_values) > 0) {
#             # Calculate 10000 independent differences (100x100)
#             n_samples <- min(50, length(bootstrapped_values), length(permuted_values))
            
#             print(paste("  Processing pathogen:", pathogen, "- Bootstrap samples:", length(bootstrapped_values), "Permuted samples:", length(permuted_values)))
            
#             # Sample independently for each difference calculation
#             for (i in 1:2500) {
#                 b <- sample(bootstrapped_values, 1)
#                 p <- sample(permuted_values, 1)
#                 differences <- rbind(differences, data.frame(Antibiotic = ab,
#                                                              Pathogen = pathogen,
#                                                              Difference = b - p))
#             }
#         }
#     }
# }
# # plot median and 95% CI of differences per antibiotic class and pathogen
# differences_summary <- differences %>%
#     group_by(Antibiotic, Pathogen) %>%
#     summarise(median_difference = median(Difference, na.rm = TRUE),
#               Lower_sd = median(Difference, na.rm = TRUE) - sd(Difference, na.rm = TRUE),
#               Upper_sd = median(Difference, na.rm = TRUE) + sd(Difference, na.rm = TRUE),
#               Lower_CI = quantile(Difference, 0.025, na.rm = TRUE),
#               Upper_CI = quantile(Difference, 0.975, na.rm = TRUE))
# print(differences_summary, n=100)
# # save
# fwrite(differences_summary, 'Outputs/bootstrapped_permuted_gradient_differences_sd_summary_downsampled.csv')

# # load
# differences_summary <- fread('Outputs/bootstrapped_permuted_gradient_differences_sd_summary_downsampled.csv')
# # print entire table
# print(differences_summary)

# library(ggplot2)
# # plot the differences
# p_diff <- ggplot(differences_summary, aes(x = Pathogen, y = median_difference)) +
#     geom_point() +
#     geom_errorbar(aes(ymin = Lower_sd, ymax = Upper_sd), width = 0.2) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
#     theme_minimal() +
#     labs(title = "Median Difference between Bootstrapped and Permuted Gradients by Antibiotic Class",
#         x = "Pathogen",
#         y = "Median Difference (Bootstrapped - Permuted)") +
#     scale_x_discrete(limits = rev) +
#     scale_y_continuous(limits = c(-10, 10)) +
#     facet_wrap(~ factor(Antibiotic, levels = c("J01M", "J01G", "J01D", "J01C", "J01F", "J01E", "J01A"),
#                        labels = ab_labels[c("J01M", "J01G", "J01D", "J01C", "J01F", "J01E", "J01A")]), 
#                ncol = 1) +
#     coord_flip()
# png('bootstrapped_permuted_gradients_all_antibiotics_median_difference_sd_downsampled.png', width = 1950, height = 2800, res = 300)
# print(p_diff)
# dev.off()








# # calculate standard deviation of gradients within drug classes, scaled by mean
# df_summary <- df %>%
#   group_by(Antibiotic) %>%
#   summarise(mean_gradient = mean(Response, na.rm = TRUE),
#             sd_gradient = sd(Response, na.rm = TRUE),
#             n = n()) %>%
#   mutate(cv_gradient = sd_gradient / abs(mean_gradient)) %>%
#   arrange(desc(cv_gradient))

# # same organised by pathogen
# df_summary_pathogen <- df %>%
#   group_by(Pathogen) %>%
#   summarise(mean_gradient = mean(Response, na.rm = TRUE),
#             sd_gradient = sd(Response, na.rm = TRUE),
#             n = n()) %>%
#   mutate(cv_gradient = sd_gradient / abs(mean_gradient)) %>%
#   arrange(desc(cv_gradient))

# # calculate standard deviation of gradients within drug classes, with bootstrapped confidence intervals
# bootstraps <- bootstraps %>%
#   group_by(Pathogen, Antibiotic) %>%
#   summarise(sd_gradient = sd(Gradient.Consumption, na.rm = TRUE))
# # bootstraps$sd_gradient <- sd()

# df_summary <- df %>%
#   group_by(Pathogen) %>%
#   summarise(mean_gradient = mean(Response, na.rm = TRUE),
#             sd_gradient = sd(Response, na.rm = TRUE),
#             n = n()) %>%
#   left_join(
#     bootstraps %>%
#       group_by(Pathogen) %>%
#       summarise(
#         sd_gradient_lower = quantile(sd_gradient, 0.025, na.rm = TRUE),
#         sd_gradient_upper = quantile(sd_gradient, 0.975, na.rm = TRUE)
#       ),
#     by = "Pathogen"
#   ) %>%
#   arrange(desc(cv_gradient))

# print(df_summary)
# print(mean(df_summary$sd_gradient, na.rm = TRUE))

# # calculate the mean standard deviation 1000 times, sampling from the bootstrapped sd gradients, to get a confidence interval on the mean sd gradient
# set.seed(42)
# mean_sds <- replicate(1000, {
#   sampled_sds <- bootstraps %>%
#     group_by(Pathogen, Antibiotic) %>%
#     summarise(sd_gradient = sample(Gradient.Consumption, 1), .groups = 'drop') %>%
#     ungroup() %>%
#     group_by(Antibiotic) %>%
#     summarise(mean_sd = mean(sd_gradient, na.rm = TRUE))
#   mean(sampled_sds$mean_sd, na.rm = TRUE)
# })
# mean_sd <- mean(mean_sds, na.rm = TRUE)
# mean_sd_lower <- quantile(mean_sds, 0.025, na.rm = TRUE)
# mean_sd_upper <- quantile(mean_sds, 0.975, na.rm = TRUE)
# print(paste0("Mean SD of gradients by antibiotic: ", round(mean_sd, 4), " (95% CI: ", round(mean_sd_lower, 4), "-", round(mean_sd_upper, 4), ")"))
# mean_sds <- replicate(1000, {
#   sampled_sds <- bootstraps %>%
#     group_by(Pathogen, Antibiotic) %>%
#     summarise(sd_gradient = sample(Gradient.Consumption, 1), .groups = 'drop') %>%
#     ungroup() %>%
#     group_by(Pathogen) %>%
#     summarise(mean_sd = mean(sd_gradient, na.rm = TRUE))
#   mean(sampled_sds$mean_sd, na.rm = TRUE)
# })
# mean_sd <- mean(mean_sds, na.rm = TRUE)
# mean_sd_lower <- quantile(mean_sds, 0.025, na.rm = TRUE)
# mean_sd_upper <- quantile(mean_sds, 0.975, na.rm = TRUE)
# print(paste0("Mean SD of gradients by pathogen: ", round(mean_sd, 4), " (95% CI: ", round(mean_sd_lower, 4), "-", round(mean_sd_upper, 4), ")"))