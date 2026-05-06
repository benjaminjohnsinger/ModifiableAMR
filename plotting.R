library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggforce)
library(svglite)
library(ggtext)
library(scales)
source("utils.R")

# # main results
# results <- read.csv("Outputs/database_gradients_pathogen_ATC3_PCA_joelike60cutoff_no0resistance.csv", stringsAsFactors = FALSE)

classes <- c("J01A","J01B","J01C","J01D","J01E","J01F","J01G","J01M")
antibiotic_names<- c(  "Tetracyclines",
  "Glycopeptides and Lipopeptides",
  "Penicillins",
  "Non-Penicillin Beta-Lactams",
  "Sulfonamides and Trimethoprim",
  "Macrolides",
  "Aminoglycosides",
  "Quinolones") 
atc_names <- setNames(antibiotic_names, classes)

# # Figure 1 - all antibiotic-pathogen pairs with drug-specific pathogen random effects models
# drug_class_gradients <- setNames(as.vector(read.csv("Outputs/database_gradients_ATC3_PCA_joelike60cutoff_no0resistance.csv", stringsAsFactors = FALSE)[,2]), 
#                             read.csv("Outputs/database_gradients_ATC3_PCA_joelike60cutoff_no0resistance.csv", stringsAsFactors = FALSE)[,1])
# drug_class_lowerCI <- setNames(as.vector(read.csv("Outputs/database_lowerCI_ATC3_PCA_joelike60cutoff_no0resistance.csv", stringsAsFactors = FALSE)[,2]),
#                             read.csv("Outputs/database_lowerCI_ATC3_PCA_joelike60cutoff_no0resistance.csv", stringsAsFactors = FALSE)[,1])
# drug_class_upperCI <- setNames(as.vector(read.csv("Outputs/database_upperCI_ATC3_PCA_joelike60cutoff_no0resistance.csv", stringsAsFactors = FALSE)[,2]),
#                             read.csv("Outputs/database_upperCI_ATC3_PCA_joelike60cutoff_no0resistance.csv", stringsAsFactors = FALSE)[,1])
# antibiotics <- unique(results$Antibiotic)


# segment_data <- data.frame(
#     Antibiotic = antibiotics,
#     Response = drug_class_gradients[antibiotics],
#     Lower_CI = drug_class_lowerCI[antibiotics],
#     Upper_CI = drug_class_upperCI[antibiotics]
# )

# # create interaction of Pathogen and Antibiotic in results to replace those columns, then combine with segment_data
# # create an extra column that is 0 for data form results and 1 for data from segment_data
# # order by antibiotic, with segment_data at the top of each section
# results$Segment <- 0
# results <- results %>%
#     select(Antibiotic, Pathogen, Response, Lower_CI, Upper_CI, Segment)
# segment_data$Pathogen <- " Overall"
# segment_data$Segment <- 1
# segment_data <- segment_data %>%
#     select(Antibiotic, Pathogen, Response, Lower_CI, Upper_CI, Segment)
# # combine results and segment_data
# combined_results <- rbind(results, segment_data)

# pathogen_classes <- c(" Overall",
# "Acinetobacter spp.","E. coli","E. faecalis","E. faecium","Enterococcus spp.",
# "H. influenzae","K. pneumoniae","Morganella spp.","N. gonorrhoeae","P. aeruginosa",
# "S. agalactiae","S. aureus","S. pneumoniae","S. pyogenes","Salmonella spp.")
# pathogen_names <- c(" Overall",
# "*Acinetobacter* spp","*E. coli*","*E. faecalis*","*E. faecium*","*Enterococcus* spp",
# "*H. influenzae*","*K. pneumoniae*","*Morganella* spp","*N. gonorrhoeae*","*P. aeruginosa*",
# "*S. agalactiae*","*S. aureus*","*S. pneumoniae*","*S. pyogenes*","*Salmonella* spp")
# pathogen_names <- setNames(pathogen_names, pathogen_classes)
# combined_results$Pathogen <- as.character(combined_results$Pathogen)  # Convert factors to characters
# for (i in 1:nrow(combined_results)) {
#     combined_results$Pathogen[i] <- pathogen_names[[combined_results$Pathogen[i]]]
# }

# # replace Antibiotic names with antibiotic_names
# combined_results$Antibiotic <- as.character(combined_results$Antibiotic)  # Convert factors to characters
# for (i in 1:nrow(combined_results)) {
#     combined_results$Antibiotic[i] <- atc_names[[combined_results$Antibiotic[i]]]
# }

# # add column with formatted response with CI
# combined_results$Antibiotic <- factor(combined_results$Antibiotic, levels = rev(antibiotic_names))
# combined_results$Pathogen <- factor(combined_results$Pathogen, levels = rev(pathogen_names))
# combined_results$Response_fmt <- paste0(round(combined_results$Response, 2), " (", round(combined_results$Lower_CI, 2), "–", round(combined_results$Upper_CI, 2), ")")
# for (i in 1:nrow(combined_results)) {
#     combined_results$Response_fmt[i] <- paste0(combined_results$Pathogen[i], "...", paste(rep(".", floor(1.49*(18 - nchar(combined_results$Response_fmt[i])))),collapse=""), combined_results$Response_fmt[i])
# }
# # combined_results$Response_fmt <- paste0(combined_results$Pathogen, "...", rep(".",18-nchar(combined_results$Response_fmt)), combined_results$Response_fmt)
# combined_results$Response_fmt <- factor(combined_results$Response_fmt, levels = rev(unique(combined_results$Response_fmt)))

# # forest plot
# # y is drug and pathogen combos
# plot <- ggplot(combined_results, aes(x = Response, y = Response_fmt, color = as.factor(Segment))) +
#     geom_point(stat = "identity") +
#     geom_errorbar(aes(xmin = Lower_CI, xmax = Upper_CI), width = 0.2, position = position_dodge(0.9)) +
#     labs(x = "Elasticity of resistance with respect to consumption", y = "") +  # Added title
#     scale_color_manual(values = c("0" = "grey", "1" = "black"), guide = FALSE) +  # Disable the color legend
#     geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#     # add title
#     # ggtitle("Elasticity of Resistance with Respect to Consumption by Antibiotic and Pathogen") +
#     theme_minimal() +
#     theme(panel.background = element_rect(fill = "white")) +
#     theme(axis.title.x = element_text(size = 10)) +
#     theme(axis.title.y = element_text(size = 10)) +
#     theme(axis.text.x = element_text(size = 10)) +
#     theme(axis.text.y = element_markdown(size = 8)) +
#     # theme(plot.title = element_text(size = 20)) +  # Set title size
#     theme(strip.text = element_text(size = 10)) +  # Set facet titles size
#     ggforce::facet_col(Antibiotic ~ ., scales = "free_y", space = "free", drop = TRUE) # Create facets for each antibiotic in a single column
# ggsave("gradients_drugs_pathogen_PCA_joelike60cutoff_no0resistance_response_fmt.png", plot, width = 6.5, height = 9)

# # Figure 2 - pathogen-specific models with drug random effects
# pathogen_gradients <- setNames(as.vector(read.csv("Outputs/database_gradients_pathogen_PCA_joelike60cutoff_no0resistance.csv", stringsAsFactors = FALSE)[,2]), 
#                             read.csv("Outputs/database_gradients_pathogen_PCA_joelike60cutoff_no0resistance.csv", stringsAsFactors = FALSE)[,1])
# pathogen_lowerCI <- setNames(as.vector(read.csv("Outputs/database_lowerCI_pathogen_PCA_joelike60cutoff_no0resistance.csv", stringsAsFactors = FALSE)[,2]),
#                             read.csv("Outputs/database_lowerCI_pathogen_PCA_joelike60cutoff_no0resistance.csv", stringsAsFactors = FALSE)[,1])
# pathogen_upperCI <- setNames(as.vector(read.csv("Outputs/database_upperCI_pathogen_PCA_joelike60cutoff_no0resistance.csv", stringsAsFactors = FALSE)[,2]),
#                             read.csv("Outputs/database_upperCI_pathogen_PCA_joelike60cutoff_no0resistance.csv", stringsAsFactors = FALSE)[,1])
# pathogens <- unique(results$Pathogen)

# segment_data_pathogen <- data.frame(
#     Pathogen = pathogens,
#     Response = pathogen_gradients[pathogens],
#     Lower_CI = pathogen_lowerCI[pathogens],
#     Upper_CI = pathogen_upperCI[pathogens]
# )

# pathogen_classes <- c(" Overall",
# "Acinetobacter spp.","E. coli","E. faecalis","E. faecium","Enterococcus spp.",
# "H. influenzae","K. pneumoniae","Morganella spp.","N. gonorrhoeae","P. aeruginosa",
# "S. agalactiae","S. aureus","S. pneumoniae","S. pyogenes","Salmonella spp.")
# pathogen_names <- c(" Overall",
# "*Acinetobacter* spp","*E. coli*","*E. faecalis*","*E. faecium*","*Enterococcus* spp",
# "*H. influenzae*","*K. pneumoniae*","*Morganella* spp","*N. gonorrhoeae*","*P. aeruginosa*",
# "*S. agalactiae*","*S. aureus*","*S. pneumoniae*","*S. pyogenes*","*Salmonella* spp")
# pathogen_names <- setNames(pathogen_names, pathogen_classes)
# segment_data_pathogen$Pathogen <- as.character(segment_data_pathogen$Pathogen)  # Convert factors to characters
# for (i in 1:nrow(segment_data_pathogen)) {
#     segment_data_pathogen$Pathogen[i] <- pathogen_names[[segment_data_pathogen$Pathogen[i]]]
# }
# segment_data_pathogen$Pathogen <- factor(segment_data_pathogen$Pathogen, levels = rev(sort(unique(segment_data_pathogen$Pathogen))))

# segment_data_pathogen$Response_fmt <- paste0(round(segment_data_pathogen$Response, 2), " (", round(segment_data_pathogen$Lower_CI, 2), "–", round(segment_data_pathogen$Upper_CI, 2), ")")
# for (i in 1:nrow(segment_data_pathogen)) {
#     segment_data_pathogen$Response_fmt[i] <- paste0(segment_data_pathogen$Pathogen[i], "...", paste(rep(".", floor(1.49*(18 - nchar(segment_data_pathogen$Response_fmt[i])))),collapse=""), segment_data_pathogen$Response_fmt[i])
# }
# segment_data_pathogen$Response_fmt <- factor(segment_data_pathogen$Response_fmt, levels = rev(unique(segment_data_pathogen$Response_fmt)))

# plot <- ggplot(segment_data_pathogen, aes(x = Response, y = fct_reorder(Response_fmt, Response))) +
#     geom_point(stat = "identity") +
#     geom_errorbar(aes(xmin = Lower_CI, xmax = Upper_CI), width = 0.2, position = position_dodge(0.9)) +
#     labs(x = "Elasticity of resistance with respect to consumption", y = "") +  # Added title
#     scale_color_manual(values = c("0" = "grey", "1" = "black"), guide = FALSE) +  # Disable the color legend
#     geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#     # add title
#     # ggtitle("Elasticity of Resistance with Respect to Consumption by Antibiotic and Pathogen") +
#     theme_minimal() +
#     theme(panel.background = element_rect(fill = "white")) +
#     theme(axis.title.x = element_text(size = 10)) +
#     theme(axis.title.y = element_text(size = 10)) +
#     theme(axis.text.x = element_text(size = 10)) +
#     theme(axis.text.y = element_markdown(size = 8)) +
#     # theme(plot.title = element_text(size = 20)) +  # Set title size
#     theme(strip.text = element_text(size = 10)) # Set facet titles size
# ggsave("gradients_pathogen_PCA_joelike60cutoff_no0resistance.png", plot, width = 6.5, height = 3)

# Figure 3 - global avertable burden
avertable_by_pathogen <- read.csv("Outputs/10pc_avertable_burden_by_pathogen_joelike_weighted_lower_region_v2.csv", stringsAsFactors = FALSE)
optimistic_by_pathogen <- read.csv("Outputs/10pc_avertable_burden_by_pathogen_joelike_weighted_upper_region_optimistic_overall.csv", stringsAsFactors = FALSE)
pessimistic_by_pathogen <- read.csv("Outputs/10pc_avertable_burden_by_pathogen_joelike_weighted_upper_region_pessimistic_overall.csv", stringsAsFactors = FALSE)

# Helper function to format the burden columns consistently
format_burden <- function(df) {
  df$burden_fmt <- paste0(
    formatC(round(df$avertable_burden, -floor(log10(abs(df$avertable_burden))) + 1), format = "f", big.mark = ",", digits = 0), 
    " (", 
    formatC(round(df$lower_bound, -floor(log10(abs(df$lower_bound))) + 1), format = "f", big.mark = ",", digits = 0), 
    " to ", 
    formatC(round(df$upper_bound, -floor(log10(abs(df$upper_bound))) + 1), format = "f", big.mark = ",", digits = 0), 
    ")"
  )
  return(df)
}

# Apply formatting
avertable_by_pathogen <- format_burden(avertable_by_pathogen)
optimistic_by_pathogen <- format_burden(optimistic_by_pathogen)
pessimistic_by_pathogen <- format_burden(pessimistic_by_pathogen)

# Add scenario labels
avertable_by_pathogen$Scenario <- "Main"
optimistic_by_pathogen$Scenario <- "Optimistic"
pessimistic_by_pathogen$Scenario <- "Pessimistic"

# Combine datasets into one
all_data <- bind_rows(avertable_by_pathogen, optimistic_by_pathogen, pessimistic_by_pathogen)

# Order the pathogens strictly based on the main avertable burden
pathogen_order <- avertable_by_pathogen$pathogen[order(avertable_by_pathogen$avertable_burden)]
# Reverse the order here because facet_wrap places the first factor level at the top
all_data$pathogen <- factor(all_data$pathogen, levels = rev(pathogen_order))

# Order the scenarios so "Main" appears at the top of the group, then optimistic, then pessimistic
all_data$Scenario <- factor(all_data$Scenario, levels = c("Pessimistic", "Optimistic", "Main"))

# Calculate uniform text X-coordinate based on the overall maximums, and ADD it to the dataframe
max_x <- max(all_data$upper_bound, na.rm = TRUE)
all_data$text_pos_x <- max_x + 0.07 * (max_x - min(all_data$lower_bound, na.rm = TRUE))


# Plot
plot <- ggplot(all_data, aes(x = avertable_burden, y = Scenario, fill = Scenario)) +
  # Map y directly to Scenario to eliminate the need for position_dodge
  geom_bar(aes(color=Scenario), stat = "identity", width = 0.8) +
  geom_errorbar(aes(xmin = lower_bound, xmax = upper_bound, color = Scenario), 
                width = 0.3) +
  
  # Map fontface to Scenario for selective bolding
  geom_text(aes(x = text_pos_x, label = burden_fmt, color = Scenario, size = Scenario, 
                fontface = Scenario),
            hjust = 0) +
  
  # FILLS: Main gets dark grey, sub-bars get white (making them look hollow)
  scale_fill_manual(values = c("Main" = "grey50", 
                               "Optimistic" = "white", 
                               "Pessimistic" = "white")) +
  
  # COLORS: Outlines and text use these. Main is black, others are coloured
  scale_color_manual(values = c("Main" = "black", 
                                "Optimistic" = "gray30",  
                                "Pessimistic" = "gray30")) + 
  
  # TEXT SIZES (Set exactly to 10 points to match the strip.text size)
  scale_size_manual(values = c("Main" = 8 / .pt, 
                               "Optimistic" = 8 / .pt, 
                               "Pessimistic" = 8 / .pt)) +
                               
  # FONT WEIGHT: Main is bold, sub-bars are standard text
  scale_discrete_manual(aesthetics = "fontface", 
                        values = c("Main" = "bold", 
                                   "Optimistic" = "plain", 
                                   "Pessimistic" = "plain")) +

  # Tightly control y-axis padding to bring labels closer to the top bar (Main)
  scale_y_discrete(expand = expansion(add = c(0.1, 0.1))) +

  scale_x_continuous(
    labels = function(x) format(x, big.mark = ",", scientific = FALSE),
    breaks = seq(0, 53000, by = 10000) # Creates evenly spaced ticks
  ) +

  # Facet the plot to move pathogen names on top of the bars
  facet_wrap(~ pathogen, ncol = 1, scales = "free_y") +

  labs(x = "Deaths averted", y = "") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 10, family = "Helvetica"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10, family = "Helvetica"),
    
    # Hide the y-axis text and ticks
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    # --- NEW: Add the x-axis line and tick marks ---
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length.x = unit(4, "points"),
    # -----------------------------------------------
    
    # Style the new top-stacked labels and remove bottom margin to bring text flush with the panel
    strip.text = element_text(size = 8, family = "Helvetica", face = "italic", hjust = 0, margin = margin(t = 0, b = 1)),
    strip.background = element_blank(),
    panel.spacing = unit(0.3, "lines"), # Keep grouped facets visually tight
    text = element_text(family = "Helvetica"),
    legend.position = "none",
    plot.margin = margin(2.5, 100, 5.5, 5.5, "points")
  ) +
  coord_cartesian(xlim = c(-100, 53000), clip = "off")

ggsave("Figure3.pdf", plot, width = 6.5, height = 9.3, units = "in")
# -----------------------------------------------------------------------------
# PowerPoint Slides Version of Figure 3 (Dynamic Ranges & Sorting)
# -----------------------------------------------------------------------------
cat("Generating PowerPoint slides versions of Figure 3...\n")

# Define the scenarios to iterate over
scenarios <- c("Main", "Optimistic", "Pessimistic")

for (scen in scenarios) {
  cat(paste("Generating slide for scenario:", scen, "...\n"))
  
  # Filter data for the specific scenario
  slide_data <- all_data[all_data$Scenario == scen, ]
  
  # 1. ORDER PATHOGENS: Greatest at the top, least at the bottom.
  # For ggplot standard y-axis, the last factor level is placed at the top.
  # Sorting ascending ensures the lowest value is bottom (first level) and highest is top (last level).
  slide_data$pathogen <- factor(slide_data$pathogen, 
                                levels = slide_data$pathogen[order(slide_data$avertable_burden)])
  
  # 2. DYNAMIC TEXT POSITIONING
  # Find the maximum upper bound for THIS scenario to align text in a neat column
  max_bound <- max(slide_data$upper_bound, na.rm = TRUE)
  
  # Align text uniformly just to the right of the longest error bar for this scenario
  slide_data$text_pos_x <- max_bound + 0.05 * max_bound
  
  # Determine fill and color based on scenario to match original aesthetic
  bar_fill <- "grey50"
  bar_color <- "black"
  font_weight <- "plain"
  
  plot_slide <- ggplot(slide_data, aes(x = avertable_burden, y = pathogen)) +
    geom_bar(stat = "identity", width = 0.7, fill = bar_fill, color = bar_color) +
    geom_errorbar(aes(xmin = lower_bound, xmax = upper_bound), 
                  width = 0.25, color = bar_color) +
    
    # Add text labels positioned uniformly to the right
    geom_text(aes(x = text_pos_x, label = burden_fmt),
              hjust = 0, size = 12 / .pt, color = bar_color, 
              fontface = font_weight) +
    
    # 3. DYNAMIC X-AXIS RANGES
    scale_x_continuous(
      labels = function(x) format(x, big.mark = ",", scientific = FALSE),
      breaks = function(x) pretty(x, n = 5) # Let R compute 5 neat tick marks automatically
    ) +
    
    labs(x = "Deaths averted", y = "", title = paste(scen, "Scenario")) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      
      axis.title.x = element_text(size = 14, family = "Helvetica"),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 12, family = "Helvetica"),
      
      # Show pathogen names and italicize them
      axis.text.y = element_text(size = 14, family = "Helvetica", face = "italic"),
      axis.ticks.y = element_blank(),
      
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.ticks.x = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length.x = unit(4, "points"),
      
      text = element_text(family = "Helvetica"),
      legend.position = "none",
      
      # Wide right margin so the text column fits outside the plot panel
      plot.margin = margin(10, 160, 10, 10, "points"),
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold", margin = margin(b = 15))
    ) +
    # Set the x-axis limits dynamically based on this scenario's max bound. 
    # clip = "off" ensures the text isn't cut off when it draws into the right margin.
    coord_cartesian(xlim = c(0, max_bound * 1.02), clip = "off")
  
  # Save slide with PowerPoint dimensions
  filename <- paste0("Figure3_", scen, "_Slide_narrow.pdf")
  ggsave(filename, plot_slide, width = 8, height = 7.5, units = "in")
  
  cat(paste("Slide saved as", filename, "\n"))
}

cat("All dynamic PowerPoint slide versions for Figure 3 generated successfully.\n")

# # Supplementary Figure 1 - drug-class effects in HICs and LMICs
# drug_class_gradients_HIC <- setNames(as.vector(read.csv("Outputs/database_gradients_ATC3_PCA_joelike60cutoff_no0resistance_HIC.csv", stringsAsFactors = FALSE)[,2]), 
#                             read.csv("Outputs/database_gradients_ATC3_PCA_joelike60cutoff_no0resistance_HIC.csv", stringsAsFactors = FALSE)[,1])
# drug_class_lowerCI_HIC <- setNames(as.vector(read.csv("Outputs/database_lowerCI_ATC3_PCA_joelike60cutoff_no0resistance_HIC.csv", stringsAsFactors = FALSE)[,2]),
#                             read.csv("Outputs/database_lowerCI_ATC3_PCA_joelike60cutoff_no0resistance_HIC.csv", stringsAsFactors = FALSE)[,1])
# drug_class_upperCI_HIC <- setNames(as.vector(read.csv("Outputs/database_upperCI_ATC3_PCA_joelike60cutoff_no0resistance_HIC.csv", stringsAsFactors = FALSE)[,2]),
#                             read.csv("Outputs/database_upperCI_ATC3_PCA_joelike60cutoff_no0resistance_HIC.csv", stringsAsFactors = FALSE)[,1])         
# drug_class_gradients_LMIC <- setNames(as.vector(read.csv("Outputs/database_gradients_ATC3_PCA_joelike60cutoff_no0resistance_LMIC.csv", stringsAsFactors = FALSE)[,2]), 
#                             read.csv("Outputs/database_gradients_ATC3_PCA_joelike60cutoff_no0resistance_LMIC.csv", stringsAsFactors = FALSE)[,1])
# drug_class_lowerCI_LMIC <- setNames(as.vector(read.csv("Outputs/database_lowerCI_ATC3_PCA_joelike60cutoff_no0resistance_LMIC.csv", stringsAsFactors = FALSE)[,2]),
#                             read.csv("Outputs/database_lowerCI_ATC3_PCA_joelike60cutoff_no0resistance_LMIC.csv", stringsAsFactors = FALSE)[,1])
# drug_class_upperCI_LMIC <- setNames(as.vector(read.csv("Outputs/database_upperCI_ATC3_PCA_joelike60cutoff_no0resistance_LMIC.csv", stringsAsFactors = FALSE)[,2]),
#                             read.csv("Outputs/database_upperCI_ATC3_PCA_joelike60cutoff_no0resistance_LMIC.csv", stringsAsFactors = FALSE)[,1])

# segment_data_HIC <- data.frame(
#     Antibiotic = antibiotics,
#     Response = drug_class_gradients_HIC[antibiotics],
#     Lower_CI = drug_class_lowerCI_HIC[antibiotics],
#     Upper_CI = drug_class_upperCI_HIC[antibiotics],
#     Income = "HICs"
# )
# segment_data_LMIC <- data.frame(
#     Antibiotic = antibiotics,
#     Response = drug_class_gradients_LMIC[antibiotics],
#     Lower_CI = drug_class_lowerCI_LMIC[antibiotics],
#     Upper_CI = drug_class_upperCI_LMIC[antibiotics],
#     Income = "LMICs"
# )
# segment_data_Income <- rbind(segment_data_HIC, segment_data_LMIC)

# # replace Antibiotic names with antibiotic_names
# segment_data_Income$Antibiotic <- as.character(segment_data_Income$Antibiotic)  # Convert factors to characters
# for (i in 1:nrow(segment_data_Income)) {
#     segment_data_Income$Antibiotic[i] <- atc_names[[segment_data_Income$Antibiotic[i]]]
# }
# segment_data_Income$Antibiotic <- factor(segment_data_Income$Antibiotic, levels = antibiotic_names)

# # forest plot
# plot <- ggplot(segment_data_Income, aes(x = Response, y = Antibiotic)) +
#     geom_point(stat = "identity") +
#     geom_errorbar(aes(xmin = Lower_CI, xmax = Upper_CI), width = 0.2, position = position_dodge(0.9)) +
#     labs(x = "Elasticity of resistance with respect to consumption", y = "") +  # Added title
#     geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#     # add title
#     # ggtitle("Elasticity of Resistance with Respect to Consumption by Antibiotic and Income Group") +
#     theme_minimal() +
#     theme(panel.background = element_rect(fill = "white")) +
#     theme(axis.title.x = element_text(size = 10)) +
#     theme(axis.title.y = element_text(size = 10)) +
#     theme(axis.text.x = element_text(size = 10)) +
#     theme(axis.text.y = element_markdown(size = 8)) +
#     # theme(plot.title = element_text(size = 20)) +  # Set title size
#     theme(strip.text = element_text(size = 10)) +  # Set facet titles size
#     facet_wrap(Income ~ .)
# ggsave("gradients_drugs_income_PCA_joelike60cutoff_no0resistance_response_fmt.png", plot, width = 6.5, height = 3)

# # Supplementary Figure 2 - meta-analysis hospital drug effects
# results_hospital <- read.csv("Outputs/Nagorsen_gradients_pathogen_ATC3_PCA_joelike_hospital_to_all.csv", stringsAsFactors = FALSE)
# class_gradients_hospital <- setNames(as.vector(read.csv("Outputs/Nagorsen_gradients_ATC3_PCA_joelike_hospital_to_all.csv", stringsAsFactors = FALSE)[,2]), 
#                             read.csv("Outputs/Nagorsen_gradients_ATC3_PCA_joelike_hospital_to_all.csv", stringsAsFactors = FALSE)[,1])
# class_lowerCI_hospital <- setNames(as.vector(read.csv("Outputs/Nagorsen_lowerCI_ATC3_PCA_joelike_hospital_to_all.csv", stringsAsFactors = FALSE)[,2]),
#                             read.csv("Outputs/Nagorsen_lowerCI_ATC3_PCA_joelike_hospital_to_all.csv", stringsAsFactors = FALSE)[,1])
# class_upperCI_hospital <- setNames(as.vector(read.csv("Outputs/Nagorsen_upperCI_ATC3_PCA_joelike_hospital_to_all.csv", stringsAsFactors = FALSE)[,2]),
#                             read.csv("Outputs/Nagorsen_upperCI_ATC3_PCA_joelike_hospital_to_all.csv", stringsAsFactors = FALSE)[,1])
# antibiotics_hospital <- unique(results_hospital$Antibiotic)

# segment_data_hospital <- data.frame(
#     Antibiotic = antibiotics_hospital,
#     Response = class_gradients_hospital[antibiotics_hospital],
#     Lower_CI = class_lowerCI_hospital[antibiotics_hospital],
#     Upper_CI = class_upperCI_hospital[antibiotics_hospital]
# )

# results_hospital$Segment <- 0
# results_hospital <- results_hospital %>%
#     select(Antibiotic, Pathogen, Response, Lower_CI, Upper_CI, Segment)
# segment_data_hospital$Pathogen <- " Overall"
# segment_data_hospital$Segment <- 1
# segment_data_hospital <- segment_data_hospital %>%
#     select(Antibiotic, Pathogen, Response, Lower_CI, Upper_CI, Segment)
# # combine results and segment_data
# hospital_results <- rbind(results_hospital, segment_data_hospital)

# hospital_results <- hospital_results %>%
#     filter(!is.na(Lower_CI) & !is.na(Upper_CI))

# pathogen_classes <- c(" Overall",
# "Acinetobacter spp.","E. coli","E. faecalis","E. faecium","Enterococcus spp.","Enterobacter spp.",
# "H. influenzae","K. pneumoniae","Morganella spp.","N. gonorrhoeae","P. aeruginosa",
# "S. agalactiae","S. aureus","S. pneumoniae","S. pyogenes","Salmonella spp.")
# pathogen_names <- c(" Overall",
# "*Acinetobacter* spp","*E. coli*","*E. faecalis*","*E. faecium*","*Enterococcus* spp","*Enterobacter* spp",
# "*H. influenzae*","*K. pneumoniae*","*Morganella* spp","*N. gonorrhoeae*","*P. aeruginosa*",
# "*S. agalactiae*","*S. aureus*","*S. pneumoniae*","*S. pyogenes*","*Salmonella* spp")
# pathogen_names <- setNames(pathogen_names, pathogen_classes)
# hospital_results$Pathogen <- as.character(hospital_results$Pathogen)  # Convert factors to characters
# for (i in 1:nrow(hospital_results)) {
#     hospital_results$Pathogen[i] <- pathogen_names[[hospital_results$Pathogen[i]]]
# }
# # replace Antibiotic names with antibiotic_names
# hospital_results$Antibiotic <- as.character(hospital_results$Antibiotic)  # Convert factors to characters
# for (i in 1:nrow(hospital_results)) {
#     hospital_results$Antibiotic[i] <- atc_names[[hospital_results$Antibiotic[i]]]
# }
# hospital_results$Antibiotic <- factor(hospital_results$Antibiotic, levels = rev(antibiotic_names))
# # add column with formatted response with CI
# hospital_results$Response_fmt <- paste0(round(hospital_results$Response, 2), " (", round(hospital_results$Lower_CI, 2), "–", round(hospital_results$Upper_CI, 2), ")")
# print(nchar(hospital_results$Response_fmt))
# for (i in 1:nrow(hospital_results)) {
#     hospital_results$Response_fmt[i] <- paste0(hospital_results$Pathogen[i], "...", paste(rep(".", floor(1.49*(19 - nchar(hospital_results$Response_fmt[i])))),collapse=""), hospital_results$Response_fmt[i])
# }
# hospital_results$Response_fmt <- factor(hospital_results$Response_fmt, levels = rev(unique(hospital_results$Response_fmt)))

# # forest plot
# plot <- ggplot(hospital_results, aes(x = Response, y = Response_fmt, color = as.factor(Segment))) +
#     geom_point(stat = "identity") +
#     geom_errorbar(aes(xmin = Lower_CI, xmax = Upper_CI), width = 0.2, position = position_dodge(0.9)) +
#     labs(x = "Elasticity of resistance with respect to consumption", y = "") +  # Added title
#     scale_color_manual(values = c("0" = "grey", "1" = "black"), guide = FALSE) +  # Disable the color legend
#     geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#     # add title
#     # ggtitle("Elasticity of Resistance with Respect to Consumption by Antibiotic and Pathogen") +
#     theme_minimal() +
#     theme(panel.background = element_rect(fill = "white")) +
#     theme(axis.title.x = element_text(size = 10)) +
#     theme(axis.title.y = element_text(size = 10)) +
#     theme(axis.text.x = element_text(size = 10)) +
#     theme(axis.text.y = element_markdown(size = 8)) +
#     # theme(plot.title = element_text(size = 20)) +  # Set title size
#     theme(strip.text = element_text(size = 10)) +  # Set facet titles size
#     ggforce::facet_col(Antibiotic ~ ., scales = "free_y", space = "free", drop = TRUE)
# ggsave("gradients_drugs_pathogen_PCA_joelike_hospital_to_all_response_fmt.png", plot, width = 6.5, height = 5)

# # Supplementary Figure XX - avertable burden by pathogen and region
# avertable_by_pathogen_region <- read.csv(
#         "Outputs/10pc_avertable_burden_by_pathogen_and_region_joelike_weighted_upper_region_v2.csv", 
#         stringsAsFactors = FALSE)
# print(avertable_by_pathogen_region)
# pathogen_classes <- c(
#         "Acinetobacter spp.", "Citrobacter spp.", "E. coli", "E. faecalis", 
#         "E. faecium", "Enterococcus spp.", "Enterobacter spp.",
#         "H. influenzae", "K. pneumoniae", "Morganella spp.", 
#         "N. gonorrhoeae", "P. aeruginosa", "Proteus spp.",
#         "S. agalactiae", "S. aureus", "S. pneumoniae", "S. pyogenes", 
#         "Salmonella spp.", "Serratia spp.", "Shigella spp.")
# pathogen_expressions <- list(
#         expression(italic("Acinetobacter") ~ "spp."), 
#         expression(italic("Citrobacter") ~ "spp."), 
#         expression(italic("E. coli")),
#         expression(italic("E. faecalis")), 
#         expression(italic("E. faecium")), 
#         expression(italic("Enterococcus") ~ "spp."),
#         expression(italic("Enterobacter") ~ "spp."), 
#         expression(italic("H. influenzae")),
#         expression(italic("K. pneumoniae")), 
#         expression(italic("Morganella") ~ "spp."),
#         expression(italic("N. gonorrhoeae")), 
#         expression(italic("P. aeruginosa")), 
#         expression(italic("Proteus") ~ "spp."),
#         expression(italic("S. agalactiae")), 
#         expression(italic("S. aureus")),
#         expression(italic("S. pneumoniae")), 
#         expression(italic("S. pyogenes")),
#         expression(italic("Salmonella") ~ "spp."), 
#         expression(italic("Serratia") ~ "spp."), 
#         expression(italic("Shigella") ~ "spp.")
# )
# print(unique(avertable_by_pathogen_region$pathogen))

# # sort pathogens by total avertable_burden taken from avertable_by_pathogen above
# avertable_by_pathogen <- avertable_by_pathogen %>%
#         arrange(desc(avertable_burden))
# pathogen_classes_sorted <- avertable_by_pathogen$pathogen
# avertable_by_pathogen_region$pathogen <- factor(
#         avertable_by_pathogen_region$pathogen, levels = pathogen_classes_sorted)
# print(unique(avertable_by_pathogen_region$pathogen))
# # get sum of avertable burden in each region and sort
# avertable_by_region <- avertable_by_pathogen_region %>%
#         group_by(region) %>%
#         summarise(avertable_burden_per_100k = sum(avertable_burden_per_100k), 
#                                                 lower_bound = sum(lower_bound_per_100k), 
#                                                 upper_bound = sum(upper_bound_per_100k))
# avertable_by_region <- avertable_by_region %>%
#         arrange(desc(avertable_burden_per_100k))
# avertable_by_pathogen_region$region <- factor(
#         avertable_by_pathogen_region$region, levels = avertable_by_region$region)

# # print sum of avertable burden in each region
# print(avertable_by_region)

# pathogen_expressions <- setNames(pathogen_expressions, pathogen_classes)
# avertable_by_pathogen_region$pathogen <- as.character(
#         avertable_by_pathogen_region$pathogen)
# for (i in seq_len(nrow(avertable_by_pathogen_region))) {
#         pathogen_key <- avertable_by_pathogen_region$pathogen[i]
#         avertable_by_pathogen_region$pathogen[i] <- 
#                 as.character(pathogen_expressions[[pathogen_key]])
# }
# # make sure the levels are preserved as expressions
# avertable_by_pathogen_region$pathogen <- factor(
#         avertable_by_pathogen_region$pathogen, 
#         levels = rev(sapply(pathogen_classes_sorted, 
#                                                                                         function(x) as.character(pathogen_expressions[[x]]))))



# plot <- ggplot(avertable_by_pathogen_region, 
#                                                          aes(x = avertable_burden_per_100k, y = pathogen)) +
#         geom_bar(stat = "identity", fill = "grey50") +
#         geom_errorbar(aes(xmin = lower_bound_per_100k, xmax = upper_bound_per_100k), 
#                                                                 width = 0.2, position = position_dodge(0.9)) +
#         labs(x = "Deaths averted per 100,000 persons", y = "") +
#         theme_minimal() +
#         theme(panel.background = element_rect(fill = "white")) +
#         theme(panel.grid = element_blank()) +
#         theme(axis.title.x = element_text(size = 10, family = "Helvetica")) +
#         theme(axis.title.y = element_text(size = 10, family = "Helvetica")) +
#         theme(axis.text.x = element_text(size = 10, family = "Helvetica")) +
#         theme(axis.text.y = element_text(size = 8, family = "Helvetica")) +
#         theme(strip.text = element_text(size = 8, family = "Helvetica")) +
#         theme(text = element_text(family = "Helvetica")) +
#         facet_wrap(~ region, ncol = 2) +
#         scale_y_discrete(labels = function(x) parse(text = x))
# ggsave("10pc_avertable_attributable_burden_per_100k_by_pathogen_region_joelike_upper_region_v2.png", 
#                          plot, width = 6.5, height = 9, units = "in")

# # Supplementary Figure XX - avertable burden by drug and region
# avertable_by_drug_region <- read.csv(
#         "Outputs/10pc_avertable_burden_by_drug_and_region_joelike_weighted_upper_region_v2.csv", 
#         stringsAsFactors = FALSE)
# avertable_by_drug_region$drug <- as.character(avertable_by_drug_region$drug)
# for (i in seq_len(nrow(avertable_by_drug_region))) {
#         if (avertable_by_drug_region$drug[i] %in% classes) {
#                 avertable_by_drug_region$drug[i] <- 
#                         atc_names[[avertable_by_drug_region$drug[i]]]
#         } else {
#                 avertable_by_drug_region$drug[i] <- paste0("Other")
#         }
# }
# # remove rows with drug == "Other"
# avertable_by_drug_region <- avertable_by_drug_region %>%
#         filter(drug != "Other")
# avertable_by_drug_region$drug <- factor(avertable_by_drug_region$drug, 
#                                                                                                                                                                 levels = antibiotic_names)

# # get sum of avertable burden in each region and sort
# avertable_by_region <- avertable_by_drug_region %>%
#         group_by(region) %>%
#         summarise(avertable_burden_per_100k = sum(avertable_burden_per_100k), 
#                                                 lower_bound = sum(lower_bound_per_100k), 
#                                                 upper_bound = sum(upper_bound_per_100k))
# avertable_by_region <- avertable_by_region %>%
#         arrange(desc(avertable_burden_per_100k))
# avertable_by_drug_region$region <- factor(avertable_by_drug_region$region, 
#                                                                                                                                                                         levels = avertable_by_region$region)
# print(avertable_by_drug_region)
# # print sum of avertable burden for each drug
# avertable_by_drug <- avertable_by_drug_region %>%
#         group_by(drug) %>%
#         summarise(avertable_burden_per_100k = sum(avertable_burden_per_100k), 
#                                                 lower_bound = sum(lower_bound_per_100k), 
#                                                 upper_bound = sum(upper_bound_per_100k))
# avertable_by_drug <- avertable_by_drug %>%
#         arrange(desc(avertable_burden_per_100k))
# print(avertable_by_drug)

# # add formatted burden with CI column
# avertable_by_drug_region$burden_fmt <- paste0(
#         format(round(avertable_by_drug_region$avertable_burden_per_100k, 2), 
#                                  nsmall = 2), 
#         " (", 
#         format(round(avertable_by_drug_region$lower_bound_per_100k, 2), 
#                                  nsmall = 2), 
#         "–", 
#         format(round(avertable_by_drug_region$upper_bound_per_100k, 2), 
#                                  nsmall = 2), 
#         ")")

# plot <- ggplot(avertable_by_drug_region, 
#                                                          aes(x = avertable_burden_per_100k, y = drug)) +
#         geom_bar(stat = "identity", fill = "grey50") +
#         geom_errorbar(aes(xmin = lower_bound_per_100k, xmax = upper_bound_per_100k), 
#                                                                 width = 0.2, position = position_dodge(0.9)) +
#         geom_text(aes(label = burden_fmt),
#                                                 hjust = 0, size = 2.5, 
#                                                 x = max(avertable_by_drug_region$upper_bound_per_100k) + 
#                                                         0.07 * (max(avertable_by_drug_region$upper_bound_per_100k) - 
#                                                                                         min(avertable_by_drug_region$lower_bound_per_100k))) +
#         labs(x = "Deaths averted per 100,000 persons", y = "") +
#         theme_minimal() +
#         theme(panel.background = element_rect(fill = "white")) +
#         theme(panel.grid = element_blank()) +
#         theme(axis.title.x = element_text(size = 10, family = "Helvetica")) +
#         theme(axis.title.y = element_text(size = 10, family = "Helvetica")) +
#         theme(axis.text.x = element_text(size = 10, family = "Helvetica")) +
#         theme(axis.text.y = element_markdown(size = 8, family = "Helvetica")) +
#         theme(strip.text = element_text(size = 8, family = "Helvetica")) +
#         theme(text = element_text(family = "Helvetica")) +
#         coord_cartesian(clip = "off") +
#         theme(plot.margin = margin(5.5, 80, 5.5, 5.5, "points")) +
#         ggforce::facet_col(region ~ ., scales = "free_y", space = "free", 
#                                                                                  drop = TRUE)
# ggsave("10pc_avertable_attributable_burden_by_drug_region_joelike_upper_region_v2.png", 
#                          plot, width = 6.5, height = 7, units = "in")
