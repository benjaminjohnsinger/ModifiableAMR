#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(grid)
})

source("utils.R")

classes <- c("J01A", "J01B", "J01C", "J01D", "J01E", "J01F", "J01G", "J01M")
antibiotic_names <- c(
  "Tetracyclines", "Glycopeptides and Lipopeptides", "Penicillins",
  "Non-Penicillin Beta-Lactams", "Sulfonamides and Trimethoprim",
  "Macrolides", "Aminoglycosides", "Quinolones"
)
atc_names_map <- setNames(antibiotic_names, classes)

out_dir <- "Outputs_forestplot"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pathogen_path <- "Outputs/Nagorsen_gradients_pathogen_ATC3_PCA_canonical_hospital_to_all_filtered.csv"
summary_path <- "Outputs/Nagorsen_gradients_ATC3_PCA_canonical_hospital_to_all_filtered.csv"
lower_path <- "Outputs/Nagorsen_lowerCI_ATC3_PCA_canonical_hospital_to_all_filtered.csv"
upper_path <- "Outputs/Nagorsen_upperCI_ATC3_PCA_canonical_hospital_to_all_filtered.csv"
bootstrap_path <- "Outputs/Nagorsen_gradients_bootstraps_ATC3_PCA_canonical_hospital_to_all_filtered.csv"

pathogen_data <- fread(pathogen_path) %>%
  filter(!is.na(Lower_CI) & !is.na(Upper_CI)) %>%
  mutate(Response_fmt = sprintf("%.2f (%.2f, %.2f)", Response, Lower_CI, Upper_CI))

summary_grad <- setNames(as.vector(read.csv(summary_path)[, 2]), read.csv(summary_path)[, 1])
summary_lower <- setNames(as.vector(read.csv(lower_path)[, 2]), read.csv(lower_path)[, 1])
summary_upper <- setNames(as.vector(read.csv(upper_path)[, 2]), read.csv(upper_path)[, 1])

antibiotics <- unique(pathogen_data$Antibiotic)
summary_data <- data.frame(
  Antibiotic = antibiotics,
  Response = summary_grad[antibiotics],
  Lower_CI = summary_lower[antibiotics],
  Upper_CI = summary_upper[antibiotics],
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(Lower_CI) & !is.na(Upper_CI)) %>%
  mutate(Response_fmt = sprintf("%.2f (%.2f, %.2f)", Response, Lower_CI, Upper_CI))

bootstrap_data <- data.frame()
if (file.exists(bootstrap_path)) {
  bootstrap_data <- fread(bootstrap_path)
  if (!"Gradient" %in% names(bootstrap_data)) {
    stop("Expected column 'Gradient' in bootstrap file: ", bootstrap_path)
  }
  bootstrap_data <- bootstrap_data %>%
    mutate(Antibiotic = atc_names_map[Antibiotic]) %>%
    filter(Antibiotic %in% atc_names_map[summary_data$Antibiotic])
}

plot_rows <- list()
current_order <- 1

for (abx_code in classes) {
  if (!abx_code %in% summary_data$Antibiotic) next

  abx_name <- atc_names_map[[abx_code]]
  path_subset <- pathogen_data %>% filter(Antibiotic == abx_code) %>% arrange(Pathogen)
  n_pathogens <- nrow(path_subset)

  plot_rows[[length(plot_rows) + 1]] <- data.frame(
    Antibiotic = abx_name,
    Plot_Label = abx_name,
    Order = current_order,
    Type = "Header",
    stringsAsFactors = FALSE
  )
  current_order <- current_order + 1

  if (n_pathogens > 0) {
    path_subset$Plot_Label <- paste0("  ", path_subset$Pathogen)
    path_subset$Order <- current_order:(current_order + n_pathogens - 1)
    path_subset$Type <- "Pathogen"
    path_subset$Antibiotic <- atc_names_map[path_subset$Antibiotic]
    plot_rows[[length(plot_rows) + 1]] <- path_subset
    current_order <- current_order + n_pathogens
  }

  if (n_pathogens > 1) {
    total_subset <- summary_data %>%
      filter(Antibiotic == abx_code) %>%
      mutate(Antibiotic = atc_names_map[Antibiotic])
    if (nrow(total_subset) > 0) {
      total_subset$Plot_Label <- "  Total"
      total_subset$Order <- current_order
      total_subset$Type <- "Total"
      plot_rows[[length(plot_rows) + 1]] <- total_subset
      current_order <- current_order + 1
    }
  }

  plot_rows[[length(plot_rows) + 1]] <- data.frame(
    Antibiotic = abx_name,
    Plot_Label = "",
    Order = current_order,
    Type = "Spacer",
    stringsAsFactors = FALSE
  )
  current_order <- current_order + 1
}

plot_df <- rbindlist(plot_rows, fill = TRUE) %>% arrange(Order)

shading_df <- plot_df %>%
  filter(Order %% 2 == 1) %>%
  mutate(ymin_rect = Order - 0.5, ymax_rect = Order + 0.5)

plot_violins <- data.frame()
if (nrow(bootstrap_data) > 0) {
  bootstrap_data$Plot_Label <- "  Total"
  plot_violins <- left_join(
    bootstrap_data,
    plot_df %>% filter(Type == "Total") %>% select(Antibiotic, Plot_Label, Order),
    by = c("Antibiotic", "Plot_Label")
  )
}

p <- ggplot() +
  geom_rect(
    data = shading_df,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin_rect, ymax = ymax_rect),
    fill = "grey95"
  ) +
  {
    if (nrow(plot_violins) > 0) {
      geom_violin(
        data = plot_violins,
        aes(x = Gradient, y = Order, group = Order),
        fill = "grey",
        alpha = 0.8,
        trim = TRUE,
        width = 3
      )
    }
  } +
  geom_pointrange(
    data = plot_df %>% filter(Type == "Pathogen"),
    aes(x = Response, y = Order, xmin = Lower_CI, xmax = Upper_CI),
    shape = 15,
    size = 0.5
  ) +
  geom_text(
    data = plot_df %>% filter(Type == "Header"),
    aes(x = -0.8, y = Order, label = Plot_Label),
    hjust = 0,
    size = 8 / .pt,
    fontface = "bold"
  ) +
  geom_text(
    data = plot_df %>% filter(Type == "Pathogen"),
    aes(x = -0.8, y = Order, label = Plot_Label),
    hjust = 0,
    size = 8 / .pt,
    fontface = "italic"
  ) +
  geom_text(
    data = plot_df %>% filter(Type == "Total"),
    aes(x = -0.8, y = Order, label = Plot_Label),
    hjust = 0,
    size = 8 / .pt,
    fontface = "bold"
  ) +
  geom_text(
    data = plot_df %>% filter(Type == "Pathogen"),
    aes(x = 1.3, y = Order, label = Response_fmt),
    hjust = 1,
    size = 8 / .pt
  ) +
  geom_text(
    data = plot_df %>% filter(Type == "Total"),
    aes(x = 1.3, y = Order, label = Response_fmt),
    hjust = 1,
    size = 8 / .pt,
    fontface = "bold"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("segment", x = -0.5, xend = 1, y = Inf, yend = Inf, color = "black", linewidth = 0.5) +
  scale_x_continuous("Elasticity", limits = c(-0.8, 1.3), breaks = c(-0.5, 0, 0.5, 1)) +
  scale_y_continuous(name = NULL, breaks = NULL, trans = "reverse", expand = expansion(mult = 0.01, add = 0.5)) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length.x = unit(4, "points"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing.y = unit(0.1, "lines")
  ) +
  coord_cartesian(clip = "off")

main_out <- file.path(out_dir, "supplementary_figure2b_violin_forest_hospital_to_all_filtered.pdf")

pdf(file = main_out, width = 6.5, height = 9.3)
print(p)
dev.off()

message("Wrote: ", main_out)
