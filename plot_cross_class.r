# Load required libraries
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(forcats)

# --- 1. Define Names and Mapping ---
# Using the mapping from your original plotting.r
classes <- c("J01A", "J01B", "J01C", "J01D", "J01E", "J01F", "J01G", "J01M")
antibiotic_names <- c(
  "Tetracyclines", "Glycopeptides and Lipopeptides", "Penicillins",
  "Non-Penicillin Beta-Lactams", "Sulfonamides and Trimethoprim",
  "Macrolides", "Aminoglycosides", "Quinolones"
)
atc_names_map <- setNames(antibiotic_names, classes)

# --- 2. Load and Prepare Data ---
cat("Loading cross-class effects data...\n")
data <- fread("Outputs/database_cross_class_effects_glmnet_weighted_glmnet_test.csv")

# Format the variables and strings
data <- data %>%
  mutate(
    Target_Name = atc_names_map[Target_Antibiotic],
    Cross_Name = atc_names_map[Cross_Class_Antibiotic],
    Response_fmt = sprintf("%.2f (%.2f, %.2f)", Coefficient, Lower_CI, Upper_CI)
  )

pathogens <- unique(data$Pathogen)

# --- 3. Generate Forest Plot per Pathogen ---
for (pathogen_name in pathogens) {
  
  # Filter data for the current pathogen
  path_data <- data %>% filter(Pathogen == pathogen_name)
  
  # Determine bounds for x-axis limits dynamically to space the text
  min_x <- min(path_data$Lower_CI, na.rm = TRUE)
  max_x <- max(path_data$Upper_CI, na.rm = TRUE)
  
  x_range <- max_x - min_x
  text_left_x <- min_x - (x_range * 0.45)
  text_right_x <- max_x + (x_range * 0.45)
  
  # Build plot layout data frame
  plot_data_list <- list()
  current_order <- 1
  
  target_abxs <- unique(path_data$Target_Name)
  target_abxs <- target_abxs[!is.na(target_abxs)]
  
  for (i in seq_along(target_abxs)) {
    target <- target_abxs[i]
    
    # 1. Header row for the focal (Target) antibiotic
    plot_data_list[[length(plot_data_list) + 1]] <- data.frame(
      Target_Name = target,
      Plot_Label = target,
      Order = current_order,
      Type = "Header",
      group_index = i,
      Coefficient = NA,
      Lower_CI = NA,
      Upper_CI = NA,
      Response_fmt = ""
    )
    current_order <- current_order + 1
    
    # 2. Rows for each cross-class antibiotic
    sub_data <- path_data %>% 
      filter(Target_Name == target) %>%
      arrange(Cross_Name)
    
    if (nrow(sub_data) > 0) {
      sub_data <- sub_data %>%
        mutate(
          Plot_Label = paste0("  ", Cross_Name), # Indent the text
          Order = current_order:(current_order + n() - 1),
          Type = "CrossClass",
          group_index = i
        )
      plot_data_list[[length(plot_data_list) + 1]] <- sub_data
      current_order <- current_order + nrow(sub_data)
    }
    
    # 3. Blank spacer row between blocks
    plot_data_list[[length(plot_data_list) + 1]] <- data.frame(
      Target_Name = target,
      Plot_Label = "",
      Order = current_order,
      Type = "Spacer",
      group_index = i,
      Coefficient = NA,
      Lower_CI = NA,
      Upper_CI = NA,
      Response_fmt = ""
    )
    current_order <- current_order + 1
  }
  
  # --- 1. Define Fixed Display Bounds ---
  view_min <- -2
  view_max <- 2
  
  # Calculate text positions with a fixed offset from the display bounds
  # Adjust the '2.0' padding if your pathogen names or estimates are very long
  text_left_x <- view_min - 2.0 
  text_right_x <- view_max + 2.0

  # Combine layout components and arrange
  plot_df <- rbindlist(plot_data_list, fill = TRUE) %>%
    arrange(Order)
    
  # Create a row-shading dataframe (odd rows only)
  shading_df <- plot_df %>%
    filter(Order %% 2 == 1) %>%
    mutate(
      ymin_rect = Order - 0.5,
      ymax_rect = Order + 0.5
    )
  
# --- 2. Update Plotting Layers ---
  p <- ggplot() +
    # [Keep geom_rect and geom_vline as they are]
    
    # Restrict the error bars: values outside -2 to 2 will be clipped at the plot edge
    geom_pointrange(
      data = plot_df %>% filter(Type == "CrossClass"),
      aes(x = Coefficient, y = Order, xmin = Lower_CI, xmax = Upper_CI),
      shape = 15, size = 0.5
    ) +
    
    # Left Text: Headers and Cross-class Labels
    geom_text(
      data = plot_df %>% filter(Type == "Header"),
      aes(x = text_left_x, y = Order, label = Plot_Label),
      hjust = 0, size = 8 / .pt, fontface = "bold"
    ) +
    geom_text(
      data = plot_df %>% filter(Type == "CrossClass"),
      aes(x = text_left_x, y = Order, label = Plot_Label),
      hjust = 0, size = 8 / .pt
    ) +
    
    # Right Text: Numeric Format Label
    geom_text(
      data = plot_df %>% filter(Type == "CrossClass"),
      aes(x = text_right_x, y = Order, label = Response_fmt),
      hjust = 1, size = 8 / .pt
    ) +
    
    # --- 3. Update the X-Axis Scale ---
    # We set the scale limits to include the text, but set 'breaks' only for the -2 to 2 range
    scale_x_continuous(
      "Coefficient",
      limits = c(text_left_x, text_right_x),
      breaks = seq(-2, 2, by = 0.5)
    ) +
    
    # Use coord_cartesian to prevent long error bars from expanding the plot beyond our text
    # clip = "on" will truncate any error bars that extend past text_left/right
    coord_cartesian(xlim = c(text_left_x, text_right_x), clip = "on") +
    
    # Theming / Layout Configuration
    scale_x_continuous(
      "Coefficient",
      limits = c(text_left_x, text_right_x)
    ) +
    scale_y_continuous(
      name = element_blank(),
      breaks = NULL,
      trans = "reverse", # Top-down ordering
      expand = expansion(mult = 0.01, add = 0.5)
    ) +
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
      panel.spacing.y = unit(0.1, "lines"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12)
    ) +
    coord_cartesian(clip = "off") +
    labs(title = paste("Cross-Class Resistance Effects -", pathogen_name))
  
  # Format a valid filename
  safe_name <- gsub("[ /]", "_", pathogen_name)
  safe_name <- gsub("\\.", "", safe_name)
  
  # Dynamically calculate plot height based on rows
  plot_height <- max(4, current_order * 0.22)
  
  # Save the plot
  png(filename = paste0("Cross_Class_Forest_", safe_name, ".png"), width = 8, height = plot_height, units = "in", res = 300)
  print(p)
  dev.off()
  
  cat("Generated and saved plot for:", pathogen_name, "\n")
}