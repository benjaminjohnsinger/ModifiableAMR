# -----------------------------------------------------------------------------
# Figure 1 (Advanced ggplot2 Version - FINAL Corrected & Refined)
# - Uses user's file names and data loading logic
# - Uses scale_y_continuous() to fix overplotting
# - Violins for summary
# - Inline drug headers
# - Right-aligned text for estimates
# - NEW: Reversed antibiotic order
# - NEW: Blank line after 'Total'
# - NEW: Smaller violin plots
# - NEW: 'Response_fmt' annotations moved further right
# - NEW: Bold text for "Total" rows
# - NEW: Text placed via geom_text to allow full-width shading and identical sizes
# -----------------------------------------------------------------------------
cat("Generating Figure 1 (Advanced ggplot2 version - FINAL Corrected & Refined)...\n")

# Load necessary libraries (add ggplot2 if not already present)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(forcats) # For factor manipulation

# --- 1. Define Names and Order ---
# Use the same mappings and order from your original script
classes <- c("J01A", "J01B", "J01C", "J01D", "J01E", "J01F", "J01G", "J01M")
antibiotic_names <- c(
  "Tetracyclines", "Glycopeptides and Lipopeptides", "Penicillins",
  "Non-Penicillin Beta-Lactams", "Sulfonamides and Trimethoprim",
  "Macrolides", "Aminoglycosides", "Quinolones"
)
atc_names_map <- setNames(antibiotic_names, classes)

antibiotic_list<- c("Quinolones", "Aminoglycosides", "Non-Penicillin Beta-Lactams", 
                                     "Penicillins", "Macrolides", "Sulfonamides and Trimethoprim", 
                                     "Tetracyclines")

# --- 2. Load and Prepare All Data Sources ---
cat("Loading data for Figure 1...\n")

# Load Pathogen-specific data
pathogen_data <- fread("Outputs/database_gradients_pathogen_ATC3_PCA_joelike_weighted_lagged.csv") %>%
  mutate(Antibiotic = atc_names_map[Antibiotic]) %>%
  filter(Antibiotic %in% antibiotic_list) %>%
  mutate(Response_fmt = sprintf("%.2f (%.2f, %.2f)", Response, Lower_CI, Upper_CI))

# Load Drug-summary data (for the text labels)
summary_data <- fread("Outputs/database_gradients_ATC3_PCA_joelike_weighted_all_lagged.csv") %>%
  mutate(Antibiotic = atc_names_map[as.character(V1)]) %>%
  filter(Antibiotic %in% antibiotic_list)
summary_data_LowerCI <- fread("Outputs/database_lowerCI_ATC3_PCA_joelike_weighted_all_lagged.csv") %>%
  mutate(Antibiotic = atc_names_map[as.character(V1)]) %>%
  filter(Antibiotic %in% antibiotic_list)
summary_data_UpperCI <- fread("Outputs/database_upperCI_ATC3_PCA_joelike_weighted_all_lagged.csv") %>%
    mutate(Antibiotic = atc_names_map[as.character(V1)]) %>%
    filter(Antibiotic %in% antibiotic_list)
summary_data <- summary_data %>%
    left_join(summary_data_LowerCI %>% select(Antibiotic, x) %>% rename(Lower_CI = x), by = "Antibiotic") %>%
    left_join(summary_data_UpperCI %>% select(Antibiotic, x) %>% rename(Upper_CI = x), by = "Antibiotic") %>%
    mutate(Response_fmt = sprintf("%.2f (%.2f, %.2f)", x, Lower_CI, Upper_CI)) %>%
    rename(Response = x)

# Load Bootstrap gradient data (for the violins)
tryCatch({
  bootstrap_data <- fread("Outputs/database_gradients_bootstraps_ATC3_PCA_joelike_weighted_all_lagged.csv") %>%
    mutate(Antibiotic = atc_names_map[Antibiotic]) %>%
    filter(Antibiotic %in% antibiotic_list)
  cat("Bootstrap data loaded for Figure 1.\n")
}, error = function(e) {
  stop("Could not load 'Outputs/database_gradients_bootstraps_ATC3_PCA_joelike_weighted_lagged.csv'. 
       Please generate this file from your modeling script first.")
})

# --- 3. Build the Custom Y-Axis Data Frame ---
cat("Building plot layout for Figure 1...\n")
plot_data_list <- list()
current_order <- 1

# --- REFINEMENT 1: Iterate through the antibiotics in the desired order (top to bottom) ---
for (i in seq_along(antibiotic_list)) { 
  abx <- antibiotic_list[i]

  # 1. Add Antibiotic Header Row
  plot_data_list[[length(plot_data_list) + 1]] <- data.frame(
    Antibiotic = abx,
    Plot_Label = abx,
    Order = current_order,
    Type = "Header",
    group_index = i
  )
  current_order <- current_order + 1
  
  # 2. Add Pathogen Rows
  path_subset <- pathogen_data %>%
    filter(Antibiotic == abx) %>%
    arrange(Pathogen)
  
  if (nrow(path_subset) > 0) {
    path_subset$Plot_Label <- paste0("  ", path_subset$Pathogen) # Indent
    path_subset$Order <- current_order:(current_order + nrow(path_subset) - 1)
    path_subset$Type = "Pathogen"
    path_subset$group_index = i
    
    plot_data_list[[length(plot_data_list) + 1]] <- path_subset
    current_order <- current_order + nrow(path_subset)
  }
  
  # 3. Add "Total" Row
  summary_subset <- summary_data %>% filter(Antibiotic == abx)
  
  if (nrow(summary_subset) > 0) {
    summary_subset$Plot_Label <- "  Total" # Indent
    summary_subset$Order <- current_order
    summary_subset$Type = "Total"
    summary_subset$group_index = i
    
    plot_data_list[[length(plot_data_list) + 1]] <- summary_subset
    current_order <- current_order + 1
  }

  # --- REFINEMENT 2: Add a blank line after "Total" ---
  plot_data_list[[length(plot_data_list) + 1]] <- data.frame(
    Antibiotic = abx, # Associate with the current antibiotic for completeness
    Plot_Label = "", # Blank label for this row
    Order = current_order,
    Type = "Spacer", # A new type to identify this row
    group_index = i
  )
  current_order <- current_order + 1
}

# Combine all pieces and sort by the unique Order
plot_df <- rbindlist(plot_data_list, fill = TRUE) %>%
  arrange(Order)

cat("Creating row-shading layout for Figure 1...\n")
shading_df <- plot_df %>%
  filter(Order %% 2 == 1) %>% # Keep only odd-numbered rows
  mutate(
    ymin_rect = Order - 0.5, # Boundary is 0.5 above the Order number
    ymax_rect = Order + 0.5  # Boundary is 0.5 below the Order number
  )
# --- 4. Prepare Violin Data ---
# Join bootstrap data with plot_df to get the unique 'Order' for each violin
cat("Preparing violin data for Figure 1...\n")
bootstrap_data$Plot_Label <- "  Total" 
plot_violins <- left_join(
  bootstrap_data, 
  plot_df %>% filter(Type == "Total") %>% select(Antibiotic, Plot_Label, Order), 
  by = c("Antibiotic", "Plot_Label")
)

# --- 6. Generate the Plot (with Corrected Y-Axis) ---
cat("Generating plot for Figure 1...\n")
p_final_plot <- ggplot() +

  geom_rect(
    data = shading_df,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin_rect, ymax = ymax_rect),
    fill = "grey95"
  ) +
  
  # Add the violin plots for "Total"
  geom_violin(
    data = plot_violins,
    aes(x = Consumption, y = Order, group = Order),
    fill = "grey", alpha = 0.8, trim = TRUE, 
    width = 3
  ) +
  
  # Add the point estimates for individual pathogens
  geom_pointrange(
    data = plot_df %>% filter(Type == "Pathogen"),
    aes(x = Response, y = Order, xmin = Lower_CI, xmax = Upper_CI),
    shape = 15,
    size = 0.5
  ) +
  
  # --- LEFT TEXT (Replaces axis.text.y) ---
  geom_text(
    data = plot_df %>% filter(Type == "Header"),
    aes(x = -2.3, y = Order, label = Plot_Label),
    hjust = 0, size = 8 / .pt, fontface = "bold"
  ) +
  geom_text(
    data = plot_df %>% filter(Type == "Pathogen"),
    aes(x = -2.3, y = Order, label = Plot_Label),
    hjust = 0, size = 8 / .pt, fontface = "italic"
  ) +
  geom_text(
    data = plot_df %>% filter(Type == "Total"),
    aes(x = -2.3, y = Order, label = Plot_Label),
    hjust = 0, size = 8 / .pt, fontface = "bold"
  ) +
  
  # --- RIGHT TEXT ---
  geom_text(
    data = plot_df %>% filter(Type == "Pathogen"),
    aes(x = 4.5, y = Order, label = Response_fmt),
    hjust = 1, size = 8 / .pt, fontface = "plain"
  ) +
  geom_text(
    data = plot_df %>% filter(Type == "Total"),
    aes(x = 4.5, y = Order, label = Response_fmt),
    hjust = 1, size = 8 / .pt, fontface = "bold"
  ) +
  
# Add a vertical line at zero
geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +

# --- NEW: Manually draw the x-axis line exactly from -1 to 3 ---
annotate("segment", x = -1, xend = 3, y = Inf, yend = Inf, color = "black", linewidth = 0.5) +

# --- Theming and Axis ---

# Expand limits drastically to the left to act as a placeholder for the text
scale_x_continuous(
  "Elasticity",
  limits = c(-2.3, 4.5), 
  breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3)
) +

# Remove the explicit labels mapped to the Y axis since geom_text handles them
scale_y_continuous(
  name = element_blank(),
  breaks = NULL,          
  trans = "reverse",         
  expand = expansion(mult = 0.01, add = 0.5)
) +

# Apply a clean theme
theme_minimal() +

# Apply custom theme elements
theme(
  axis.text.y = element_blank(), 
  axis.ticks.y = element_blank(),
  
  # --- NEW: Turn off default full line, but keep x-axis ticks ---
  axis.line.x = element_blank(), 
  axis.ticks.x = element_line(color = "black", linewidth = 0.5),
  axis.ticks.length.x = unit(4, "points"),
  # --------------------------------------------------------------
  
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.spacing.y = unit(0.1, "lines") 
) +
coord_cartesian(clip = "off")

# --- 7. Save the Plot ---
pdf(file = "Figure1_lagged.pdf", width = 6.5, height = 9.3)
print(p_final_plot)
dev.off()


# -----------------------------------------------------------------------------
# PowerPoint Slides Version of Figure 1
# -----------------------------------------------------------------------------
cat("Generating PowerPoint slides version of Figure 1...\n")

# Define antibiotic groupings for slides
slide_groups <- list(
  slide1 = c("Quinolones", "Aminoglycosides"),
  slide2 = c("Non-Penicillin Beta-Lactams", "Penicillins"),
  slide3 = c("Macrolides", "Sulfonamides and Trimethoprim", "Tetracyclines")
)

# Function to create filtered plot data for specific antibiotics
create_slide_data <- function(selected_antibiotics) {
  filtered_df <- plot_df %>%
    filter(Antibiotic %in% selected_antibiotics | Type %in% c("Spacer"))
  
  filtered_df$Order <- seq_len(nrow(filtered_df))
  
  shading_df_slide <- filtered_df %>%
    filter(Order %% 2 == 1) %>%
    mutate(
      ymin_rect = Order - 0.5,
      ymax_rect = Order + 0.5
    )
  
  violin_data_slide <- plot_violins %>%
    filter(Antibiotic %in% selected_antibiotics) %>%
    left_join(
      filtered_df %>% 
        filter(Type == "Total") %>% 
        select(Antibiotic, Plot_Label, Order), 
      by = c("Antibiotic", "Plot_Label"),
      suffix = c("", "_new")
    ) %>%
    mutate(Order = Order_new) %>%
    select(-Order_new)
  
  return(list(
    plot_data = filtered_df,
    shading_data = shading_df_slide,
    violin_data = violin_data_slide
  ))
}

# Function to create plot for a slide
create_slide_plot <- function(slide_data, title) {
  
  ggplot() +
    geom_rect(
      data = slide_data$shading_data,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin_rect, ymax = ymax_rect),
      fill = "grey95"
    ) +
    
    geom_violin(
      data = slide_data$violin_data,
      aes(x = Consumption, y = Order, group = Order),
      fill = "grey", alpha = 0.8, trim = TRUE, 
      width = 3
    ) +
    
    geom_pointrange(
      data = slide_data$plot_data %>% filter(Type == "Pathogen"),
      aes(x = Response, y = Order, xmin = Lower_CI, xmax = Upper_CI),
      shape = 15,
      size = 0.7
    ) +
    
    # --- LEFT TEXT --- (Set specifically to size 12 for slides)
    geom_text(data = slide_data$plot_data %>% filter(Type == "Header"),
              aes(x = -2.6, y = Order, label = Plot_Label),
              hjust = 0, size = 12 / .pt, fontface = "bold") +
    geom_text(data = slide_data$plot_data %>% filter(Type == "Pathogen"),
              aes(x = -2.6, y = Order, label = Plot_Label),
              hjust = 0, size = 12 / .pt, fontface = "italic") +
    geom_text(data = slide_data$plot_data %>% filter(Type == "Total"),
              aes(x = -2.6, y = Order, label = Plot_Label),
              hjust = 0, size = 12 / .pt, fontface = "bold") +
    
    # --- RIGHT TEXT ---
    geom_text(data = slide_data$plot_data %>% filter(Type == "Pathogen"),
              aes(x = 4.5, y = Order, label = Response_fmt),
              hjust = 1, size = 12 / .pt, fontface = "plain") +
    geom_text(data = slide_data$plot_data %>% filter(Type == "Total"),
              aes(x = 4.5, y = Order, label = Response_fmt),
              hjust = 1, size = 12 / .pt, fontface = "bold") +
    
geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    
    # --- NEW: Manually draw the x-axis line exactly from -1 to 3 ---
    annotate("segment", x = -1, xend = 3, y = Inf, yend = Inf, color = "black", linewidth = 0.5) +
    
    scale_x_continuous(
      "Elasticity",
      limits = c(-2.6, 4.5),
      breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3)
    ) +
    
    scale_y_continuous(
      name = element_blank(),
      breaks = NULL,
      trans = "reverse",
      expand = expansion(mult = 0.01, add = 0.5)
    ) +
    
    theme_minimal() +
    
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      
      # --- NEW: Turn off default full line, but keep x-axis ticks ---
      axis.line.x = element_blank(), 
      axis.ticks.x = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length.x = unit(4, "points"),
      # --------------------------------------------------------------
      
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.spacing.y = unit(0.1, "lines"),
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold")
    ) +
    coord_cartesian(clip = "off") +
    ggtitle(title)
}

# Generate slides
for (i in 1:3) {
  cat(paste("Generating slide", i, "...\n"))
  
  slide_antibiotics <- slide_groups[[i]]
  slide_data <- create_slide_data(slide_antibiotics)
  
  slide_title <- paste("Estimated Elasticities (Part", i, ")")
  slide_plot <- create_slide_plot(slide_data, slide_title)
  
  # Save slide with PowerPoint dimensions
  filename <- paste0("Figure1_Slide", i, "_narrow_spillover_lagged.pdf")
  pdf(file = filename, width = 8, height = 7.5)
  print(slide_plot)
  dev.off()
  
  cat(paste("Slide", i, "saved as", filename, "\n"))
}

cat("All PowerPoint slides generated successfully.\n")
cat("Figure 1 (Advanced ggplot2 version - FINAL Corrected & Refined) saved to 'Figure1_lagged.pdf'.\n\n")
# -----------------------------------------------------------------------------
# Supplementary Figure 2b: Meta-analysis of Hospital Drug Effects (ggplot2)
# - Replicates the style of Figure 1 (violins for summary, custom y-axis)
# - Uses the *Nagorsen* data sources from 'Gemini_plotting.r'
# - REMOVED: forestplot elements (diamond, error bar) from "Total" row.
# - NEW: Don't plot "Total" when there's just one pathogen species
# - NEW: Matches Figure 1 formatting with full-width stripes and overlapping
#   drug class names
# -----------------------------------------------------------------------------
cat("Generating Supplementary Figure 2b (ggplot2 version)...\n")

# --- 1. Load and Prepare All Data Sources for Supp Fig 2b ---
cat("Loading data for Supp Fig 2b...\n")

# Load Pathogen-specific data (hospital)
pathogen_data_hosp <- fread(
  "Outputs/Nagorsen_gradients_pathogen_ATC3_PCA_joelike_hospital_to_all_filtered.csv"
) %>%
  filter(!is.na(Lower_CI) & !is.na(Upper_CI)) %>%
  mutate(
    Response_fmt = sprintf("%.2f (%.2f, %.2f)", Response, Lower_CI, Upper_CI)
  )

# Load Drug-summary data (hospital, for the text labels)
hosp_summary_grad <- setNames(
  as.vector(read.csv(
    "Outputs/Nagorsen_gradients_ATC3_PCA_joelike_hospital_to_all_filtered.csv"
  )[, 2]),
  read.csv(
    "Outputs/Nagorsen_gradients_ATC3_PCA_joelike_hospital_to_all_filtered.csv"
  )[, 1]
)
hosp_summary_lower <- setNames(
  as.vector(read.csv(
    "Outputs/Nagorsen_lowerCI_ATC3_PCA_joelike_hospital_to_all_filtered.csv"
  )[, 2]),
  read.csv(
    "Outputs/Nagorsen_lowerCI_ATC3_PCA_joelike_hospital_to_all_filtered.csv"
  )[, 1]
)
hosp_summary_upper <- setNames(
  as.vector(read.csv(
    "Outputs/Nagorsen_upperCI_ATC3_PCA_joelike_hospital_to_all_filtered.csv"
  )[, 2]),
  read.csv(
    "Outputs/Nagorsen_upperCI_ATC3_PCA_joelike_hospital_to_all_filtered.csv"
  )[, 1]
)

antibiotics_hosp <- unique(pathogen_data_hosp$Antibiotic)
summary_data_hosp <- data.frame(
  Antibiotic = antibiotics_hosp,
  Response = hosp_summary_grad[antibiotics_hosp],
  Lower_CI = hosp_summary_lower[antibiotics_hosp],
  Upper_CI = hosp_summary_upper[antibiotics_hosp]
) %>%
  filter(!is.na(Lower_CI) & !is.na(Upper_CI)) %>%
  mutate(
    Response_fmt = sprintf(
      "%.2f (%.2f, %.2f)",
      Response,
      Lower_CI,
      Upper_CI
    )
  )

# Load Bootstrap gradient data (for the violins)
tryCatch({
  bootstrap_data_hosp <- fread(
    "Outputs/Nagorsen_gradients_bootstraps_ATC3_PCA_joelike_hospital_to_all_filtered.csv"
  )
  bootstrap_data_hosp <- bootstrap_data_hosp %>%
    mutate(Antibiotic = atc_names_map[Antibiotic]) %>%
    filter(Antibiotic %in% atc_names_map[summary_data_hosp$Antibiotic])
  cat("Bootstrap data loaded for Supp Fig 2b.\n")
}, error = function(e) {
  cat(paste(
    "Warning: Could not load bootstrap data.",
    "Violin plots will be missing.\n"
  ))
  bootstrap_data_hosp <- data.frame()
})

# --- 2. Build the Custom Y-Axis Data Frame for Supp Fig 2b ---
cat("Building plot layout for Supp Fig 2b...\n")
plot_data_list_hosp <- list()
current_order_hosp <- 1

for (abx_code in classes) {
  if (abx_code %in% summary_data_hosp$Antibiotic) {
    abx_name <- atc_names_map[[abx_code]]

    # Count how many pathogens exist for this antibiotic
    path_subset_hosp <- pathogen_data_hosp %>%
      filter(Antibiotic == abx_code) %>%
      arrange(Pathogen)

    n_pathogens <- nrow(path_subset_hosp)

    # 1. Add Antibiotic Header Row
    plot_data_list_hosp[[length(plot_data_list_hosp) + 1]] <- data.frame(
      Antibiotic = abx_name,
      Plot_Label = abx_name,
      Order = current_order_hosp,
      Type = "Header"
    )
    current_order_hosp <- current_order_hosp + 1

    # 2. Add Pathogen Rows
    if (nrow(path_subset_hosp) > 0) {
      path_subset_hosp$Plot_Label <- paste0("  ", path_subset_hosp$Pathogen)
      path_subset_hosp$Order <- current_order_hosp:(
        current_order_hosp + nrow(path_subset_hosp) - 1
      )
      path_subset_hosp$Type <- "Pathogen"
      path_subset_hosp$Antibiotic <- atc_names_map[path_subset_hosp$Antibiotic]

      plot_data_list_hosp[[length(plot_data_list_hosp) + 1]] <- path_subset_hosp
      current_order_hosp <- current_order_hosp + nrow(path_subset_hosp)
    }

    # 3. Add "Total" Row - ONLY if there's more than one pathogen
    if (n_pathogens > 1) {
      summary_subset_hosp <- summary_data_hosp %>%
        filter(Antibiotic == abx_code) %>%
        mutate(Antibiotic = atc_names_map[Antibiotic])

      if (nrow(summary_subset_hosp) > 0) {
        summary_subset_hosp$Plot_Label <- "  Total"
        summary_subset_hosp$Order <- current_order_hosp
        summary_subset_hosp$Type <- "Total"

        plot_data_list_hosp[[length(plot_data_list_hosp) + 1]] <- (
          summary_subset_hosp
        )
        current_order_hosp <- current_order_hosp + 1
      }
    }

    # 4. Add a blank line after "Total" (or after pathogens if no Total)
    plot_data_list_hosp[[length(plot_data_list_hosp) + 1]] <- data.frame(
      Antibiotic = abx_name,
      Plot_Label = "",
      Order = current_order_hosp,
      Type = "Spacer"
    )
    current_order_hosp <- current_order_hosp + 1
  }
}

# Combine all pieces and sort by the unique Order
plot_df_hosp <- rbindlist(plot_data_list_hosp, fill = TRUE) %>%
  arrange(Order)

cat("Creating row-shading layout for Supp Fig 2b...\n")
shading_df_hosp <- plot_df_hosp %>%
  filter(Order %% 2 == 1) %>%
  mutate(
    ymin_rect = Order - 0.5,
    ymax_rect = Order + 0.5
  )

# --- 3. Prepare Violin Data for Supp Fig 2b ---
cat("Preparing violin data for Supp Fig 2b...\n")
plot_violins_hosp <- data.frame()
if (nrow(bootstrap_data_hosp) > 0) {
  bootstrap_data_hosp$Plot_Label <- "  Total"
  plot_violins_hosp <- left_join(
    bootstrap_data_hosp,
    plot_df_hosp %>%
      filter(Type == "Total") %>%
      select(Antibiotic, Plot_Label, Order),
    by = c("Antibiotic", "Plot_Label")
  )
}

# --- 4. Generate the Plot for Supp Fig 2b ---
cat("Generating plot for Supp Fig 2b...\n")

p_supp_fig_2b <- ggplot() +
  geom_rect(
    data = shading_df_hosp,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin_rect, ymax = ymax_rect),
    fill = "grey95"
  ) +

  # Add the violin plots for "Total"
  {
    if (nrow(plot_violins_hosp) > 0) {
      geom_violin(
        data = plot_violins_hosp,
        aes(x = Consumption, y = Order, group = Order),
        fill = "grey",
        alpha = 0.8,
        trim = TRUE,
        width = 3
      )
    }
  } +

  # Add the point estimates for individual pathogens
  geom_pointrange(
    data = plot_df_hosp %>% filter(Type == "Pathogen"),
    aes(x = Response, y = Order, xmin = Lower_CI, xmax = Upper_CI),
    shape = 15,
    size = 0.5
  ) +

  # --- LEFT TEXT (Drug class headers and pathogen names) ---
  geom_text(
    data = plot_df_hosp %>% filter(Type == "Header"),
    aes(x = -0.8, y = Order, label = Plot_Label),
    hjust = 0,
    size = 8 / .pt,
    fontface = "bold"
  ) +
  geom_text(
    data = plot_df_hosp %>% filter(Type == "Pathogen"),
    aes(x = -0.8, y = Order, label = Plot_Label),
    hjust = 0,
    size = 8 / .pt,
    fontface = "italic"
  ) +
  geom_text(
    data = plot_df_hosp %>% filter(Type == "Total"),
    aes(x = -0.8, y = Order, label = Plot_Label),
    hjust = 0,
    size = 8 / .pt,
    fontface = "bold"
  ) +

  # --- RIGHT TEXT ---
  geom_text(
    data = plot_df_hosp %>% filter(Type == "Pathogen"),
    aes(x = 1.3, y = Order, label = Response_fmt),
    hjust = 1,
    size = 8 / .pt,
    fontface = "plain"
  ) +
  geom_text(
    data = plot_df_hosp %>% filter(Type == "Total"),
    aes(x = 1.3, y = Order, label = Response_fmt),
    hjust = 1,
    size = 8 / .pt,
    fontface = "bold"
  ) +

  # Add a vertical line at zero
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +

  # Manually draw the x-axis line
  annotate(
    "segment",
    x = -0.5,
    xend = 1,
    y = Inf,
    yend = Inf,
    color = "black",
    linewidth = 0.5
  ) +

  # Set the x-axis label and limits
  scale_x_continuous(
    "Elasticity",
    limits = c(-0.8, 1.3),
    breaks = c(-0.5, 0, 0.5, 1)
  ) +

  # Use scale_y_continuous() for numeric Order with reverse direction
  scale_y_continuous(
    name = element_blank(),
    breaks = NULL,
    trans = "reverse",
    expand = expansion(mult = 0.01, add = 0.5)
  ) +

  # Apply a clean theme
  theme_minimal() +

  # Apply custom theme elements
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

# --- 5. Save the Plot ---
pdf(
  file = "Outputs_forestplot/supp_fig2b_violin_forest_plot.pdf",
  width = 6.5,
  height = 9.3
)
print(p_supp_fig_2b)
dev.off()

# --- 6. Generate PowerPoint Slide Version ---
cat("Generating slide version of Supp Fig 2b...\n")

p_supp_fig_2b_slide <- ggplot() +
  geom_rect(
    data = shading_df_hosp,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin_rect, ymax = ymax_rect),
    fill = "grey95"
  ) +

  # Add the violin plots for "Total"
  {
    if (nrow(plot_violins_hosp) > 0) {
      geom_violin(
        data = plot_violins_hosp,
        aes(x = Consumption, y = Order, group = Order),
        fill = "grey",
        alpha = 0.8,
        trim = TRUE,
        width = 3
      )
    }
  } +

  # Add the point estimates for individual pathogens
  geom_pointrange(
    data = plot_df_hosp %>% filter(Type == "Pathogen"),
    aes(x = Response, y = Order, xmin = Lower_CI, xmax = Upper_CI),
    shape = 15,
    size = 0.5
  ) +

  # --- LEFT TEXT (Drug class headers and pathogen names) ---
  geom_text(
    data = plot_df_hosp %>% filter(Type == "Header"),
    aes(x = -0.8, y = Order, label = Plot_Label),
    hjust = 0,
    size = 12 / .pt,
    fontface = "bold"
  ) +
  geom_text(
    data = plot_df_hosp %>% filter(Type == "Pathogen"),
    aes(x = -0.8, y = Order, label = Plot_Label),
    hjust = 0,
    size = 12 / .pt,
    fontface = "italic"
  ) +
  geom_text(
    data = plot_df_hosp %>% filter(Type == "Total"),
    aes(x = -0.8, y = Order, label = Plot_Label),
    hjust = 0,
    size = 12 / .pt,
    fontface = "bold"
  ) +

  # --- RIGHT TEXT ---
  geom_text(
    data = plot_df_hosp %>% filter(Type == "Pathogen"),
    aes(x = 1.5, y = Order, label = Response_fmt),
    hjust = 1,
    size = 12 / .pt,
    fontface = "plain"
  ) +
  geom_text(
    data = plot_df_hosp %>% filter(Type == "Total"),
    aes(x = 1.5, y = Order, label = Response_fmt),
    hjust = 1,
    size = 12 / .pt,
    fontface = "bold"
  ) +

  # Add a vertical line at zero
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +

  # Manually draw the x-axis line
  annotate(
    "segment",
    x = -0.5,
    xend = 1,
    y = Inf,
    yend = Inf,
    color = "black",
    linewidth = 0.5
  ) +

  # Set the x-axis label and limits
  scale_x_continuous(
    "Elasticity",
    limits = c(-0.8, 1.5),
    breaks = c(-0.5, 0, 0.5, 1)
  ) +

  # Use scale_y_continuous() for numeric Order with reverse direction
  scale_y_continuous(
    name = element_blank(),
    breaks = NULL,
    trans = "reverse",
    expand = expansion(mult = 0.01, add = 0.5)
  ) +

  # Apply a clean theme
  theme_minimal() +

  # Apply custom theme elements with larger text for slides
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
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

# Save slide version with PowerPoint dimensions
pdf(
  file = "Outputs_forestplot/supp_fig2b_slide_narrow.pdf",
  width = 8,
  height = 7.5
)
print(p_supp_fig_2b_slide)
dev.off()

cat(
  "Supp Fig 2b saved to 'supp_fig2b_violin_forest_plot.pdf'.\n"
)
cat(
  "Supp Fig 2b slide version saved to 'supp_fig2b_slide_narrow.pdf'.\n\n"
)
# -----------------------------------------------------------------------------
# End of Supplementary Figure 2b
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Figure 2: Meta-analysis of Pathogen Elasticities (VIOLIN PLOT Version)
# - Uses the 15-pathogen list from the user's provided image
# - CORRECTED: Loads *database* CSVs, matching Gemini_plotting.r
# -----------------------------------------------------------------------------
cat("Generating Figure 2 (Pathogen Violin Plot version, 15 pathogens)...\n")

# --- 1. Define Pathogen Names and Order (from user image) ---
fig2_pathogen_codes <- c(
  "Acinetobacter spp.", "Salmonella spp.", "S. pneumoniae", "S. aureus", 
  "E. coli", "P. aeruginosa", "Enterococcus spp.", "K. pneumoniae", 
  "E. faecalis", "E. faecium", "Morganella spp.", "N. gonorrhoeae"
)
fig2_pathogen_names <- c(
  "Acinetobacter spp.", "Salmonella spp.", "S. pneumoniae", "S. aureus", 
  "E. coli", "P. aeruginosa", "Enterococcus spp.", "K. pneumoniae", 
  "E. faecalis", "E. faecium", "Morganella spp.", "N. gonorrhoeae"
)
pathogen_map_fig2 <- setNames(fig2_pathogen_names, fig2_pathogen_codes)

# --- 2. Load and Prepare Data for Figure 2 ---
cat("Loading pathogen-level summary and bootstrap data for Figure 2...\n")

# Load Summary Data (for text labels)
tryCatch({
  # CORRECTED file paths to use 'database' prefix
  grad_path <- fread("Outputs/database_gradients_pathogen_PCA_joelike_weighted.csv", 
                     col.names = c("Pathogen_Code", "Response"))
  lower_path <- fread("Outputs/database_lowerCI_pathogen_PCA_joelike_weighted.csv", 
                      col.names = c("Pathogen_Code", "Lower_CI"))
  upper_path <- fread("Outputs/database_upperCI_pathogen_PCA_joelike_weighted.csv", 
                      col.names = c("Pathogen_Code", "Upper_CI"))
  
  # Merge into one data.table for text labels
  plot_data_fig2_summary <- grad_path %>%
    merge(lower_path, by = "Pathogen_Code") %>%
    merge(upper_path, by = "Pathogen_Code") %>%
    mutate(Pathogen_Display = pathogen_map_fig2[Pathogen_Code]) %>% 
    filter(Pathogen_Code %in% fig2_pathogen_codes) %>% 
    mutate(
      Pathogen_Display = factor(Pathogen_Display, 
                               levels = rev(fig2_pathogen_names)), 
      Response_fmt = sprintf("%.2f (%.2f, %.2f)", Response, Lower_CI, Upper_CI)
    )
  
}, error = function(e) {
  stop("Could not load pathogen-level summary data.")
})

# Load Bootstrap Data (for violins)
tryCatch({
  # This file path was already correct
  bootstrap_data_fig2 <- fread("Outputs/database_gradients_bootstraps_pathogen_PCA_joelike_weighted.csv") %>%
    mutate(Pathogen_Display = pathogen_map_fig2[Pathogen]) %>% 
    filter(Pathogen %in% fig2_pathogen_codes) %>% 
    mutate(Pathogen_Display = factor(Pathogen_Display, 
                                    levels = rev(fig2_pathogen_names)))
  
  cat("Bootstrap data loaded for Figure 2.\n")
}, error = function(e) {
  stop("Could not load pathogen bootstrap data.")
})

# --- 3. Create Shading Dataframe for Figure 2 ---
shading_data_fig2 <- data.frame(Pathogen_Display = fig2_pathogen_names) %>% 
  mutate(
    ymin = row_number() - 0.5,
    ymax = row_number() + 0.5,
    is_shaded = row_number() %% 2 != 0 # Shade odd rows
  ) %>%
  filter(is_shaded)

# --- 4. Generate the Plot for Figure 2 ---
cat("Generating plot for Figure 2...\n")
text_x_pos_fig2 <- 2.5 # Set right-hand position for text
plot_x_min <- -1.1     # Set left-hand limit for plot

p_fig2_violin <- ggplot() +
  
  # Add zebra striping first
  geom_rect(
  data = shading_data_fig2,
  aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf, y = NULL),
  fill = "grey95",
  inherit.aes = FALSE
  ) +

  # Add the violin plots
  geom_violin(
  data = bootstrap_data_fig2,
  aes(x = Consumption, y = Pathogen_Display, group = Pathogen_Display),
  fill = "grey", alpha = 0.8, trim = TRUE,
  width = 2 # Standard violin width
  ) +
  
  # Add rectangles behind the text
  geom_rect(
    data = plot_data_fig2_summary,
    aes(ymin = as.numeric(Pathogen_Display) - 0.4, 
      ymax = as.numeric(Pathogen_Display) + 0.4,
      xmin = 1.83, xmax = text_x_pos_fig2),
    fill = ifelse(as.numeric(plot_data_fig2_summary$Pathogen_Display) %% 2 != 0, 
                 "grey95", "white"), 
    alpha = 0.9
  ) +

  geom_point(
  data = plot_data_fig2_summary,
  aes(x = Response, y = Pathogen_Display),
  shape = 18,
  size = 8 / .pt
  ) +

  # Add the Response_fmt text on the right
  geom_text(
  data = plot_data_fig2_summary,
  aes(x = text_x_pos_fig2, y = Pathogen_Display, label = Response_fmt),
  hjust = 1, # Right-justify
  size = 8 / .pt
  ) +
  
  # Add a vertical line at zero
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +

  # x axis running from -1 to 2
  annotate("segment", x = -1.005, xend = 2.005, y = -Inf, yend = -Inf, color = "black", linewidth = 1) +

  
  # --- Theming and Axis ---
  scale_x_continuous(
  "Elasticity",
  limits = c(plot_x_min, text_x_pos_fig2), # Adjusted limits for violins + text
  breaks = c(-1,-0.5, 0, 0.5, 1, 1.5, 2)
  ) +
  
  scale_y_discrete(
  name = element_blank(),
  expand = expansion(mult = 0.01, add = 0.5)
  ) +

  # Apply a clean theme
  theme_minimal() +
  
  # Apply custom theme elements
  theme(
  axis.text.y = element_text(face = "italic", size = 8, hjust = 0), 
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  axis.line.x = element_blank(),
  axis.ticks.x = element_line(color = "black", linewidth = 0.5),
  axis.ticks.length.x = unit(4, "points")
  )

# --- 5. Save the Plot ---
# Regular version
pdf(file = "Figure2.pdf", width = 6.5, height = 4.7)
print(p_fig2_violin)
dev.off()

# --- 6. Generate PowerPoint slide version ---
cat("Generating PowerPoint slide version of Figure 2...\n")

p_fig2_slide <- ggplot() +
  
  # Add zebra striping first
  geom_rect(
  data = shading_data_fig2,
  aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf, y = NULL),
  fill = "grey95",
  inherit.aes = FALSE
  ) +

  # Add the violin plots
  geom_violin(
  data = bootstrap_data_fig2,
  aes(x = Consumption, y = Pathogen_Display, group = Pathogen_Display),
  fill = "grey", alpha = 0.8, trim = TRUE,
  width = 2
  ) +
  
  # Add rectangles behind the text
  geom_rect(
    data = plot_data_fig2_summary,
    aes(ymin = as.numeric(Pathogen_Display) - 0.4, 
      ymax = as.numeric(Pathogen_Display) + 0.4,
      xmin = 1.72, xmax = text_x_pos_fig2),
    fill = ifelse(as.numeric(plot_data_fig2_summary$Pathogen_Display) %% 2 != 0, 
                 "grey95", "white"), 
    alpha = 0.9
  ) +

  geom_point(
  data = plot_data_fig2_summary,
  aes(x = Response, y = Pathogen_Display),
  shape = 18,
  size = 3
  ) +

  # Add the Response_fmt text on the right
  geom_text(
  data = plot_data_fig2_summary,
  aes(x = text_x_pos_fig2, y = Pathogen_Display, label = Response_fmt),
  hjust = 1,
  size = 4
  ) +
  
  # Add a vertical line at zero
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  
  # --- Theming and Axis ---
  scale_x_continuous(
  "Elasticity",
  limits = c(plot_x_min, text_x_pos_fig2),
  breaks = c(-1,-0.5, 0, 0.5, 1, 1.5, 2)
  ) +
  
  scale_y_discrete(
  name = element_blank(),
  expand = expansion(mult = 0.01, add = 0.5)
  ) +

  # Apply a clean theme
  theme_minimal() +
  
  # Apply custom theme elements for slide
  theme(
  axis.text.y = element_text(face = "italic", size = 14, hjust = 0), 
  axis.text.x = element_text(size = 14),
  axis.title.x = element_text(size = 16),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.grid.minor.x = element_blank()
  )

# Save slide version with PowerPoint dimensions
pdf(file = "Figure2_Slide_narrow.pdf", width = 8, height = 7.5)
print(p_fig2_slide)
dev.off()

cat("Figure 2 (Pathogen Violin Plot version, 15 pathogens) saved to 'Figure2.pdf'.\n")
cat("Figure 2 slide version saved to 'Figure2_Slide_narrow.pdf'.\n\n")
# -----------------------------------------------------------------------------
# End of Figure 2 (Corrected Violin Version)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Supplementary Figure 1: Separate panels for HICs and LMICs
# - Uses the same style as Figure 2 (violins)
# - Loads *database* CSVs, matching Gemini_plotting.r
# - NOW: Creates two separate panels side by side
# -----------------------------------------------------------------------------
cat("Generating Supplementary Figure 1 (HIC vs LMIC Separate Panels)...\n")

# Load required library for panel plots
library(patchwork)

# Toggle for slide version text size
slide_version <- TRUE  # Set to FALSE for regular version

# --- 1. Load and Prepare Data for Supp Fig 1 ---
cat("Loading pathogen-level summary and bootstrap data for Supp Fig 1...\n")
drug_class_gradients_hic <- setNames(
  as.vector(read.csv("Outputs/database_gradients_ATC3_PCA_joelike_weighted_HIC.csv", 
                     stringsAsFactors = FALSE)[, 2]), 
  read.csv("Outputs/database_gradients_ATC3_PCA_joelike_weighted_HIC.csv", 
           stringsAsFactors = FALSE)[, 1]
)
drug_class_lower_ci_hic <- setNames(
  as.vector(read.csv("Outputs/database_lowerCI_ATC3_PCA_joelike_weighted_HIC.csv", 
                     stringsAsFactors = FALSE)[, 2]),
  read.csv("Outputs/database_lowerCI_ATC3_PCA_joelike_weighted_HIC.csv", 
           stringsAsFactors = FALSE)[, 1]
)
drug_class_upper_ci_hic <- setNames(
  as.vector(read.csv("Outputs/database_upperCI_ATC3_PCA_joelike_weighted_HIC.csv", 
                     stringsAsFactors = FALSE)[, 2]),
  read.csv("Outputs/database_upperCI_ATC3_PCA_joelike_weighted_HIC.csv", 
           stringsAsFactors = FALSE)[, 1]
)
drug_class_bootstraps_hic <- fread(
  "Outputs/database_gradients_bootstraps_ATC3_PCA_joelike_weighted_HIC.csv"
)
drug_class_gradients_lmic <- setNames(
  as.vector(read.csv("Outputs/database_gradients_ATC3_PCA_joelike_weighted_LMIC.csv", 
                     stringsAsFactors = FALSE)[, 2]), 
  read.csv("Outputs/database_gradients_ATC3_PCA_joelike_weighted_LMIC.csv", 
           stringsAsFactors = FALSE)[, 1]
)
drug_class_lower_ci_lmic <- setNames(
  as.vector(read.csv("Outputs/database_lowerCI_ATC3_PCA_joelike_weighted_LMIC.csv", 
                     stringsAsFactors = FALSE)[, 2]),
  read.csv("Outputs/database_lowerCI_ATC3_PCA_joelike_weighted_LMIC.csv", 
           stringsAsFactors = FALSE)[, 1]
)
drug_class_upper_ci_lmic <- setNames(
  as.vector(read.csv("Outputs/database_upperCI_ATC3_PCA_joelike_weighted_LMIC.csv", 
                     stringsAsFactors = FALSE)[, 2]),
  read.csv("Outputs/database_upperCI_ATC3_PCA_joelike_weighted_LMIC.csv", 
           stringsAsFactors = FALSE)[, 1]
)
drug_class_bootstraps_lmic <- fread(
  "Outputs/database_gradients_bootstraps_ATC3_PCA_joelike_weighted_LMIC.csv"
)

# Prepare summary data frames for HIC and LMIC
plot_data_fig_s1_summary_hic <- data.frame(
  Antibiotic = names(drug_class_gradients_hic),
  Response = drug_class_gradients_hic,
  Lower_CI = drug_class_lower_ci_hic[names(drug_class_gradients_hic)],
  Upper_CI = drug_class_upper_ci_hic[names(drug_class_gradients_hic)]
) %>%
  mutate(Antibiotic = atc_names_map[Antibiotic]) %>%
  filter(Antibiotic %in% antibiotic_names) %>%
  mutate(Response_fmt = sprintf("%.2f (%.2f to %.2f)", Response, Lower_CI, Upper_CI))

plot_data_fig_s1_summary_lmic <- data.frame(
  Antibiotic = names(drug_class_gradients_lmic),
  Response = drug_class_gradients_lmic,
  Lower_CI = drug_class_lower_ci_lmic[names(drug_class_gradients_lmic)],
  Upper_CI = drug_class_upper_ci_lmic[names(drug_class_gradients_lmic)]
) %>%
  mutate(Antibiotic = atc_names_map[Antibiotic]) %>%
  filter(Antibiotic %in% antibiotic_names) %>%
  mutate(Response_fmt = sprintf("%.2f (%.2f to %.2f)", Response, Lower_CI, Upper_CI))

# --- 2. Create shading data frame for Supp Fig 1 ---
shading_data_fig_s1 <- data.frame(Antibiotic = antibiotic_names) %>% 
  mutate(
    ymin = row_number() - 0.5,
    ymax = row_number() + 0.5,
    is_shaded = row_number() %% 2 != 0 # Shade odd rows
  ) %>%
  filter(is_shaded)

# --- 3. Generate HIC Panel ---
cat("Generating HIC panel for Supp Fig 1...\n")
text_x_pos_fig_s1 <- 2.5 # Set right-hand position for text
plot_x_min_fig_s1 <- -1.2 # Set left-hand limit for plot

p_hic <- ggplot() +
  
  # Add zebra striping first
  geom_rect(
    data = shading_data_fig_s1,
    aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf, y = NULL),
    fill = "grey95",
    inherit.aes = FALSE
  ) +

  # Add the violin plots for HIC
  geom_violin(
    data = drug_class_bootstraps_hic %>%
      mutate(Antibiotic = atc_names_map[Antibiotic]) %>%
      filter(Antibiotic %in% antibiotic_names),
    aes(x = Consumption, y = Antibiotic, group = Antibiotic),
    fill = "grey", alpha = 0.6, trim = TRUE,
    width = 2
  ) +
  
  # Add rectangles behind the text
  geom_rect(
    data = plot_data_fig_s1_summary_hic,
    aes(ymin = as.numeric(factor(Antibiotic, levels = rev(antibiotic_names))) - 0.4, 
        ymax = as.numeric(factor(Antibiotic, levels = rev(antibiotic_names))) + 0.4,
        xmin = 1.35, xmax = text_x_pos_fig_s1),
    fill = ifelse(as.numeric(factor(plot_data_fig_s1_summary_hic$Antibiotic, 
                                   levels = rev(antibiotic_names))) %% 2 != 0, 
                  "grey95", "white"), 
    alpha = 0.9
  ) +

  geom_point(
    data = plot_data_fig_s1_summary_hic,
    aes(x = Response, y = Antibiotic),
    shape = 18,
    size = if (slide_version) 3 else 2
  ) +
  
  
  # Add the Response_fmt text on the right for HIC
  geom_text(
    data = plot_data_fig_s1_summary_hic,
    aes(x = text_x_pos_fig_s1, y = Antibiotic, label = Response_fmt),
    hjust = 1, # Right-justify
    size = if (slide_version) 4 else 2.5,
    color = "black"
  ) +
  
  # Add a vertical line at zero
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  
  # --- Theming and Axis ---
  scale_x_continuous(
    "Elasticity",
    limits = c(plot_x_min_fig_s1, text_x_pos_fig_s1),
    breaks = c(-1, -0.5, 0, 0.5, 1, 1.5)
  ) +
  
  scale_y_discrete(
    name = element_blank(),
    expand = expansion(mult = 0.01, add = 0.5)
  ) +
  
  # Apply a clean theme
  theme_minimal() +
  
  # Apply custom theme elements
  theme(
    axis.text.y = element_text(size = if (slide_version) 14 else 9),
    axis.text.x = element_text(size = if (slide_version) 14 else 9),
    axis.title.x = element_text(size = if (slide_version) 16 else 11),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(hjust = 0.5, size = if (slide_version) 16 else 11)
  ) +
  
  ggtitle("High-Income Countries")

# --- 4. Generate LMIC Panel ---
cat("Generating LMIC panel for Supp Fig 1...\n")

p_lmic <- ggplot() +
  
  # Add zebra striping first
  geom_rect(
    data = shading_data_fig_s1,
    aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf, y = NULL),
    fill = "grey95",
    inherit.aes = FALSE
  ) +

  # Add the violin plots for LMIC
  geom_violin(
    data = drug_class_bootstraps_lmic %>%
      mutate(Antibiotic = atc_names_map[Antibiotic]) %>%
      filter(Antibiotic %in% antibiotic_names),
    aes(x = Consumption, y = Antibiotic, group = Antibiotic),
    fill = "grey", alpha = 0.6, trim = TRUE,
    width = 2
  ) +

  # Add rectangles behind the text
  geom_rect(
    data = plot_data_fig_s1_summary_hic,
    aes(ymin = as.numeric(factor(Antibiotic, levels = rev(antibiotic_names))) - 0.4, 
        ymax = as.numeric(factor(Antibiotic, levels = rev(antibiotic_names))) + 0.4,
        xmin = 1.35, xmax = text_x_pos_fig_s1),
    fill = ifelse(as.numeric(factor(plot_data_fig_s1_summary_hic$Antibiotic, 
                                   levels = rev(antibiotic_names))) %% 2 != 0, 
                  "grey95", "white"), 
    alpha = 0.9
  ) +

  geom_point(
    data = plot_data_fig_s1_summary_lmic,
    aes(x = Response, y = Antibiotic),
    shape = 18,
    size = if (slide_version) 3 else 2
  ) +
  
  # Add the Response_fmt text on the right for LMIC
  geom_text(
    data = plot_data_fig_s1_summary_lmic,
    aes(x = text_x_pos_fig_s1, y = Antibiotic, label = Response_fmt),
    hjust = 1, # Right-justify
    size = if (slide_version) 4 else 2.5,
    color = "black"
  ) +
  
  # Add a vertical line at zero
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  
  # --- Theming and Axis ---
  scale_x_continuous(
    "Elasticity",
    limits = c(plot_x_min_fig_s1, text_x_pos_fig_s1),
    breaks = c(-1, -0.5, 0, 0.5, 1, 1.5)
  ) +
  
  scale_y_discrete(
    name = element_blank(),
    expand = expansion(mult = 0.01, add = 0.5)
  ) +
  
  # Apply a clean theme
  theme_minimal() +
  
  # Apply custom theme elements
  theme(
    axis.text.y = element_text(size = if (slide_version) 14 else 9),
    axis.text.x = element_text(size = if (slide_version) 14 else 9),
    axis.title.x = element_text(size = if (slide_version) 16 else 11),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(hjust = 0.5, size = if (slide_version) 16 else 11)
  ) +
  
  ggtitle("Low- and Middle-Income Countries")

# --- 5. Combine panels ---
cat("Combining panels for Supp Fig 1...\n")
p_supp_fig1_combined <- p_hic + p_lmic + 
  plot_layout(ncol = 1)

# --- 6. Save the Plot ---
# 8 by 7.5 inch pdf version
pdf(file = "Supplementary_Figure1_Slide_narrow.pdf", width = 8, height = 7.5)
print(p_supp_fig1_combined)
dev.off()

cat("Supplementary Figure 1 (HIC vs LMIC Separate Panels) saved to 'Supplementary_Figure1_Slide_narrow.pdf'.\n\n")
# -----------------------------------------------------------------------------
# End of Supplementary Figure 1
# -----------------------------------------------------------------------------
