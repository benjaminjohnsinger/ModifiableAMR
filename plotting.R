slide_version = FALSE
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

# Load config for output directory constants (sourced from run_generate_figures.R if available)
if (!exists("AMR_CONFIG")) source("config.R")
# Slide outputs go to a separate directory from manuscript PDFs
SLIDES_DIR <- AMR_CONFIG$output_dirs$slides
dir.create(SLIDES_DIR, recursive = TRUE, showWarnings = FALSE)

# Load necessary libraries (add ggplot2 if not already present)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(forcats) # For factor manipulation

# Global toggle: set AMR_EXCLUDE_N_GONORRHOEAE=1 when running make figures.
exclude_n_gonorrhoeae <- tolower(trimws(
  Sys.getenv("AMR_EXCLUDE_N_GONORRHOEAE", "0")
)) %in% c("1", "true", "t", "yes", "y")

filter_excluded_pathogen_rows <- function(df, data_name = "data") {
  if (!exclude_n_gonorrhoeae || !("Pathogen" %in% names(df))) {
    return(df)
  }

  before_n <- nrow(df)
  pathogen_values <- tolower(trimws(as.character(df[["Pathogen"]])))
  df <- df[!(pathogen_values %in% c(
    "n. gonorrhoeae",
    "neisseria gonorrhoeae"
  )), , drop = FALSE]
  removed_n <- before_n - nrow(df)

  if (removed_n > 0) {
    message(
      "[plotting] Excluded ", removed_n,
      " N. gonorrhoeae rows from ", data_name, "."
    )
  }

  df
}

if (exclude_n_gonorrhoeae) {
  message("[plotting] AMR_EXCLUDE_N_GONORRHOEAE enabled.")
}

# --- 1. Define Names and Order ---
# Use the same mappings and order from your original script
classes <- c("J01A", "J01B", "J01C", "J01D", "J01E", "J01F", "J01G", "J01M")
antibiotic_names <- c(
  "Tetracyclines", "Glycopeptides and Lipopeptides", "Penicillins",
  "Non-Penicillin Beta-Lactams", "Sulfonamides and Trimethoprim",
  "Macrolides", "Aminoglycosides", "Quinolones"
)
atc_names_map <- setNames(antibiotic_names, classes)

antibiotic_list<- c("Quinolones", "Macrolides", "Aminoglycosides", "Non-Penicillin Beta-Lactams", 
                                     "Penicillins", "Sulfonamides and Trimethoprim", 
                                     "Tetracyclines")

# --- 2. Load and Prepare All Data Sources ---
cat("Loading data for Figure 1...\n")

# Load Pathogen-specific data
pathogen_data <- fread("Outputs/database_gradients_pathogen_ATC3_PCA_canonical_weighted_main.csv") %>%
  mutate(Antibiotic = atc_names_map[Antibiotic]) %>%
  filter(Antibiotic %in% antibiotic_list) %>%
  mutate(Response_fmt = sprintf("%.2f (%.2f, %.2f)", Response, Lower_CI, Upper_CI)) %>%
  filter_excluded_pathogen_rows("pathogen_data")

# Number of bug-drug pairs with significant positive relationships
num_significant_pairs <- pathogen_data %>%
  filter(Lower_CI > 0) %>%
  nrow()
cat(paste("Number of significant positive bug-drug pairs:", num_significant_pairs, "\n"))
total_pairs <- nrow(pathogen_data)

# Use lower.tail = FALSE to avoid the (1 - p) precision loss
p_value <- pbinom(num_significant_pairs - 1, total_pairs, 0.025, lower.tail = FALSE)

# Use format with a high digit count or sprintf for cleaner scientific notation
cat(paste("Binomial test p-value:", format(p_value, scientific = TRUE, digits = 10), "\n"))

# Load Drug-summary data (for the text labels)
summary_data <- fread("Outputs/database_gradients_ATC3_PCA_canonical_weighted_all.csv") %>%
  mutate(Antibiotic = atc_names_map[as.character(V1)]) %>%
  filter(Antibiotic %in% antibiotic_list)
summary_data_LowerCI <- fread("Outputs/database_lowerCI_ATC3_PCA_canonical_weighted_all.csv") %>%
  mutate(Antibiotic = atc_names_map[as.character(V1)]) %>%
  filter(Antibiotic %in% antibiotic_list)
summary_data_UpperCI <- fread("Outputs/database_upperCI_ATC3_PCA_canonical_weighted_all.csv") %>%
    mutate(Antibiotic = atc_names_map[as.character(V1)]) %>%
    filter(Antibiotic %in% antibiotic_list)
summary_data <- summary_data %>%
    left_join(summary_data_LowerCI %>% select(Antibiotic, x) %>% rename(Lower_CI = x), by = "Antibiotic") %>%
    left_join(summary_data_UpperCI %>% select(Antibiotic, x) %>% rename(Upper_CI = x), by = "Antibiotic") %>%
    mutate(Response_fmt = sprintf("%.2f (%.2f, %.2f)", x, Lower_CI, Upper_CI)) %>%
    rename(Response = x)
# Load Bootstrap gradient data (for the violins)
tryCatch({
  bootstrap_data <- fread("Outputs/database_gradients_bootstraps_ATC3_PCA_canonical_weighted_all.csv") %>%
    mutate(Antibiotic = atc_names_map[Antibiotic]) %>%
    filter(Antibiotic %in% antibiotic_list)
  cat("Bootstrap data loaded for Figure 1.\n")
}, error = function(e) {
  stop("Could not load 'Outputs/database_gradients_bootstraps_ATC3_PCA_canonical_weighted_all.csv'. 
       Please generate this file from your modeling script first.")
})

# print standard deviation of summary_data gradients with confidence intervals
cat("Standard deviation of summary drug gradients (for reference):\n")
sd_stats <- summary_data %>%
        summarise(
          SD = sd(Response),
          n = n(),
          SE_SD = sqrt(sum((Response - mean(Response))^2) / (n * (n - 1))) * sqrt(2 / (pi * (n - 1))),
          Lower_CI_SD = SD - 1.96 * SE_SD,
          Upper_CI_SD = SD + 1.96 * SE_SD
        ) %>%
        mutate(SD_fmt = sprintf("%.4f (95%% CI: %.4f - %.4f)", SD, Lower_CI_SD, Upper_CI_SD))
print(sd_stats %>% select(SD_fmt))


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
    aes(x = Gradient, y = Order, group = Order),
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
    aes(x = 4.1, y = Order, label = Response_fmt),
    hjust = 1, size = 8 / .pt, fontface = "plain"
  ) +
  geom_text(
    data = plot_df %>% filter(Type == "Total"),
    aes(x = 4.1, y = Order, label = Response_fmt),
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
  limits = c(-2.3, 4.1), 
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
pdf(file = "Figure1.pdf", width = 6.5, height = 9.3)
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
      aes(x = Gradient, y = Order, group = Order),
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
  
  # Save slide with PowerPoint dimensions (in slides output directory)
  filename <- file.path(SLIDES_DIR, paste0("Figure1_Slide", i, "_narrow_spillover_lagged.pdf"))
  pdf(file = filename, width = 8, height = 7.5)
  print(slide_plot)
  dev.off()
  
  cat(paste("Slide", i, "saved as", filename, "\n"))
}

cat("All PowerPoint slides generated successfully.\n")
cat("Figure 1 (Advanced ggplot2 version - FINAL Corrected & Refined) saved to 'Figure1.pdf'.\n\n")
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
.suppfig2b_path <- "Outputs/Nagorsen_gradients_pathogen_ATC3_PCA_canonical_hospital_to_all_filtered.csv"
if (!file.exists(.suppfig2b_path)) {
  message("[plotting] Supp Fig 2b: Nagorsen file not found — skipping Supplementary Figure 2b.")
} else {

# --- 1. Load and Prepare All Data Sources for Supp Fig 2b ---
cat("Loading data for Supp Fig 2b...\n")

# Load Pathogen-specific data (hospital)
pathogen_data_hosp <- fread(
  "Outputs/Nagorsen_gradients_pathogen_ATC3_PCA_canonical_hospital_to_all_filtered.csv"
) %>%
  filter(!is.na(Lower_CI) & !is.na(Upper_CI)) %>%
  mutate(
    Response_fmt = sprintf("%.2f (%.2f, %.2f)", Response, Lower_CI, Upper_CI)
  ) %>%
  filter_excluded_pathogen_rows("pathogen_data_hosp")

# Load Drug-summary data (hospital, for the text labels)
hosp_summary_grad <- setNames(
  as.vector(read.csv(
    "Outputs/Nagorsen_gradients_ATC3_PCA_canonical_hospital_to_all_filtered.csv"
  )[, 2]),
  read.csv(
    "Outputs/Nagorsen_gradients_ATC3_PCA_canonical_hospital_to_all_filtered.csv"
  )[, 1]
)
hosp_summary_lower <- setNames(
  as.vector(read.csv(
    "Outputs/Nagorsen_lowerCI_ATC3_PCA_canonical_hospital_to_all_filtered.csv"
  )[, 2]),
  read.csv(
    "Outputs/Nagorsen_lowerCI_ATC3_PCA_canonical_hospital_to_all_filtered.csv"
  )[, 1]
)
hosp_summary_upper <- setNames(
  as.vector(read.csv(
    "Outputs/Nagorsen_upperCI_ATC3_PCA_canonical_hospital_to_all_filtered.csv"
  )[, 2]),
  read.csv(
    "Outputs/Nagorsen_upperCI_ATC3_PCA_canonical_hospital_to_all_filtered.csv"
  )[, 1]
)

antibiotics_hosp <- unique(pathogen_data_hosp$Antibiotic)
summary_data_hosp <- data.frame(
  Antibiotic = antibiotics_hosp,
  Response = unname(hosp_summary_grad[antibiotics_hosp]),
  Lower_CI = unname(hosp_summary_lower[antibiotics_hosp]),
  Upper_CI = unname(hosp_summary_upper[antibiotics_hosp])
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
    "Outputs/Nagorsen_gradients_bootstraps_ATC3_PCA_canonical_hospital_to_all_filtered.csv"
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
  # Check if this drug class has any pathogen data
  path_subset_hosp <- pathogen_data_hosp %>%
    filter(Antibiotic == abx_code) %>%
    arrange(Pathogen)

  # Only include drug class if it has pathogen data
  if (nrow(path_subset_hosp) > 0) {
    abx_name <- atc_names_map[[abx_code]]
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
    path_subset_hosp$Plot_Label <- paste0("  ", path_subset_hosp$Pathogen)
    path_subset_hosp$Order <- current_order_hosp:(
      current_order_hosp + nrow(path_subset_hosp) - 1
    )
    path_subset_hosp$Type <- "Pathogen"
    path_subset_hosp$Antibiotic <- atc_names_map[path_subset_hosp$Antibiotic]

    plot_data_list_hosp[[length(plot_data_list_hosp) + 1]] <- path_subset_hosp
    current_order_hosp <- current_order_hosp + nrow(path_subset_hosp)

    # 3. Add "Total" Row - ONLY if random-effects model exists (in summary data)
    #    and there's more than one pathogen
    if (n_pathogens > 1 && abx_code %in% summary_data_hosp$Antibiotic) {
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
        aes(x = Gradient, y = Order, group = Order),
        fill = "grey",
        alpha = 0.8,
        trim = TRUE,
        width = 1
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
    aes(x = -1.0, y = Order, label = Plot_Label),
    hjust = 0,
    size = 8 / .pt,
    fontface = "bold"
  ) +
  geom_text(
    data = plot_df_hosp %>% filter(Type == "Pathogen"),
    aes(x = -1.0, y = Order, label = Plot_Label),
    hjust = 0,
    size = 8 / .pt,
    fontface = "italic"
  ) +
  geom_text(
    data = plot_df_hosp %>% filter(Type == "Total"),
    aes(x = -1.0, y = Order, label = Plot_Label),
    hjust = 0,
    size = 8 / .pt,
    fontface = "bold"
  ) +

  # --- RIGHT TEXT ---
  geom_text(
    data = plot_df_hosp %>% filter(Type == "Pathogen"),
    aes(x = 1.7, y = Order, label = Response_fmt),
    hjust = 1,
    size = 8 / .pt,
    fontface = "plain"
  ) +
  geom_text(
    data = plot_df_hosp %>% filter(Type == "Total"),
    aes(x = 1.7, y = Order, label = Response_fmt),
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
    limits = c(-1.0, 1.7),
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
  file = "Supplementary_Figure_4.pdf",
  width = 6.5,
  height = 6.5
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
        aes(x = Gradient, y = Order, group = Order),
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

# Save slide version with PowerPoint dimensions (in slides output directory)
pdf(
  file = file.path(SLIDES_DIR, "supp_fig2b_slide_narrow.pdf"),
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
} # end if Nagorsen file present
# -----------------------------------------------------------------------------
# End of Supplementary Figure 2b
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Figure 2: Meta-analysis of Pathogen Elasticities (VIOLIN PLOT Version)
# - Uses the 15-pathogen list from the user's provided image
# - CORRECTED: Loads *database* CSVs, matching Gemini_plotting.r
# -----------------------------------------------------------------------------
cat("Generating Figure 2 (Pathogen Violin Plot version, 15 pathogens)...\n")

# --- 1. Load and Prepare Data for Figure 2 ---
cat("Loading pathogen-level summary and bootstrap data for Figure 2...\n")

# Load Summary Data (for text labels)
tryCatch({
  # Load data to detect pathogens and order by elasticity
  grad_path <- fread(
    "Outputs/database_gradients_pathogen_PCA_canonical_weighted_all.csv",
    col.names = c("Pathogen_Code", "Response")
  )
  lower_path <- fread(
    "Outputs/database_lowerCI_pathogen_PCA_canonical_weighted_all.csv",
    col.names = c("Pathogen_Code", "Lower_CI")
  )
  upper_path <- fread(
    "Outputs/database_upperCI_pathogen_PCA_canonical_weighted_all.csv",
    col.names = c("Pathogen_Code", "Upper_CI")
  )

  # Merge into one data.table for text labels
  plot_data_fig2_summary <- grad_path %>%
    merge(lower_path, by = "Pathogen_Code") %>%
    merge(upper_path, by = "Pathogen_Code") %>%
    {
      if (exclude_n_gonorrhoeae) {
        .[!(tolower(trimws(as.character(Pathogen_Code))) %in% c(
          "n. gonorrhoeae",
          "neisseria gonorrhoeae"
        )), ]
      } else {
        .
      }
    } %>%
    arrange(desc(Response)) %>%
    mutate(
      Pathogen_Display = Pathogen_Code,
      Pathogen_Display = factor(Pathogen_Display,
                               levels = Pathogen_Code),
      Response_fmt = sprintf("%.2f (%.2f, %.2f)",
                            Response, Lower_CI, Upper_CI)
    )

  # Extract ordered pathogen names from data
  fig2_pathogen_codes <- plot_data_fig2_summary$Pathogen_Code
  fig2_pathogen_names <- plot_data_fig2_summary$Pathogen_Code

}, error = function(e) {
  stop("Could not load pathogen-level summary data.")
})

# print standard deviation of the Response variable with CI
response_sd <- sd(plot_data_fig2_summary$Response)
response_n <- length(plot_data_fig2_summary$Response)
response_se_sd <- response_sd / sqrt(2 * (response_n - 1))
response_ci_lower <- response_sd - 1.96 * response_se_sd
response_ci_upper <- response_sd + 1.96 * response_se_sd
cat("Standard deviation of pathogen Response values:",
    sprintf("%.2f (95%% CI: %.2f, %.2f)", response_sd,
            response_ci_lower, response_ci_upper), "\n")

# Load Bootstrap Data (for violins)
tryCatch({
  # This file path was already correct
  bootstrap_data_fig2 <- fread("Outputs/database_gradients_bootstraps_pathogen_PCA_canonical_weighted_all.csv") %>%
    mutate(Pathogen_Display = Pathogen) %>% 
    filter(Pathogen %in% fig2_pathogen_names) %>% 
    mutate(Pathogen_Display = factor(Pathogen_Display, 
                                    levels = rev(fig2_pathogen_names))) %>%
    filter_excluded_pathogen_rows("bootstrap_data_fig2")
  
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
  aes(x = Gradient, y = Pathogen_Display, group = Pathogen_Display),
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
  aes(x = Gradient, y = Pathogen_Display, group = Pathogen_Display),
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

# Save slide version with PowerPoint dimensions (in slides output directory)
pdf(file = file.path(SLIDES_DIR, "Figure2_Slide_narrow.pdf"), width = 8, height = 7.5)
print(p_fig2_slide)
dev.off()

cat("Figure 2 (Pathogen Violin Plot version, 15 pathogens) saved to 'Figure2.pdf'.\n")
cat("Figure 2 slide version saved to Outputs/slides/Figure2_Slide_narrow.pdf.\n\n")
# -----------------------------------------------------------------------------
# End of Figure 2 (Corrected Violin Version)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Supplementary Figure S1: Pathogen-level elasticities, HIC vs LMIC panels
# - Variation on main Figure 2: pathogen random-effects violin plot
# - Two panels: High-Income Countries and LMICs
# - Requires HIC and LMIC model outputs from make models-hic + make models-lmic
# -----------------------------------------------------------------------------
message("[plotting] Generating Supplementary Figure S2 (pathogen elasticities, HIC vs LMIC)...")
.suppfigs1_path <- "Outputs/database_gradients_pathogen_PCA_canonical_weighted_HIC.csv"
if (!file.exists(.suppfigs1_path)) {
  message("[plotting] Supp Fig S1: HIC pathogen model file not found — skipping Supplementary Figure S1.")
} else {

library(patchwork)

# Keep pathogen ordering consistent with main Figure 2.
.s1_pathogen_codes <- fig2_pathogen_codes
.s1_pathogen_names <- fig2_pathogen_names
.s1_pathogen_map <- setNames(.s1_pathogen_names, .s1_pathogen_codes)

.load_suppfigs1_fig2_data <- function(suffix) {
  grad <- fread(
    paste0("Outputs/database_gradients_pathogen_PCA_canonical_weighted_", suffix, ".csv"),
    col.names = c("Pathogen_Code", "Response")
  )
  lo <- fread(
    paste0("Outputs/database_lowerCI_pathogen_PCA_canonical_weighted_", suffix, ".csv"),
    col.names = c("Pathogen_Code", "Lower_CI")
  )
  hi <- fread(
    paste0("Outputs/database_upperCI_pathogen_PCA_canonical_weighted_", suffix, ".csv"),
    col.names = c("Pathogen_Code", "Upper_CI")
  )
  summary_df <- grad %>%
    merge(lo, by = "Pathogen_Code") %>%
    merge(hi, by = "Pathogen_Code") %>%
    filter(Pathogen_Code %in% .s1_pathogen_codes) %>%
    mutate(
      Pathogen_Display = .s1_pathogen_map[Pathogen_Code],
      Pathogen_Display = factor(Pathogen_Display, levels = rev(.s1_pathogen_names)),
      Response_fmt = sprintf("%.2f (%.2f, %.2f)", Response, Lower_CI, Upper_CI)
    )

  boot_df <- fread(
    paste0("Outputs/database_gradients_bootstraps_pathogen_PCA_canonical_weighted_", suffix, ".csv")
  ) %>%
    filter(Pathogen %in% .s1_pathogen_codes) %>%
    mutate(
      Pathogen_Display = .s1_pathogen_map[Pathogen],
      Pathogen_Display = factor(Pathogen_Display, levels = rev(.s1_pathogen_names))
    ) %>%
    filter_excluded_pathogen_rows(paste0("suppfig_s1_bootstrap_", suffix))

  list(summary = summary_df, bootstrap = boot_df)
}

.build_suppfigs1_fig2_panel <- function(summary_df, bootstrap_df, pathogen_levels, panel_title) {
  shading_df <- data.frame(Pathogen_Display = pathogen_levels) %>%
    mutate(ymin = row_number() - 0.5,
           ymax = row_number() + 0.5,
           is_shaded = row_number() %% 2 != 0) %>%
    filter(is_shaded)

  text_x <- 2.5
  plot_x_min <- -1.1

  ggplot() +
    geom_rect(
      data = shading_df,
      aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf, y = NULL),
      fill = "grey95",
      inherit.aes = FALSE
    ) +
    geom_violin(
      data = bootstrap_df,
      aes(x = Gradient, y = Pathogen_Display, group = Pathogen_Display),
      fill = "grey", alpha = 0.8, trim = TRUE, width = 2
    ) +
    geom_rect(
      data = summary_df,
      aes(ymin = as.numeric(Pathogen_Display) - 0.4,
          ymax = as.numeric(Pathogen_Display) + 0.4,
          xmin = 1.85, xmax = text_x),
      fill = ifelse(as.numeric(summary_df$Pathogen_Display) %% 2 != 0,
                    "grey95", "white"),
      alpha = 0.9
    ) +
    geom_point(
      data = summary_df,
      aes(x = Response, y = Pathogen_Display),
      shape = 18,
      size = if (slide_version) 3 else 2
    ) +
    geom_text(
      data = summary_df,
      aes(x = text_x, y = Pathogen_Display, label = Response_fmt),
      hjust = 1,
      size = if (slide_version) 4 else 2.5
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_x_continuous(
      "Elasticity",
      limits = c(plot_x_min, text_x),
      breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2)
    ) +
    scale_y_discrete(
      name = element_blank(),
      expand = expansion(mult = 0.01, add = 0.5)
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(face = "italic", size = if (slide_version) 14 else 9, hjust = 0),
      axis.text.x = element_text(size = if (slide_version) 14 else 9),
      axis.title.x = element_text(size = if (slide_version) 16 else 11),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = if (slide_version) 16 else 11)
    ) +
    ggtitle(panel_title)
}

message("[plotting] Supp Fig S2: loading random-effects HIC/LMIC data...")
s1_hic <- .load_suppfigs1_fig2_data("HIC")
s1_lmic <- .load_suppfigs1_fig2_data("LMIC")

# Restrict to like-for-like comparison: only pathogens estimated in both groups.
common_pathogen_codes <- intersect(s1_hic$summary$Pathogen_Code, s1_lmic$summary$Pathogen_Code)
common_pathogen_codes <- common_pathogen_codes[common_pathogen_codes %in% .s1_pathogen_codes]
common_pathogen_levels <- rev(.s1_pathogen_names[.s1_pathogen_names %in% .s1_pathogen_map[common_pathogen_codes]])

s1_hic$summary <- s1_hic$summary %>%
  filter(Pathogen_Code %in% common_pathogen_codes) %>%
  mutate(Pathogen_Display = factor(Pathogen_Display, levels = common_pathogen_levels))
s1_lmic$summary <- s1_lmic$summary %>%
  filter(Pathogen_Code %in% common_pathogen_codes) %>%
  mutate(Pathogen_Display = factor(Pathogen_Display, levels = common_pathogen_levels))
s1_hic$bootstrap <- s1_hic$bootstrap %>%
  filter(Pathogen %in% common_pathogen_codes) %>%
  mutate(Pathogen_Display = factor(Pathogen_Display, levels = common_pathogen_levels))
s1_lmic$bootstrap <- s1_lmic$bootstrap %>%
  filter(Pathogen %in% common_pathogen_codes) %>%
  mutate(Pathogen_Display = factor(Pathogen_Display, levels = common_pathogen_levels))

message("[plotting] Supp Fig S1: using ", length(common_pathogen_codes), " pathogens with estimates in both HIC and LMIC.")

p_suppfigs1_hic <- .build_suppfigs1_fig2_panel(
  s1_hic$summary, s1_hic$bootstrap, common_pathogen_levels, "High-Income Countries"
)
p_suppfigs1_lmic <- .build_suppfigs1_fig2_panel(
  s1_lmic$summary, s1_lmic$bootstrap, common_pathogen_levels, "Low- and Middle-Income Countries"
)

p_suppfigs1_combined <- p_suppfigs1_hic + p_suppfigs1_lmic + plot_layout(ncol = 1)

# --- 7. Save ---
pdf(file = file.path("Supplementary_Figure_S2.pdf"),
  width = 6.5, height = 6.5)
print(p_suppfigs1_combined)
dev.off()
message("[plotting] Supp Fig S2 saved to Supplementary_Figure_S2.pdf")

} # end if HIC pathogen file present
# -----------------------------------------------------------------------------
# End of Supplementary Figure S2
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Supplementary Figure 1: Separate panels for HICs and LMICs
# - Uses the same style as Figure 2 (violins)
# - Loads *database* CSVs, matching Gemini_plotting.r
# - NOW: Creates two separate panels side by side
# -----------------------------------------------------------------------------
cat("Generating Supplementary Figure 1 (HIC vs LMIC Separate Panels)...\n")
.suppfig1_hic_path <- "Outputs/database_gradients_ATC3_PCA_canonical_weighted_HIC.csv"
if (!file.exists(.suppfig1_hic_path)) {
  message("[plotting] Supp Fig 1: HIC model file not found — skipping Supplementary Figure 1.")
} else {

# Load required library for panel plots
library(patchwork)

# Toggle for slide version text size
slide_version <- FALSE  # Set to FALSE for regular version

# --- 1. Load and Prepare Data for Supp Fig 1 ---
cat("Loading pathogen-level summary and bootstrap data for Supp Fig 1...\n")
drug_class_gradients_hic <- setNames(
  as.vector(read.csv("Outputs/database_gradients_ATC3_PCA_canonical_weighted_HIC.csv", 
                     stringsAsFactors = FALSE)[, 2]), 
  read.csv("Outputs/database_gradients_ATC3_PCA_canonical_weighted_HIC.csv", 
           stringsAsFactors = FALSE)[, 1]
)
drug_class_lower_ci_hic <- setNames(
  as.vector(read.csv("Outputs/database_lowerCI_ATC3_PCA_canonical_weighted_HIC.csv", 
                     stringsAsFactors = FALSE)[, 2]),
  read.csv("Outputs/database_lowerCI_ATC3_PCA_canonical_weighted_HIC.csv", 
           stringsAsFactors = FALSE)[, 1]
)
drug_class_upper_ci_hic <- setNames(
  as.vector(read.csv("Outputs/database_upperCI_ATC3_PCA_canonical_weighted_HIC.csv", 
                     stringsAsFactors = FALSE)[, 2]),
  read.csv("Outputs/database_upperCI_ATC3_PCA_canonical_weighted_HIC.csv", 
           stringsAsFactors = FALSE)[, 1]
)
drug_class_bootstraps_hic <- fread(
  "Outputs/database_gradients_bootstraps_ATC3_PCA_canonical_weighted_HIC.csv"
)
drug_class_gradients_lmic <- setNames(
  as.vector(read.csv("Outputs/database_gradients_ATC3_PCA_canonical_weighted_LMIC.csv", 
                     stringsAsFactors = FALSE)[, 2]), 
  read.csv("Outputs/database_gradients_ATC3_PCA_canonical_weighted_LMIC.csv", 
           stringsAsFactors = FALSE)[, 1]
)
drug_class_lower_ci_lmic <- setNames(
  as.vector(read.csv("Outputs/database_lowerCI_ATC3_PCA_canonical_weighted_LMIC.csv", 
                     stringsAsFactors = FALSE)[, 2]),
  read.csv("Outputs/database_lowerCI_ATC3_PCA_canonical_weighted_LMIC.csv", 
           stringsAsFactors = FALSE)[, 1]
)
drug_class_upper_ci_lmic <- setNames(
  as.vector(read.csv("Outputs/database_upperCI_ATC3_PCA_canonical_weighted_LMIC.csv", 
                     stringsAsFactors = FALSE)[, 2]),
  read.csv("Outputs/database_upperCI_ATC3_PCA_canonical_weighted_LMIC.csv", 
           stringsAsFactors = FALSE)[, 1]
)
drug_class_bootstraps_lmic <- fread(
  "Outputs/database_gradients_bootstraps_ATC3_PCA_canonical_weighted_LMIC.csv"
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
    aes(x = Gradient, y = Antibiotic, group = Antibiotic),
    fill = "grey", alpha = 0.6, trim = TRUE,
    width = 2
  ) +
  
  # Add rectangles behind the text
  geom_rect(
    data = plot_data_fig_s1_summary_hic,
    aes(ymin = as.numeric(factor(Antibiotic, levels = rev(antibiotic_names))) - 0.4, 
        ymax = as.numeric(factor(Antibiotic, levels = rev(antibiotic_names))) + 0.4,
        xmin = 1.75, xmax = text_x_pos_fig_s1),
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
    aes(x = Gradient, y = Antibiotic, group = Antibiotic),
    fill = "grey", alpha = 0.6, trim = TRUE,
    width = 2
  ) +

  # Add rectangles behind the text
  geom_rect(
    data = plot_data_fig_s1_summary_hic,
    aes(ymin = as.numeric(factor(Antibiotic, levels = rev(antibiotic_names))) - 0.4, 
        ymax = as.numeric(factor(Antibiotic, levels = rev(antibiotic_names))) + 0.4,
        xmin = 1.75, xmax = text_x_pos_fig_s1),
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
# 8 by 7.5 inch pdf version (slide output directory)
pdf(file = file.path("Supplementary_Figure1.pdf"), width = 6.5, height = 6.5)
print(p_supp_fig1_combined)
dev.off()

cat("Supplementary Figure 1 (HIC vs LMIC Separate Panels) saved to Supplementary_Figure1.pdf.\n\n")
} # end if HIC file present
# -----------------------------------------------------------------------------
# End of Supplementary Figure 1
# -----------------------------------------------------------------------------
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggforce)
library(svglite)
library(ggtext)
library(scales)
source("utils.R")

# # main results
# results <- read.csv("Outputs/database_gradients_pathogen_ATC3_PCA_canonical60cutoff_no0resistance.csv", stringsAsFactors = FALSE)

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


# Figure 3 - global avertable burden
# Input paths declared in config.R AMR_CONFIG$burden_inputs for single-point editing
if (!exists("AMR_CONFIG")) source("config.R")
.fig3_paths <- c(AMR_CONFIG$burden_inputs$figure3_pathogen,
                  AMR_CONFIG$burden_inputs$figure3_optimistic,
                  AMR_CONFIG$burden_inputs$figure3_pessimistic)
.fig3_missing <- !all(file.exists(.fig3_paths))
if (.fig3_missing) {
  message("[plotting] Figure 3 burden input files not found — skipping Figure 3.")
  message("[plotting] Missing: ", paste(.fig3_paths[!file.exists(.fig3_paths)], collapse = ", "))
} else {
avertable_by_pathogen <- read.csv(AMR_CONFIG$burden_inputs$figure3_pathogen, stringsAsFactors = FALSE)
optimistic_by_pathogen <- read.csv(AMR_CONFIG$burden_inputs$figure3_optimistic, stringsAsFactors = FALSE)
pessimistic_by_pathogen <- read.csv(AMR_CONFIG$burden_inputs$figure3_pessimistic, stringsAsFactors = FALSE)

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
# print total avertable burden in main, optimistic, and pessimistic
# scenarios
message(
  "Total avertable burden in Main scenario: ",
  formatC(sum(avertable_by_pathogen$avertable_burden),
          big.mark = ",", format = "f", digits = 0),
  " deaths (95% CI: ",
  formatC(sum(avertable_by_pathogen$lower_bound),
          big.mark = ",", format = "f", digits = 0),
  " to ",
  formatC(sum(avertable_by_pathogen$upper_bound),
          big.mark = ",", format = "f", digits = 0),
  ")"
)
message(
  "Total avertable burden in Optimistic scenario: ",
  formatC(sum(optimistic_by_pathogen$avertable_burden),
          big.mark = ",", format = "f", digits = 0),
  " deaths (95% CI: ",
  formatC(sum(optimistic_by_pathogen$lower_bound),
          big.mark = ",", format = "f", digits = 0),
  " to ",
  formatC(sum(optimistic_by_pathogen$upper_bound),
          big.mark = ",", format = "f", digits = 0),
  ")"
)
message(
  "Total avertable burden in Pessimistic scenario: ",
  formatC(sum(pessimistic_by_pathogen$avertable_burden),
          big.mark = ",", format = "f", digits = 0),
  " deaths (95% CI: ",
  formatC(sum(pessimistic_by_pathogen$lower_bound),
          big.mark = ",", format = "f", digits = 0),
  " to ",
  formatC(sum(pessimistic_by_pathogen$upper_bound),
          big.mark = ",", format = "f", digits = 0),
  ")"
)

ggsave("Figure3.pdf", plot, width = 6.5, height = 9.3, units = "in")
# -----------------------------------------------------------------------------
# PowerPoint Slides Version of Figure 3 (Dynamic Ranges & Sorting)
# Slide outputs are written to Outputs/slides/ (separate from manuscript PDFs)
# -----------------------------------------------------------------------------
cat("Generating PowerPoint slides versions of Figure 3...\n")
slides_dir <- if (exists("AMR_CONFIG")) AMR_CONFIG$output_dirs$slides else "Outputs/slides"
dir.create(slides_dir, recursive = TRUE, showWarnings = FALSE)

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
  
  # Save slide with PowerPoint dimensions (in slides output directory)
  filename <- file.path(slides_dir, paste0("Figure3_", scen, "_Slide_narrow.pdf"))
  ggsave(filename, plot_slide, width = 8, height = 7.5, units = "in")
  
  cat(paste("Slide saved as", filename, "\n"))
}

cat("All dynamic PowerPoint slide versions for Figure 3 generated successfully.\n")
} # end if (!.fig3_missing)

# -----------------------------------------------------------------------------
# Supplementary Figure 5 - Avertable burden by drug and region (per 100k)
# -----------------------------------------------------------------------------
cat("Generating Supplementary Figure 5...\n")

# Load data
avertable_by_drug_region <- read.csv("Outputs/10pc_avertable_burden_by_drug_and_region_canonical_weighted_upper_region_main_overall.csv", stringsAsFactors = FALSE)

# Recode ATC3 class codes to full antibiotic class names for figure labels.
s5_drug_name_map <- c(
  "J01A" = "Tetracyclines",
  "J01B" = "Glycopeptides and Lipopeptides",
  "J01C" = "Penicillins",
  "J01D" = "Non-Penicillin Beta-Lactams",
  "J01E" = "Sulfonamides and Trimethoprim",
  "J01F" = "Macrolides",
  "J01G" = "Aminoglycosides",
  "J01M" = "Quinolones"
)

avertable_by_drug_region <- avertable_by_drug_region %>%
  mutate(drug = dplyr::recode(drug, !!!s5_drug_name_map, .default = drug))

# Filter out unwanted drug categories
avertable_by_drug_region <- avertable_by_drug_region %>%
  filter(!drug %in% c("Other", "J01X", "Multi-drug resistance in Salmonella Typhi and Paratyphi"))

# Helper function tailored for smaller rate numbers (preserves 2 decimal places)
format_burden_rate <- function(df) {
  df$burden_fmt <- paste0(
    formatC(df$avertable_burden_per_100k, format = "f", big.mark = ",", digits = 2), 
    " (", 
    formatC(df$lower_bound_per_100k, format = "f", big.mark = ",", digits = 2), 
    " to ", 
    formatC(df$upper_bound_per_100k, format = "f", big.mark = ",", digits = 2), 
    ")"
  )
  return(df)
}

# Apply formatting for rates
avertable_by_drug_region <- format_burden_rate(avertable_by_drug_region)

# Sort drugs based on total avertable rate overall to ensure consistent ordering across facets
drug_order <- avertable_by_drug_region %>% 
  group_by(drug) %>% 
  summarise(total_burden = sum(avertable_burden_per_100k, na.rm = TRUE)) %>% 
  arrange(total_burden) %>% 
  pull(drug)

avertable_by_drug_region$drug <- factor(avertable_by_drug_region$drug, levels = drug_order)

# Calculate text position for uniform formatting based on the rate upper bounds
max_x_s5 <- max(avertable_by_drug_region$upper_bound_per_100k, na.rm = TRUE)
avertable_by_drug_region$text_pos_x <- max_x_s5 + 0.05 * (max_x_s5 - min(avertable_by_drug_region$lower_bound_per_100k, na.rm = TRUE))

plot_s5 <- ggplot(avertable_by_drug_region, aes(x = avertable_burden_per_100k, y = drug)) +
  geom_bar(stat = "identity", fill = "grey50", color = "black", width = 0.7) +
  geom_errorbar(aes(xmin = lower_bound_per_100k, xmax = upper_bound_per_100k), width = 0.25, color = "black") +
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  facet_wrap(~ region, ncol = 1, scales = "fixed") +
  labs(x = "Deaths averted per 100,000", y = "") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.4),
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 10, family = "Helvetica"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10, family = "Helvetica"),
    axis.text.y = element_text(size = 10, family = "Helvetica"),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length.x = unit(4, "points"),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 9, family = "Helvetica", face = "bold", hjust = 0),
    strip.background = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "none"
  ) +
  coord_cartesian(clip = "off")

ggsave("Supplementary_Figure_5.pdf", plot_s5, width = 8, height = 10, units = "in")


# -----------------------------------------------------------------------------
# Supplementary Figure 6 - Avertable burden by pathogen and region (per 100k)
# -----------------------------------------------------------------------------
cat("Generating Supplementary Figure 6...\n")

# Load data 
avertable_by_pathogen_region <- read.csv("Outputs/10pc_avertable_burden_by_pathogen_and_region_canonical_weighted_upper_region_main_overall.csv", stringsAsFactors = FALSE)

# Sort pathogens strictly based on their overall avertable rate sum
pathogen_order <- avertable_by_pathogen_region %>% 
  group_by(pathogen) %>% 
  summarise(total_burden = sum(avertable_burden_per_100k, na.rm = TRUE)) %>% 
  arrange(total_burden) %>% 
  pull(pathogen)

avertable_by_pathogen_region$pathogen <- factor(avertable_by_pathogen_region$pathogen, levels = pathogen_order)

plot_s6 <- ggplot(avertable_by_pathogen_region, aes(x = avertable_burden_per_100k, y = pathogen)) +
  geom_bar(stat = "identity", fill = "grey50", color = "black", width = 0.7) +
  geom_errorbar(aes(xmin = lower_bound_per_100k, xmax = upper_bound_per_100k), width = 0.25, color = "black") +
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  # Use two columns and free x-scales (y-scales remain fixed/shared across rows)
  facet_wrap(~ region, ncol = 2, scales = "fixed") +
  labs(x = "Deaths averted per 100,000", y = "") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.4),
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 10, family = "Helvetica"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10, family = "Helvetica", angle = 45, hjust = 1), 
    
    # Italicize pathogen names
    axis.text.y = element_text(size = 10, family = "Helvetica", face = "italic"), 
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length.x = unit(4, "points"),
    axis.ticks.y = element_blank(),
    
    # Bold region names
    strip.text = element_text(size = 9, family = "Helvetica", face = "bold", hjust = 0),
    strip.background = element_blank(),
    panel.spacing.x = unit(1, "lines"), 
    panel.spacing.y = unit(0.5, "lines"),
    legend.position = "none",
    
    plot.margin = margin(10, 15, 10, 10, "points") 
  ) +
  coord_cartesian(clip = "off")

ggsave("Supplementary_Figure_6.pdf", plot_s6, width = 8, height = 12, units = "in")

# =============================================================================
# Figure 4 - Proportion of avertable burden by region, GDP, and antibiotic use
# Reads intermediate CSVs written by the burden stage (make burden).
# Input paths declared in config.R AMR_CONFIG$burden_inputs.
# =============================================================================
generate_figure4 <- function(
    avertable_region_path = AMR_CONFIG$burden_inputs$figure4_region,
    gdp_path              = AMR_CONFIG$burden_inputs$figure4_gdp,
    use_path              = AMR_CONFIG$burden_inputs$figure4_use,
    burden_region_path    = AMR_CONFIG$burden_inputs$lower_burden_region
) {
  library(metafor)

  avertable_by_region <- read.csv(avertable_region_path)
  gdp_by_lower_ihme_region <- read.csv(gdp_path)
  use_by_lower_ihme_region <- read.csv(use_path)

  # Both burden and GDP/use files are at the same lower-IHME-region level.
  # Merge directly on region name.
  avertable_by_region <- merge(avertable_by_region,
                               gdp_by_lower_ihme_region,
                               by.x = "region",
                               by.y = "lower_ihme_region")

  # Order regions by proportion avertable (descending)
  avertable_by_region <- avertable_by_region[
                         order(avertable_by_region$proportion_avertable,
                               decreasing = TRUE), ]

  # Merge with lower-region antibiotic use data
  avertable_by_region <- merge(avertable_by_region,
                               use_by_lower_ihme_region,
                               by.x = "region",
                               by.y = "lower_ihme_region")

  # Filter out regions with no use data (use_2018 == 0 or NA)
  avertable_by_region_with_use <- avertable_by_region[
                                  !is.na(avertable_by_region$use_2018) &
                                  avertable_by_region$use_2018 > 0, ]

# Meta-regression for totals vs GDP (linear and quadratic)
meta_totals_gdp_linear <- rma(yi = avertable_by_region$avertable_burden,
                vi = avertable_by_region$variance,
                mods = ~ gdp_2018,
                data = avertable_by_region)

meta_totals_gdp <- rma(yi = avertable_by_region$avertable_burden,
              vi = avertable_by_region$variance,
              mods = ~ gdp_2018 + I(gdp_2018^2),
              data = avertable_by_region)

# Meta-regression for totals vs use (linear and quadratic)
meta_totals_use_linear <- rma(yi = avertable_by_region_with_use$avertable_burden,
                vi = avertable_by_region_with_use$variance,
                mods = ~ use_2018,
                data = avertable_by_region_with_use)

meta_totals_use <- rma(yi = avertable_by_region_with_use$avertable_burden,
              vi = avertable_by_region_with_use$variance,
              mods = ~ use_2018 + I(use_2018^2),
              data = avertable_by_region_with_use)

# Meta-regression for proportions vs GDP (linear and quadratic)
meta_prop_gdp_linear <- rma(yi = avertable_by_region$proportion_avertable,
              vi = avertable_by_region$proportion_avertable_variance,
              mods = ~ gdp_2018,
              data = avertable_by_region)

meta_prop_gdp <- rma(yi = avertable_by_region$proportion_avertable,
            vi = avertable_by_region$proportion_avertable_variance,
            mods = ~ gdp_2018 + I(gdp_2018^2),
            data = avertable_by_region)

# Meta-regression for proportions vs use (linear and quadratic)
meta_prop_use_linear <- rma(yi = avertable_by_region_with_use$proportion_avertable,
              vi = avertable_by_region_with_use$
                proportion_avertable_variance,
              mods = ~ use_2018,
              data = avertable_by_region_with_use)

meta_prop_use <- rma(yi = avertable_by_region_with_use$proportion_avertable,
            vi = avertable_by_region_with_use$
              proportion_avertable_variance,
            mods = ~ use_2018 + I(use_2018^2),
            data = avertable_by_region_with_use)

# Meta-regression for per 100k vs GDP (linear and quadratic)
meta_per100k_gdp_linear <- rma(yi = avertable_by_region$avertable_burden_per_100k,
                  vi = avertable_by_region$variance_per_100k,
                  mods = ~ gdp_2018,
                  data = avertable_by_region)

meta_per100k_gdp <- rma(yi = avertable_by_region$avertable_burden_per_100k,
            vi = avertable_by_region$variance_per_100k,
            mods = ~ gdp_2018 + I(gdp_2018^2),
            data = avertable_by_region)

# Meta-regression for per 100k vs use (linear and quadratic)
meta_per100k_use_linear <- rma(
  yi = avertable_by_region_with_use$avertable_burden_per_100k,
  vi = avertable_by_region_with_use$variance_per_100k,
  mods = ~ use_2018,
  data = avertable_by_region_with_use)

meta_per100k_use <- rma(
  yi = avertable_by_region_with_use$avertable_burden_per_100k,
  vi = avertable_by_region_with_use$variance_per_100k,
  mods = ~ use_2018 + I(use_2018^2),
  data = avertable_by_region_with_use)

# Print BIC comparisons
cat("BIC Comparisons:\n")
cat("Totals vs GDP - Linear:", meta_totals_gdp_linear$fit.stats["BIC", "ML"], 
  "vs Quadratic:", meta_totals_gdp$fit.stats["BIC", "ML"], "\n")
cat("Totals vs Use - Linear:", meta_totals_use_linear$fit.stats["BIC", "ML"], 
  "vs Quadratic:", meta_totals_use$fit.stats["BIC", "ML"], "\n")
cat("Proportions vs GDP - Linear:", meta_prop_gdp_linear$fit.stats["BIC", "ML"], 
  "vs Quadratic:", meta_prop_gdp$fit.stats["BIC", "ML"], "\n")
cat("Proportions vs Use - Linear:", meta_prop_use_linear$fit.stats["BIC", "ML"], 
  "vs Quadratic:", meta_prop_use$fit.stats["BIC", "ML"], "\n")
cat("Per 100k vs GDP - Linear:", meta_per100k_gdp_linear$fit.stats["BIC", "ML"], 
  "vs Quadratic:", meta_per100k_gdp$fit.stats["BIC", "ML"], "\n")
cat("Per 100k vs Use - Linear:", meta_per100k_use_linear$fit.stats["BIC", "ML"], 
  "vs Quadratic:", meta_per100k_use$fit.stats["BIC", "ML"], "\n")

# Create prediction data for smooth curves
gdp_range <- seq(min(avertable_by_region$gdp_2018),
          max(avertable_by_region$gdp_2018),
          length.out = 100)
use_range <- seq(min(avertable_by_region_with_use$use_2018),
          max(avertable_by_region_with_use$use_2018),
          length.out = 100)

# Predictions for proportions vs GDP
pred_prop_gdp <- predict(meta_prop_gdp,
              newmods = cbind(gdp_range, gdp_range^2))
pred_df_prop_gdp <- data.frame(
  gdp = gdp_range,
  pred = pred_prop_gdp$pred,
  ci_lower = pred_prop_gdp$ci.lb,
  ci_upper = pred_prop_gdp$ci.ub
)

# Predictions for proportions vs use
pred_prop_use <- predict(meta_prop_use,
              newmods = cbind(use_range, use_range^2))
pred_df_prop_use <- data.frame(
  use = use_range,
  pred = pred_prop_use$pred,
  ci_lower = pred_prop_use$ci.lb,
  ci_upper = pred_prop_use$ci.ub
)

# Predictions for per 100k vs GDP
pred_per100k_gdp <- predict(meta_per100k_gdp,
              newmods = cbind(gdp_range, gdp_range^2))
pred_df_per100k_gdp <- data.frame(
  gdp = gdp_range,
  pred = pred_per100k_gdp$pred,
  ci_lower = pred_per100k_gdp$ci.lb,
  ci_upper = pred_per100k_gdp$ci.ub
)

# Predictions for per 100k vs use
pred_per100k_use <- predict(meta_per100k_use,
              newmods = cbind(use_range, use_range^2))
pred_df_per100k_use <- data.frame(
  use = use_range,
  pred = pred_per100k_use$pred,
  ci_lower = pred_per100k_use$ci.lb,
  ci_upper = pred_per100k_use$ci.ub
)

# Where do the predicted proportion curves peak?
gdp_seq_fine <- seq(100, 60000, by = 100)
pred_prop_gdp_fine <- predict(meta_prop_gdp,
                              newmods = cbind(gdp_seq_fine, gdp_seq_fine^2))
pred_df_prop_gdp_fine <- data.frame(
  gdp = gdp_seq_fine,
  pred = pred_prop_gdp_fine$pred,
  ci_lower = pred_prop_gdp_fine$ci.lb,
  ci_upper = pred_prop_gdp_fine$ci.ub
)
max_index_gdp <- which.max(pred_df_prop_gdp_fine$pred)
peak_gdp <- pred_df_prop_gdp_fine$gdp[max_index_gdp]
peak_proportion <- pred_df_prop_gdp_fine$pred[max_index_gdp]
print(paste0("Peak proportion avertable at GDP: ", peak_gdp,
             " with proportion: ", peak_proportion))
use_seq_fine <- seq(1, 30, by = 0.1)
pred_prop_use_fine <- predict(meta_prop_use,
                              newmods = cbind(use_seq_fine, use_seq_fine^2))
pred_df_prop_use_fine <- data.frame(
  use = use_seq_fine,
  pred = pred_prop_use_fine$pred,
  ci_lower = pred_prop_use_fine$ci.lb,
  ci_upper = pred_prop_use_fine$ci.ub
)
max_index_use <- which.max(pred_df_prop_use_fine$pred)
peak_use <- pred_df_prop_use_fine$use[max_index_use]
peak_proportion_use <- pred_df_prop_use_fine$pred[max_index_use]
print(paste0("Peak proportion avertable at use: ", peak_use,
             " with proportion: ", peak_proportion_use))

# Population-weighted average of proportion avertable for key regions
total_burden_by_region <- read.csv(burden_region_path)
# divide total_burden by population and multiply by 100k to get per 100k
total_burden_by_region$burden_per_100k <- (total_burden_by_region$total_burden / total_burden_by_region$population) * 100000
selected_regions <- c("Eastern Europe", "Southern Latin America",
            "South Asia", "North Africa and Middle East")
# selected_regions <- c("Central Sub-Saharan Africa", "Eastern Sub-Saharan Africa",
#                       "Southern Sub-Saharan Africa", "Western Sub-Saharan Africa")
pops <- numeric(length(selected_regions))
proportions <- numeric(length(selected_regions))
lower_bounds <- numeric(length(selected_regions))
upper_bounds <- numeric(length(selected_regions))
for (i in seq_along(selected_regions)) {
  region <- selected_regions[i]
  pops[i] <- total_burden_by_region[
  total_burden_by_region$region == region, "population"]
  proportions[i] <- avertable_by_region[
  avertable_by_region$region == region, "proportion_avertable"]
  lower_bounds[i] <- avertable_by_region[
  avertable_by_region$region == region, 
  "proportion_avertable_lower_bound"]
  upper_bounds[i] <- avertable_by_region[
  avertable_by_region$region == region, 
  "proportion_avertable_upper_bound"]
}
weighted_average_proportion <- sum(pops * proportions) / sum(pops)
weighted_average_lower <- sum(pops * lower_bounds) / sum(pops)
weighted_average_upper <- sum(pops * upper_bounds) / sum(pops)
print(paste0("Population-weighted average proportion avertable for ",
       "selected regions: ", round(weighted_average_proportion, 4), 
       " (", round(weighted_average_lower, 4), "-", 
       round(weighted_average_upper, 4), ")"))
for (i in seq_along(selected_regions)) {
  region <- selected_regions[i]
  gdp <- avertable_by_region[
    avertable_by_region$region == region, "gdp_2018"]
  print(paste0("GDP per capita for ", region, ": ", round(gdp, 2)))
}
print(avertable_by_region[c("region", "avertable_burden", "lower_bound", "upper_bound")])
print(avertable_by_region[c("region", "avertable_burden_per_100k", "lower_bound_per_100k", "upper_bound_per_100k")])

cat("Percent of total burden avertable:", sum(avertable_by_region$avertable_burden)/sum(avertable_by_region$total_burden), "(95% CI: ", sum(avertable_by_region$lower_bound)/sum(avertable_by_region$total_burden), " to ", sum(avertable_by_region$upper_bound)/sum(avertable_by_region$total_burden), ")\n")
# Precompute Panel A label text and x-position so spacing can be tuned explicitly
# for manuscript print layout (instead of inheriting slide spacing).
.panelA_xmax <- max(avertable_by_region$proportion_avertable_upper_bound * 100)
.panelA_xmin <- min(avertable_by_region$proportion_avertable_lower_bound * 100)
.panelA_range <- .panelA_xmax - .panelA_xmin
avertable_by_region$panelA_label <- paste0(
  format(round(avertable_by_region$proportion_avertable * 100, 2), nsmall = 2),
  "% (",
  format(round(avertable_by_region$proportion_avertable_lower_bound * 100, 2), nsmall = 2),
  "-",
  format(round(avertable_by_region$proportion_avertable_upper_bound * 100, 2), nsmall = 2),
  "%)"
)
avertable_by_region$panelA_label_x <- .panelA_xmax + 0.05 * .panelA_range
.panelA_x_upper <- max(avertable_by_region$panelA_label_x) - 0.02 * .panelA_range

# Panel A: Proportion avertable by region (bar chart)
proportions_by_region <- ggplot(avertable_by_region,
                 aes(x = proportion_avertable * 100,
                   y = reorder(region, 
                         proportion_avertable))) +
  geom_bar(stat = "identity", fill = "grey50") +
  geom_errorbar(aes(xmin = proportion_avertable_lower_bound * 100,
            xmax = proportion_avertable_upper_bound * 100), 
          width = 0.2) +
  geom_text(aes(label = panelA_label,
                x = panelA_label_x),
        hjust = 0, size = 2.8, family = "Helvetica",
        ) +
  labs(x = "Percentage of bacterial mortality avertible (%)", y = "Region", 
     tag = "A") +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.title.x = element_text(size = 10, family = "Helvetica"),
      axis.title.y = element_text(size = 10, family = "Helvetica"),
      axis.text.x = element_text(size = 8, family = "Helvetica"),
      axis.text.y = element_text(size = 8, family = "Helvetica"),
      plot.title = element_text(size = 12, family = "Helvetica"),
      axis.ticks = element_line(color = "black"),
      plot.tag = element_text(size = 12, face = "bold", family = "Helvetica"),
      text = element_text(family = "Helvetica")) +
  coord_cartesian(xlim = c(0, .panelA_x_upper), clip = "off") +
  theme(plot.margin = margin(5.5, 90, 5.5, 5.5, "points"))

proportions_vs_gdp <- ggplot(avertable_by_region,
              aes(x = gdp_2018, 
                y = proportion_avertable * 100)) +
  geom_ribbon(data = pred_df_prop_gdp, aes(x = gdp, ymin = ci_lower * 100,
                       ymax = ci_upper * 100),
        inherit.aes = FALSE, alpha = 0.2, fill = "black") +
  geom_line(data = pred_df_prop_gdp, aes(x = gdp, y = pred * 100),
        inherit.aes = FALSE, color = "black", linewidth = 0.7, alpha = 0.7) +
  geom_point(size = 2, color = "black") +
  geom_errorbar(aes(ymin = proportion_avertable_lower_bound * 100,
            ymax = proportion_avertable_upper_bound * 100), width = 0.2,
          color = "black") +
  labs(x = "GDP per capita (USD PPP)", y = "Percentage avertible (%)", tag = "B") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.title.x = element_text(size = 10, family = "Helvetica"),
      axis.title.y = element_text(size = 10, family = "Helvetica"),
      axis.text.x = element_text(size = 8, family = "Helvetica"),
      axis.text.y = element_text(size = 8, family = "Helvetica"),
      plot.title = element_text(size = 12, family = "Helvetica"),
      axis.ticks = element_line(color = "black"),
      plot.tag = element_text(size = 12, face = "bold", family = "Helvetica"),
      text = element_text(family = "Helvetica"))

proportions_vs_use <- ggplot(avertable_by_region_with_use,
              aes(x = use_2018, 
                y = proportion_avertable * 100)) +
  geom_ribbon(data = pred_df_prop_use, aes(x = use, ymin = ci_lower * 100,
                       ymax = ci_upper * 100),
        inherit.aes = FALSE, alpha = 0.2, fill = "black") +
  geom_line(data = pred_df_prop_use, aes(x = use, y = pred * 100),
        inherit.aes = FALSE, color = "black", linewidth = 0.7, alpha = 0.7) +
  geom_point(size = 2, color = "black") +
  geom_errorbar(aes(ymin = proportion_avertable_lower_bound * 100,
            ymax = proportion_avertable_upper_bound * 100), width = 0.2,
          color = "black") +
  labs(x = "Antibiotic use (DDD/1000 person-days)",
     y = "Percentage avertible (%)", tag = "C") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.title.x = element_text(size = 10, family = "Helvetica"),
      axis.title.y = element_text(size = 10, family = "Helvetica"),
      axis.text.x = element_text(size = 8, family = "Helvetica"),
      axis.text.y = element_text(size = 8, family = "Helvetica"),
      plot.title = element_text(size = 12, family = "Helvetica"),
      axis.ticks = element_line(color = "black"),
      plot.tag = element_text(size = 12, face = "bold", family = "Helvetica"),
      text = element_text(family = "Helvetica"))

# Panel D: Avertable per 100k vs GDP
per100k_vs_gdp <- ggplot(avertable_by_region,
            aes(x = gdp_2018, y = avertable_burden_per_100k)) +
  geom_ribbon(data = pred_df_per100k_gdp, aes(x = gdp, ymin = ci_lower,
                        ymax = ci_upper), 
        inherit.aes = FALSE, alpha = 0.2, fill = "black") +
  geom_line(data = pred_df_per100k_gdp, aes(x = gdp, y = pred),
        inherit.aes = FALSE, color = "black", linewidth = 0.7, alpha = 0.7) +
  geom_point(size = 2, color = "black") +
  geom_errorbar(aes(ymin = lower_bound_per_100k,
            ymax = upper_bound_per_100k), width = 0.2,
          color = "black") +
  labs(x = "GDP per capita (USD PPP)", y = "Avertible per 100k", tag = "D") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.title.x = element_text(size = 10, family = "Helvetica"),
      axis.title.y = element_text(size = 10, family = "Helvetica"),
      axis.text.x = element_text(size = 8, family = "Helvetica"),
      axis.text.y = element_text(size = 8, family = "Helvetica"),
      plot.title = element_text(size = 12, family = "Helvetica"),
      axis.ticks = element_line(color = "black"),
      plot.tag = element_text(size = 12, face = "bold", family = "Helvetica"),
      text = element_text(family = "Helvetica"))

# Panel E: Avertable per 100k vs antibiotic use
per100k_vs_use <- ggplot(avertable_by_region_with_use,
            aes(x = use_2018, y = avertable_burden_per_100k)) +
  geom_ribbon(data = pred_df_per100k_use, aes(x = use, ymin = ci_lower,
                        ymax = ci_upper), 
        inherit.aes = FALSE, alpha = 0.2, fill = "black") +
  geom_line(data = pred_df_per100k_use, aes(x = use, y = pred),
        inherit.aes = FALSE, color = "black", linewidth = 0.7, alpha = 0.7) +
  geom_point(size = 2, color = "black") +
  geom_errorbar(aes(ymin = lower_bound_per_100k,
            ymax = upper_bound_per_100k), width = 0.2,
          color = "black") +
  labs(x = "Antibiotic use (DDD/1000 person-days)",
     y = "Avertible per 100k", tag = "E") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.title.x = element_text(size = 10, family = "Helvetica"),
      axis.title.y = element_text(size = 10, family = "Helvetica"),
      axis.text.x = element_text(size = 8, family = "Helvetica"),
      axis.text.y = element_text(size = 8, family = "Helvetica"),
      plot.title = element_text(size = 12, family = "Helvetica"),
      axis.ticks = element_line(color = "black"),
      plot.tag = element_text(size = 12, face = "bold", family = "Helvetica"),
      text = element_text(family = "Helvetica"))

# Manuscript panels keep compact sizing consistent with other main figures.
panel_A <- proportions_by_region
panels_B_C <- gridExtra::arrangeGrob(
  proportions_vs_gdp,
  proportions_vs_use,
  ncol = 2, nrow = 1
)
panels_D_E <- gridExtra::arrangeGrob(
  per100k_vs_gdp,
  per100k_vs_use,
  ncol = 2, nrow = 1
)
panels_B_E <- gridExtra::arrangeGrob(
  proportions_vs_gdp,
  proportions_vs_use,
  per100k_vs_gdp,
  per100k_vs_use,
  ncol = 2, nrow = 2
)

# Slide variants use larger typography for projected presentation.
panel_A_slide <- proportions_by_region +
  labs(tag = "") +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(size = 24))
panels_B_C_slide <- gridExtra::arrangeGrob(
  proportions_vs_gdp + theme(axis.title.x = element_text(size = 16),
                            axis.title.y = element_text(size = 16),
                            axis.text.x = element_text(size = 14),
                            axis.text.y = element_text(size = 14)),
  proportions_vs_use + theme(axis.title.x = element_text(size = 16),
                            axis.title.y = element_text(size = 16),
                            axis.text.x = element_text(size = 14),
                            axis.text.y = element_text(size = 14)),
  ncol = 2, nrow = 1
)
panels_D_E_slide <- gridExtra::arrangeGrob(
  per100k_vs_gdp + theme(axis.title.x = element_text(size = 16),
                        axis.title.y = element_text(size = 16),
                        axis.text.x = element_text(size = 14),
                        axis.text.y = element_text(size = 14)),
  per100k_vs_use + theme(axis.title.x = element_text(size = 16),
                        axis.title.y = element_text(size = 16),
                        axis.text.x = element_text(size = 14),
                        axis.text.y = element_text(size = 14)),
  ncol = 2, nrow = 1
)
panels_B_E_slide <- gridExtra::arrangeGrob(
  proportions_vs_gdp + theme(axis.title.x = element_text(size = 16),
                            axis.title.y = element_text(size = 16),
                            axis.text.x = element_text(size = 14),
                            axis.text.y = element_text(size = 14)),
  proportions_vs_use + theme(axis.title.x = element_text(size = 16),
                            axis.title.y = element_text(size = 16),
                            axis.text.x = element_text(size = 14),
                            axis.text.y = element_text(size = 14)),
  per100k_vs_gdp + theme(axis.title.x = element_text(size = 16),
                        axis.title.y = element_text(size = 16),
                        axis.text.x = element_text(size = 14),
                        axis.text.y = element_text(size = 14)),
  per100k_vs_use + theme(axis.title.x = element_text(size = 16),
                        axis.title.y = element_text(size = 16),
                        axis.text.x = element_text(size = 14),
                        axis.text.y = element_text(size = 14)),
  ncol = 2, nrow = 2
)

# Full Figure 4: panel A on top, panels B-E in 2×2 grid below
figure4 <- gridExtra::grid.arrange(
  panel_A,
  panels_B_E,
  nrow = 2,
  heights = c(1.2, 1.8)
)

# Save manuscript Figure 4
ggsave("Figure4.pdf", figure4,
     width = 6.5, height = 7, units = "in")

# Save slide versions to slides output directory
.f4_slides_dir <- if (exists("AMR_CONFIG")) AMR_CONFIG$output_dirs$slides else "Outputs/slides"
dir.create(.f4_slides_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(.f4_slides_dir, "Figure4_Panels_B_C_narrow.pdf"), panels_B_C_slide,
       width = 9, height = 7.5, units = "in")
ggsave(file.path(.f4_slides_dir, "Figure4_Panels_D_E_narrow.pdf"), panels_D_E_slide,
       width = 9, height = 7.5, units = "in")
ggsave(file.path(.f4_slides_dir, "Figure4_Panels_B_E_narrow.pdf"), panels_B_E_slide,
  width = 9, height = 7.5, units = "in")
ggsave(file.path(.f4_slides_dir, "Figure4_Panel_A_narrow.pdf"), panel_A_slide,
       width = 9, height = 7.5, units = "in")
} # end generate_figure4()

# Invoke Figure 4 if all required input CSVs are present.
if (!exists("AMR_CONFIG")) source("config.R")
.fig4_paths <- c(
    AMR_CONFIG$burden_inputs$figure4_region,
    AMR_CONFIG$burden_inputs$figure4_gdp,
    AMR_CONFIG$burden_inputs$figure4_use,
    AMR_CONFIG$burden_inputs$lower_burden_region
)
.fig4_missing <- !all(file.exists(.fig4_paths))
if (.fig4_missing) {
  message("[plotting] Figure 4 burden input files not found — skipping Figure 4.")
  message("[plotting] Missing: ",
          paste(.fig4_paths[!file.exists(.fig4_paths)], collapse = ", "))
} else {
  .smoke_mode <- isTRUE(getOption("amr_smoke_mode", FALSE)) ||
                 identical(Sys.getenv("AMR_DEV_SMOKE", "0"), "1")
  if (.smoke_mode) {
    message("[plotting] Skipping Figure 4 in smoke mode (rma() requires adequate sample size).")
  } else {
    message("[plotting] Generating Figure 4...")
    generate_figure4()
  }
} # end if (!.fig4_missing)

