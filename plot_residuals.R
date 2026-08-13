library(tidyverse)
library(ggh4x)

# Build residual scatter plots from per-combination files in Outputs/residuals.
residual_dir <- "Outputs/residuals"
residual_files <- list.files(
  residual_dir,
  pattern = "^residuals_.*\\.csv$",
  full.names = TRUE
)

if (length(residual_files) == 0) {
  stop("No residual files found in Outputs/residuals.")
}

atc_names <- c(
  J01A = "Tetracyclines",
  J01B = "Glycopeptides and Lipopeptides",
  J01C = "Penicillins",
  J01D = "Non-Penicillin Beta-Lactams",
  J01E = "Sulfonamides and Trimethoprim",
  J01F = "Macrolides",
  J01G = "Aminoglycosides",
  J01M = "Quinolones"
)

read_one_residual <- function(path) {
  base <- basename(path)
  parsed <- stringr::str_match(base, "^residuals_(J01[A-Z])_(.*)\\.csv$")
  if (is.na(parsed[1, 1])) {
    stop(paste0("Unexpected filename format: ", base))
  }

  antibiotic_code <- parsed[1, 2]
  pathogen_name <- parsed[1, 3]

  df <- read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("Residuals", "Consumption") %in% names(df))) {
    stop(paste0("Missing required columns in file: ", base))
  }

  df %>%
    mutate(
      Pathogen = pathogen_name,
      AntibioticCode = antibiotic_code,
      Antibiotic = antibiotic_code
    )
}

residual_df <- bind_rows(lapply(residual_files, read_one_residual)) %>%
  mutate(
    Pathogen = factor(Pathogen, levels = sort(unique(Pathogen))),
    AntibioticCode = factor(AntibioticCode, levels = names(atc_names)),
    Antibiotic = factor(Antibiotic, levels = names(atc_names))
  )

y_min <- -100
y_max <- 100

annotation_df <- residual_df %>%
  group_by(Pathogen, Antibiotic) %>%
  summarise(
    n_under = sum(Residuals < y_min, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_under > 0) %>%
  mutate(label = paste0("N. under -100: ", n_under))

pathogens_with_j01f <- residual_df %>%
  filter(AntibioticCode == "J01F") %>%
  pull(Pathogen) %>%
  unique() %>%
  as.character()

residual_with_j01f <- residual_df %>%
  filter(as.character(Pathogen) %in% pathogens_with_j01f) %>%
  droplevels()
residual_without_j01f <- residual_df %>%
  filter(!(as.character(Pathogen) %in% pathogens_with_j01f)) %>%
  droplevels()

annotation_with_j01f <- annotation_df %>%
  filter(as.character(Pathogen) %in% pathogens_with_j01f) %>%
  droplevels()
annotation_without_j01f <- annotation_df %>%
  filter(!(as.character(Pathogen) %in% pathogens_with_j01f)) %>%
  droplevels()

build_plot <- function(data_in, ann_in) {
  ggplot(data_in, aes(x = Consumption, y = Residuals)) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
    geom_point(alpha = 0.5, size = 0.5, color = "black") +
    facet_grid2(
      Pathogen ~ Antibiotic,
      scales = "free_x",
      drop = TRUE,
      independent = "x",
      render_empty = FALSE
    ) +
    geom_label(
      data = ann_in,
      aes(x = Inf, y = Inf, label = label),
      hjust = 1.05,
      vjust = 1.1,
      size = 2.7,
      color = "red",
      fill = scales::alpha("white", 0.7),
      linewidth = 0,
      inherit.aes = FALSE
    ) +
    scale_y_continuous(limits = c(y_min, y_max)) +
    labs(
      x = "log(Consumption)",
      y = "Residuals",
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold", size = 8),
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7)
    )
}

p_with_j01f <- build_plot(residual_with_j01f, annotation_with_j01f)
p_without_j01f <- build_plot(residual_without_j01f, annotation_without_j01f)

ggsave(
  filename = file.path(residual_dir, "residuals_vs_consumption_with_j01f.png"),
  plot = p_with_j01f,
  width = 8,
  height = 10,
  units = "in",
  dpi = 300
)

ggsave(
  filename = file.path(residual_dir, "residuals_vs_consumption_without_j01f.png"),
  plot = p_without_j01f,
  width = 8,
  height = 10,
  units = "in",
  dpi = 300
)

message("Saved residual plots to Outputs/residuals/.")