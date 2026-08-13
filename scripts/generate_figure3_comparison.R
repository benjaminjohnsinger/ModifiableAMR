# Generate a comparison Figure 3 using lower-region IHME for opt/pes scenarios.
# This is a one-off comparison script; it does not affect production outputs.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(forcats)
})

main <- read.csv("Outputs/10pc_avertable_burden_by_pathogen_canonical_weighted_lower_region_v2.csv",
                 stringsAsFactors = FALSE)
opt  <- read.csv("Outputs/10pc_avertable_burden_by_pathogen_lower_region_optimistic_comparison.csv",
                 stringsAsFactors = FALSE)
pes  <- read.csv("Outputs/10pc_avertable_burden_by_pathogen_lower_region_pessimistic_comparison.csv",
                 stringsAsFactors = FALSE)

cat("--- Lower-region opt/pes values ---\n")
for (p in c("E. coli", "Acinetobacter spp.", "K. pneumoniae")) {
  cat(sprintf("%-25s  main: %6d  opt: %6d  pes: %6d\n", p,
      round(main[main$pathogen == p, "avertable_burden"]),
      round(opt[opt$pathogen == p, "avertable_burden"]),
      round(pes[pes$pathogen == p, "avertable_burden"])))
}
cat("--- Upper-region opt/pes values (current production) ---\n")
opt2 <- read.csv("Outputs/10pc_avertable_burden_by_pathogen_canonical_weighted_upper_region_optimistic_overall.csv",
                 stringsAsFactors = FALSE)
pes2 <- read.csv("Outputs/10pc_avertable_burden_by_pathogen_canonical_weighted_upper_region_pessimistic_overall.csv",
                 stringsAsFactors = FALSE)
for (p in c("E. coli", "Acinetobacter spp.", "K. pneumoniae")) {
  cat(sprintf("%-25s  main: %6d  opt: %6d  pes: %6d\n", p,
      round(main[main$pathogen == p, "avertable_burden"]),
      round(opt2[opt2$pathogen == p, "avertable_burden"]),
      round(pes2[pes2$pathogen == p, "avertable_burden"])))
}

# ---------- Build the comparison figure ----------
main$Scenario <- "Main"
opt$Scenario  <- "Optimistic"
pes$Scenario  <- "Pessimistic"
all_data <- bind_rows(main, opt, pes)

pathogen_order <- main$pathogen[order(main$avertable_burden)]
all_data$pathogen <- factor(all_data$pathogen, levels = rev(pathogen_order))
all_data$Scenario <- factor(all_data$Scenario,
                            levels = c("Pessimistic", "Optimistic", "Main"))

max_x <- max(all_data$upper_bound, na.rm = TRUE)
all_data$text_pos_x <- max_x + 0.07 * (max_x - min(all_data$lower_bound, na.rm = TRUE))

fmt_n <- function(x) formatC(round(x, -floor(log10(abs(x))) + 1),
                              format = "f", big.mark = ",", digits = 0)
all_data$burden_fmt <- paste0(fmt_n(all_data$avertable_burden),
                               " (", fmt_n(all_data$lower_bound),
                               " to ", fmt_n(all_data$upper_bound), ")")

plot <- ggplot(all_data, aes(x = avertable_burden, y = Scenario, fill = Scenario)) +
  geom_bar(aes(color = Scenario), stat = "identity", width = 0.8) +
  geom_errorbar(aes(xmin = lower_bound, xmax = upper_bound, color = Scenario), width = 0.3) +
  geom_text(aes(x = text_pos_x, label = burden_fmt, color = Scenario,
                size = Scenario, fontface = Scenario), hjust = 0) +
  scale_fill_manual(values = c("Main" = "grey50", "Optimistic" = "white", "Pessimistic" = "white")) +
  scale_color_manual(values = c("Main" = "black", "Optimistic" = "gray30", "Pessimistic" = "gray30")) +
  scale_size_manual(values = c("Main" = 8 / .pt, "Optimistic" = 8 / .pt, "Pessimistic" = 8 / .pt)) +
  scale_discrete_manual(aesthetics = "fontface",
                        values = c("Main" = "bold", "Optimistic" = "plain", "Pessimistic" = "plain")) +
  scale_y_discrete(expand = expansion(add = c(0.1, 0.1))) +
  scale_x_continuous(
    labels = function(x) format(x, big.mark = ",", scientific = FALSE),
    breaks = seq(0, max(all_data$upper_bound, na.rm = TRUE) * 1.05, by = 10000)
  ) +
  facet_wrap(~ pathogen, ncol = 1, scales = "free_y") +
  labs(x = "Deaths averted", y = "",
       title = "Figure 3 comparison — lower-region IHME for opt/pes") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 10, family = "Helvetica"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10, family = "Helvetica"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length.x = unit(4, "points"),
    strip.text = element_text(size = 8, family = "Helvetica", face = "italic",
                              hjust = 0, margin = margin(t = 0, b = 1)),
    strip.background = element_blank(),
    panel.spacing = unit(0.3, "lines"),
    text = element_text(family = "Helvetica"),
    legend.position = "none",
    plot.margin = margin(2.5, 100, 5.5, 5.5, "points")
  ) +
  coord_cartesian(xlim = c(-100, max_x * 1.15), clip = "off")

out <- "Outputs/Figure3_lower_region_comparison.pdf"
ggsave(out, plot, width = 6.5, height = 9.3, units = "in")
cat("Saved:", out, "\n")
