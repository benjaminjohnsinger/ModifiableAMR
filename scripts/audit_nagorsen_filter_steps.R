library(dplyr)
source("utils.R")

raw <- read.csv("Nagorsen_clean.csv", colClasses = c("units" = "character"), na.strings = c("NA"))
raw$row_id <- seq_len(nrow(raw))

record_after <- function(df) sort(unique(df$row_id))
dropped_ids <- function(before_ids, after_ids) setdiff(before_ids, after_ids)

# ----- Legacy path (pre-refactor behavior) -----
legacy <- raw
legacy_ids <- list()
legacy_drop <- list()

legacy_ids$start <- record_after(legacy)

before <- record_after(legacy)
legacy <- legacy[!is.na(legacy$amt_consumed), ]
after <- record_after(legacy)
legacy_ids$na_amt <- after
legacy_drop$na_amt <- dropped_ids(before, after)

before <- record_after(legacy)
legacy <- legacy[!is.na(legacy$units), ]
after <- record_after(legacy)
legacy_ids$na_units <- after
legacy_drop$na_units <- dropped_ids(before, after)

before <- record_after(legacy)
legacy <- legacy[legacy$amt_consumed < 10000, ]
after <- record_after(legacy)
legacy_ids$amt_lt_10000 <- after
legacy_drop$amt_lt_10000 <- dropped_ids(before, after)

before <- record_after(legacy)
legacy <- legacy[!is.na(legacy$class_for_resistance), ]
after <- record_after(legacy)
legacy_ids$na_class <- after
legacy_drop$na_class <- dropped_ids(before, after)

legacy$pathogen <- vapply(legacy$pathogen, get_bacteria_name, character(1))
for (atc_code in names(atc_mapping)) {
  legacy[legacy$class_for_resistance %in% atc_mapping[[atc_code]], "class_for_resistance"] <- atc_code
}

before <- record_after(legacy)
legacy <- legacy[!legacy$class_for_resistance %in% c("J01X", "Other"), ]
after <- record_after(legacy)
legacy_ids$class_exclusion <- after
legacy_drop$class_exclusion <- dropped_ids(before, after)

# Legacy behavior: two mutate calls without assignment (no effect)
legacy %>% filter(units == "DDD/100 bed days") %>% mutate(amt_consumed = amt_consumed / 10, units = "DDD/1000 bed days")
legacy %>% filter(units == "DDD/1000 women/year") %>% mutate(amt_consumed = amt_consumed / 365, units = "DDD/1000 women/day")
legacy[legacy$units == "DDD/inhabitants/day", "units"] <- "DDD/1000 inhabitants/day"
legacy[legacy$units == "DDD/1000 inhabitants", "units"] <- "DDD/1000 inhabitants/day"

before <- record_after(legacy)
legacy <- legacy[grepl("DDD", legacy$units) & grepl("1000", legacy$units) & grepl("day", legacy$units), ]
after <- record_after(legacy)
legacy_ids$unit_pattern <- after
legacy_drop$unit_pattern <- dropped_ids(before, after)

before <- record_after(legacy)
legacy <- legacy[!grepl("community", legacy$ab_setting), ]
after <- record_after(legacy)
legacy_ids$community_filter <- after
legacy_drop$community_filter <- dropped_ids(before, after)

legacy$ISO3 <- iso3_ihme_mapping$iso3[match(legacy$country, iso3_ihme_mapping$country_name)]

df_pc <- read.csv("Chungman/Chungman_pca_renamed.csv")
idx <- match(paste(legacy$ISO3, legacy$end_year), paste(df_pc$ISO3, df_pc$Year))
legacy$PC1 <- df_pc$PC1[idx]
legacy$PC2 <- df_pc$PC2[idx]
legacy$PC3 <- df_pc$PC3[idx]
legacy$GDP <- df_pc$GDP[idx]

legacy <- legacy %>%
  select(
    row_id,
    Consumption = amt_consumed,
    Resistance = percent_isolates_resistant,
    Pathogen = pathogen,
    DOI = doi,
    Antibiotic = class_for_resistance,
    Weight = end_year,
    ISO3,
    PC1,
    PC2,
    PC3,
    GDP,
    Year = end_year
  )
legacy$Weight <- 1

before <- record_after(legacy)
combo_counts <- table(paste(legacy$Pathogen, legacy$Antibiotic))
combo_remove <- names(combo_counts[combo_counts <= 20])
legacy <- legacy[!paste(legacy$Pathogen, legacy$Antibiotic) %in% combo_remove, ]
after <- record_after(legacy)
legacy_ids$combo_min_20 <- after
legacy_drop$combo_min_20 <- dropped_ids(before, after)

before <- record_after(legacy)
legacy <- na.omit(legacy)
after <- record_after(legacy)
legacy_ids$na_omit <- after
legacy_drop$na_omit <- dropped_ids(before, after)

# ----- Current path (data_processing + hospital_nagorsen settings) -----
current <- raw
current_ids <- list()
current_drop <- list()

current_ids$start <- record_after(current)

before <- record_after(current)
current <- current[
  !is.na(current$amt_consumed) &
    !is.na(current$units) &
    !is.na(current$class_for_resistance) &
    !is.na(current$pathogen),
]
after <- record_after(current)
current_ids$na_core <- after
current_drop$na_core <- dropped_ids(before, after)

before <- record_after(current)
current <- current[current$amt_consumed < 10000, ]
after <- record_after(current)
current_ids$amt_lt_10000 <- after
current_drop$amt_lt_10000 <- dropped_ids(before, after)

current$pathogen <- vapply(current$pathogen, get_bacteria_name, character(1))
for (atc_code in names(atc_mapping)) {
  current[current$class_for_resistance %in% atc_mapping[[atc_code]], "class_for_resistance"] <- atc_code
}

before <- record_after(current)
current <- current[!current$class_for_resistance %in% c("J01X", "Other"), ]
after <- record_after(current)
current_ids$class_exclusion <- after
current_drop$class_exclusion <- dropped_ids(before, after)

# Current behavior: assignments are applied
current$amt_consumed[current$units == "DDD/100 bed days"] <- current$amt_consumed[current$units == "DDD/100 bed days"] / 10
current$units[current$units == "DDD/100 bed days"] <- "DDD/1000 bed days"
current$amt_consumed[current$units == "DDD/1000 women/year"] <- current$amt_consumed[current$units == "DDD/1000 women/year"] / 365
current$units[current$units == "DDD/1000 women/year"] <- "DDD/1000 women/day"
current$units[current$units == "DDD/inhabitants/day"] <- "DDD/1000 inhabitants/day"
current$units[current$units == "DDD/1000 inhabitants"] <- "DDD/1000 inhabitants/day"

before <- record_after(current)
current <- current[grepl("DDD", current$units) & grepl("1000", current$units) & grepl("day", current$units), ]
after <- record_after(current)
current_ids$unit_pattern <- after
current_drop$unit_pattern <- dropped_ids(before, after)

before <- record_after(current)
current <- current[!grepl("community", current$ab_setting), ]
after <- record_after(current)
current_ids$community_filter <- after
current_drop$community_filter <- dropped_ids(before, after)

current$ISO3 <- iso3_ihme_mapping$iso3[match(current$country, iso3_ihme_mapping$country_name)]
idx <- match(paste(current$ISO3, current$end_year), paste(df_pc$ISO3, df_pc$Year))
current$PC1 <- df_pc$PC1[idx]
current$PC2 <- df_pc$PC2[idx]
current$PC3 <- df_pc$PC3[idx]
current$GDP <- df_pc$GDP[idx]

current <- current %>%
  select(
    row_id,
    Consumption = amt_consumed,
    Resistance = percent_isolates_resistant,
    Pathogen = pathogen,
    DOI = doi,
    Antibiotic = class_for_resistance,
    Weight = end_year,
    ISO3,
    PC1,
    PC2,
    PC3,
    GDP,
    Year = end_year
  )
current$Weight <- 1

before <- record_after(current)
current <- current[complete.cases(current), ]
after <- record_after(current)
current_ids$complete_cases <- after
current_drop$complete_cases <- dropped_ids(before, after)

before <- record_after(current)
combo_counts <- table(paste(current$Pathogen, current$Antibiotic))
combo_remove <- names(combo_counts[combo_counts <= 20])
current <- current[!paste(current$Pathogen, current$Antibiotic) %in% combo_remove, ]
after <- record_after(current)
current_ids$combo_min_20 <- after
current_drop$combo_min_20 <- dropped_ids(before, after)

# ----- Compare semantically similar filter stages -----
step_map <- list(
  initial_missing = list(legacy = c("na_amt", "na_units", "na_class"), current = c("na_core")),
  amt_lt_10000 = list(legacy = c("amt_lt_10000"), current = c("amt_lt_10000")),
  class_exclusion = list(legacy = c("class_exclusion"), current = c("class_exclusion")),
  unit_pattern = list(legacy = c("unit_pattern"), current = c("unit_pattern")),
  community_filter = list(legacy = c("community_filter"), current = c("community_filter")),
  complete_cases = list(legacy = c("na_omit"), current = c("complete_cases")),
  combo_min_20 = list(legacy = c("combo_min_20"), current = c("combo_min_20"))
)

legacy_ref <- raw %>% select(row_id, country, pathogen, class_for_resistance, units, ab_setting, end_year, doi)
current_ref <- legacy_ref

summary_rows <- list()
detail_rows <- list()

for (step_name in names(step_map)) {
  legacy_step_ids <- unique(unlist(legacy_drop[step_map[[step_name]]$legacy]))
  current_step_ids <- unique(unlist(current_drop[step_map[[step_name]]$current]))

  legacy_only_ids <- setdiff(legacy_step_ids, current_step_ids)
  current_only_ids <- setdiff(current_step_ids, legacy_step_ids)
  overlap_ids <- intersect(legacy_step_ids, current_step_ids)

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    step = step_name,
    legacy_dropped_n = length(legacy_step_ids),
    current_dropped_n = length(current_step_ids),
    overlap_dropped_n = length(overlap_ids),
    legacy_only_dropped_n = length(legacy_only_ids),
    current_only_dropped_n = length(current_only_ids)
  )

  if (length(legacy_only_ids) > 0) {
    legacy_only_df <- legacy_ref %>%
      filter(row_id %in% legacy_only_ids) %>%
      mutate(step = step_name, dropped_by = "legacy_only")
    detail_rows[[length(detail_rows) + 1]] <- legacy_only_df
  }

  if (length(current_only_ids) > 0) {
    current_only_df <- current_ref %>%
      filter(row_id %in% current_only_ids) %>%
      mutate(step = step_name, dropped_by = "current_only")
    detail_rows[[length(detail_rows) + 1]] <- current_only_df
  }
}

summary_df <- bind_rows(summary_rows) %>% arrange(step)
detail_df <- if (length(detail_rows) > 0) bind_rows(detail_rows) else data.frame()

dir.create("Outputs", showWarnings = FALSE, recursive = TRUE)
write.csv(summary_df, "Outputs/nagorsen_filter_step_comparison_summary.csv", row.names = FALSE)
write.csv(detail_df, "Outputs/nagorsen_filter_step_comparison_details.csv", row.names = FALSE)

cat("Wrote Outputs/nagorsen_filter_step_comparison_summary.csv\n")
cat("Wrote Outputs/nagorsen_filter_step_comparison_details.csv\n")
print(summary_df)
