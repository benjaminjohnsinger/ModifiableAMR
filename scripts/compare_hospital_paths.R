library(dplyr)
source("utils.R")

legacy_prepare <- function() {
  data <- read.csv("Nagorsen_clean.csv", colClasses = c("units" = "character"), na.strings = c("NA"))
  data <- data[!is.na(data$amt_consumed), ]
  data <- data[!is.na(data$units), ]
  data <- data[data$amt_consumed < 10000, ]
  data <- data[!is.na(data$class_for_resistance), ]

  for (i in seq_len(nrow(data))) {
    data$pathogen[i] <- get_bacteria_name(data$pathogen[i])
  }

  for (atc_code in names(atc_mapping)) {
    data[data$class_for_resistance %in% atc_mapping[[atc_code]], "class_for_resistance"] <- atc_code
  }
  data <- data[!data$class_for_resistance %in% c("J01X", "Other"), ]

  df.pc <- read.csv("Chungman/Chungman_pca_renamed.csv")
  data$ISO3 <- iso3_ihme_mapping$iso3[match(data$country, iso3_ihme_mapping$country_name)]
  for (i in seq_len(nrow(data))) {
    year <- data$end_year[i]
    iso3 <- data$ISO3[i]
    if (iso3 %in% df.pc$ISO3) {
      row_index <- which(df.pc$ISO3 == iso3 & df.pc$Year == year)
      if (length(row_index) > 0) {
        data$PC1[i] <- df.pc$PC1[row_index[1]]
        data$PC2[i] <- df.pc$PC2[row_index[1]]
        data$PC3[i] <- df.pc$PC3[row_index[1]]
        data$GDP[i] <- df.pc$GDP[row_index[1]]
      }
    }
  }

  # Pre-refactor code called mutate() without assignment for these two conversions.
  data %>% filter(units == "DDD/100 bed days") %>% mutate(amt_consumed = amt_consumed / 10, units = "DDD/1000 bed days")
  data %>% filter(units == "DDD/1000 women/year") %>% mutate(amt_consumed = amt_consumed / 365, units = "DDD/1000 women/day")

  data[data$units == "DDD/inhabitants/day", "units"] <- "DDD/1000 inhabitants/day"
  data[data$units == "DDD/1000 inhabitants", "units"] <- "DDD/1000 inhabitants/day"

  data <- data[grepl("DDD", data$units) & grepl("1000", data$units) & grepl("day", data$units), ]
  data <- data[!grepl("community", data$ab_setting), ]

  data <- data %>%
    select(
      Consumption = amt_consumed,
      Resistance = percent_isolates_resistant,
      Pathogen = pathogen,
      DOI = doi,
      Antibiotic = class_for_resistance,
      PC1 = PC1,
      PC2 = PC2,
      PC3 = PC3,
      GDP = GDP,
      Year = end_year,
      Weight = 1,
      ISO3 = ISO3
    )

  pathogen_drug_counts <- table(paste(data$Pathogen, data$Antibiotic))
  pathogen_drug_to_remove <- names(pathogen_drug_counts[pathogen_drug_counts <= 20])
  data <- data[!paste(data$Pathogen, data$Antibiotic) %in% pathogen_drug_to_remove, ]
  data <- na.omit(data)

  data
}

current_prepare <- function() {
  read.csv("merged_data_Nagorsen_hospital_to_all_filtered.csv")
}

normalize_for_audit <- function(df) {
  numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  df[numeric_cols] <- lapply(df[numeric_cols], function(x) round(x, 10))
  df
}

scale_transform <- function(df) {
  global <- read.csv("antibiotic_consumption_by_ATC3.csv") %>%
    filter(Location == "Global", Year == "2018") %>%
    rename(Antibiotic = ATC.level.3.class, Global.Consumption = Antibiotic.consumption..DDD.1.000.day.) %>%
    select(Antibiotic, Global.Consumption)
  global$Antibiotic <- sub("-.*", "", global$Antibiotic)

  out <- df %>%
    left_join(global, by = c("Antibiotic")) %>%
    mutate(Consumption = Consumption / Global.Consumption) %>%
    select(-Global.Consumption)
  out <- na.omit(out)
  out$Consumption <- log(out$Consumption + 1)
  out$Resistance <- log(out$Resistance + 1)
  out$Weight <- out$Weight / max(out$Weight, na.rm = TRUE)
  out
}

fit_combined <- function(df) {
  abs <- sort(unique(df$Antibiotic))
  bugs <- sort(unique(df$Pathogen))
  rows <- list()

  for (a in abs) {
    for (p in bugs) {
      d <- df[df$Antibiotic == a & df$Pathogen == p, ]
      d <- d[complete.cases(d), ]
      d <- d[!is.infinite(d$Consumption) & !is.infinite(d$Resistance), ]
      if (nrow(d) <= 1) next

      m <- lm(Resistance ~ Consumption + PC1 + PC2 + PC3 + GDP + Year, data = d, weights = Weight)
      g <- coef(m)["Consumption"]
      if (is.na(g)) next

      rows[[length(rows) + 1]] <- data.frame(
        Antibiotic = a,
        Pathogen = p,
        Response = as.numeric(g)
      )
    }
  }

  if (length(rows) == 0) return(data.frame())
  bind_rows(rows)
}

legacy_raw <- legacy_prepare()
current_raw <- current_prepare()

legacy_raw <- normalize_for_audit(legacy_raw)
current_raw <- normalize_for_audit(current_raw)

legacy_keys <- sort(unique(paste(legacy_raw$Pathogen, legacy_raw$Antibiotic)))
current_keys <- sort(unique(paste(current_raw$Pathogen, current_raw$Antibiotic)))

cat("raw_rows legacy=", nrow(legacy_raw), " current=", nrow(current_raw), "\n", sep = "")
cat("raw_combos legacy=", length(legacy_keys), " current=", length(current_keys), "\n", sep = "")
cat("combos_only_legacy=", length(setdiff(legacy_keys, current_keys)), "\n", sep = "")
cat("combos_only_current=", length(setdiff(current_keys, legacy_keys)), "\n", sep = "")

legacy_counts <- legacy_raw %>% count(Pathogen, Antibiotic, name = "legacy_n")
current_counts <- current_raw %>% count(Pathogen, Antibiotic, name = "current_n")
count_diff <- full_join(legacy_counts, current_counts, by = c("Pathogen", "Antibiotic")) %>%
  mutate(
    legacy_n = ifelse(is.na(legacy_n), 0L, legacy_n),
    current_n = ifelse(is.na(current_n), 0L, current_n),
    delta_n = current_n - legacy_n
  ) %>%
  arrange(desc(abs(delta_n)))

cat("raw_count_diff_by_combo:\n")
print(count_diff)

audit_key_cols <- c("Pathogen", "Antibiotic", "Year", "ISO3", "DOI")
legacy_keys_df <- legacy_raw %>% distinct(across(all_of(audit_key_cols)))
current_keys_df <- current_raw %>% distinct(across(all_of(audit_key_cols)))

current_only_keys <- anti_join(current_keys_df, legacy_keys_df, by = audit_key_cols)
legacy_only_keys <- anti_join(legacy_keys_df, current_keys_df, by = audit_key_cols)

legacy_value_summary <- legacy_raw %>%
  group_by(across(all_of(audit_key_cols))) %>%
  summarise(
    legacy_n = n(),
    legacy_consumption_mean = mean(Consumption, na.rm = TRUE),
    legacy_resistance_mean = mean(Resistance, na.rm = TRUE),
    .groups = "drop"
  )

current_value_summary <- current_raw %>%
  group_by(across(all_of(audit_key_cols))) %>%
  summarise(
    current_n = n(),
    current_consumption_mean = mean(Consumption, na.rm = TRUE),
    current_resistance_mean = mean(Resistance, na.rm = TRUE),
    .groups = "drop"
  )

shared_key_value_changes <- inner_join(legacy_value_summary, current_value_summary, by = audit_key_cols) %>%
  mutate(
    delta_n = current_n - legacy_n,
    delta_consumption_mean = current_consumption_mean - legacy_consumption_mean,
    delta_resistance_mean = current_resistance_mean - legacy_resistance_mean
  ) %>%
  arrange(desc(abs(delta_consumption_mean) + abs(delta_resistance_mean) + abs(delta_n)))

combined_key_audit <- bind_rows(
  mutate(current_only_keys, audit_origin = "current_only"),
  mutate(legacy_only_keys, audit_origin = "legacy_only")
) %>%
  select(audit_origin, everything())

dir.create("Outputs", showWarnings = FALSE, recursive = TRUE)
write.csv(current_only_keys, "Outputs/hospital_path_row_audit_current_only.csv", row.names = FALSE)
write.csv(legacy_only_keys, "Outputs/hospital_path_row_audit_legacy_only.csv", row.names = FALSE)
write.csv(combined_key_audit, "Outputs/hospital_path_row_audit_combined.csv", row.names = FALSE)
write.csv(shared_key_value_changes, "Outputs/hospital_path_row_audit_shared_key_value_changes.csv", row.names = FALSE)

cat("row_audit_current_only_n=", nrow(current_only_keys), "\n", sep = "")
cat("row_audit_legacy_only_n=", nrow(legacy_only_keys), "\n", sep = "")
cat("row_audit_shared_keys_n=", nrow(shared_key_value_changes), "\n", sep = "")

if (length(setdiff(legacy_keys, current_keys)) > 0) {
  cat("only_legacy:\n")
  print(setdiff(legacy_keys, current_keys))
}
if (length(setdiff(current_keys, legacy_keys)) > 0) {
  cat("only_current:\n")
  print(setdiff(current_keys, legacy_keys))
}

legacy_model_df <- scale_transform(legacy_raw)
current_model_df <- scale_transform(current_raw)

cat("model_rows_after_scale legacy=", nrow(legacy_model_df), " current=", nrow(current_model_df), "\n", sep = "")

legacy_out <- fit_combined(legacy_model_df)
current_out <- fit_combined(current_model_df)

cat("grad_rows legacy=", nrow(legacy_out), " current=", nrow(current_out), "\n", sep = "")
cat("grad_mean legacy=", mean(legacy_out$Response), " current=", mean(current_out$Response), "\n", sep = "")

keys_legacy <- paste(legacy_out$Pathogen, legacy_out$Antibiotic)
keys_current <- paste(current_out$Pathogen, current_out$Antibiotic)
shared <- intersect(keys_legacy, keys_current)

if (length(shared) > 0) {
  l <- legacy_out[match(shared, keys_legacy), ]
  c <- current_out[match(shared, keys_current), ]
  delta <- c$Response - l$Response
  comp <- data.frame(
    Pathogen = c$Pathogen,
    Antibiotic = c$Antibiotic,
    legacy = l$Response,
    current = c$Response,
    delta = delta
  )
  comp <- comp[order(-abs(comp$delta)), ]

  cat("shared_grad_rows=", nrow(comp), "\n", sep = "")
  print(comp)
}
