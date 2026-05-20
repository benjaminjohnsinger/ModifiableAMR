## BJS May 2025
## Code to merge datasets on AMR

library(dplyr)
library(tidyr)
source("utils.R")

prepare_nagorsen_hospital_regression_data <- function(
    nagorsen_path = "Nagorsen_clean.csv",
    pca_path = "Chungman/Chungman_pca_renamed.csv",
    output_path = "merged_data_Nagorsen_hospital_to_all_filtered.csv",
    sums_output_path = "merged_data_sums_Nagorsen_hospital_to_all_filtered.csv",
    min_entries_per_combo = 20
) {
    data <- read.csv(nagorsen_path, colClasses = c("units" = "character"), na.strings = c("NA"))

    data <- data[
        !is.na(data$amt_consumed) &
            !is.na(data$units) &
            !is.na(data$class_for_resistance) &
            !is.na(data$pathogen),
    ]
    data <- data[data$amt_consumed < 10000, ]

    data$pathogen <- vapply(data$pathogen, get_bacteria_name, character(1))
    for (atc_code in names(atc_mapping)) {
        data[data$class_for_resistance %in% atc_mapping[[atc_code]], "class_for_resistance"] <- atc_code
    }
    data <- data[!data$class_for_resistance %in% c("J01X", "Other"), ]

    data$amt_consumed[data$units == "DDD/100 bed days"] <- data$amt_consumed[data$units == "DDD/100 bed days"] / 10
    data$units[data$units == "DDD/100 bed days"] <- "DDD/1000 bed days"

    data$amt_consumed[data$units == "DDD/1000 women/year"] <- data$amt_consumed[data$units == "DDD/1000 women/year"] / 365
    data$units[data$units == "DDD/1000 women/year"] <- "DDD/1000 women/day"

    data$units[data$units == "DDD/inhabitants/day"] <- "DDD/1000 inhabitants/day"
    data$units[data$units == "DDD/1000 inhabitants"] <- "DDD/1000 inhabitants/day"

    data <- data[
        grepl("DDD", data$units) &
            grepl("1000", data$units) &
            grepl("day", data$units),
    ]

    # Hospital-to-all filtered analysis excludes rows explicitly labeled as community.
    data <- data[!grepl("community", data$ab_setting), ]

    data$ISO3 <- iso3_ihme_mapping$iso3[match(data$country, iso3_ihme_mapping$country_name)]

    df_pc <- read.csv(pca_path)
    idx <- match(paste(data$ISO3, data$end_year), paste(df_pc$ISO3, df_pc$Year))
    data$PC1 <- df_pc$PC1[idx]
    data$PC2 <- df_pc$PC2[idx]
    data$PC3 <- df_pc$PC3[idx]
    data$GDP <- df_pc$GDP[idx]

    data <- data %>%
        select(
            Consumption = amt_consumed,
            Resistance = percent_isolates_resistant,
            Pathogen = pathogen,
            DOI = doi,
            Antibiotic = class_for_resistance,
            Weight = end_year,
            ISO3 = ISO3,
            PC1 = PC1,
            PC2 = PC2,
            PC3 = PC3,
            GDP = GDP,
            Year = end_year
        )

    data$Weight <- 1
    data <- data[complete.cases(data), ]

    pathogen_drug_counts <- table(paste(data$Pathogen, data$Antibiotic))
    pathogen_drug_to_remove <- names(pathogen_drug_counts[pathogen_drug_counts <= min_entries_per_combo])
    data <- data[!paste(data$Pathogen, data$Antibiotic) %in% pathogen_drug_to_remove, ]

    write.csv(data, output_path, row.names = FALSE)

    data_sums <- data %>%
        group_by(Pathogen, Antibiotic) %>%
        summarise(Total.Isolates = n(), .groups = "drop")
    write.csv(data_sums, sums_output_path, row.names = FALSE)

    message("[data_processing] wrote ", output_path, " (rows=", nrow(data), ")")
    message("[data_processing] wrote ", sums_output_path, " (rows=", nrow(data_sums), ")")
}

prepare_main_regression_data <- function(
    joe_path = "pathogen_abx_analysis_all_variables_(class-specific).csv",
    atlas_path = "ATLAS_data/ATLAS_data_renamed.csv",
    atlas2_path = "ATLAS_more/ATLAS_more_renamed.csv",
    atlase_path = "ATLAS_Enterococcus/ATLAS_Enterococcus_renamed.csv",
    consumption_path = "antibiotic_consumption_by_ATC3.csv",
    pca_path = "Chungman/Chungman_pca_renamed.csv",
    no_covariates_path = "merged_data_new_no_covariates.csv",
    output_path = "merged_data_new.csv",
    sums_output_path = "merged_data_sums_new.csv",
    year_cutoff = 2018
) {
  ## Merge Joe's data with ATLAS
  JOE <- read.csv(joe_path)
  ATLAS <- read.csv(atlas_path)
  ATLAS2 <- read.csv(atlas2_path)
  ATLASE <- read.csv(atlase_path)
  # concatenate the three ATLAS datasets
  ATLAS <- rbind(ATLAS, ATLAS2, ATLASE)
  # remove duplicated rows
  ATLAS <- ATLAS[!duplicated(ATLAS), ]

  JOE <- JOE[,c("ISO3", "Year", "Pathogen", "ATC.Class", "Percent.Resistant.Isolates", "Total.Isolates")]
  ATLAS <- ATLAS[,c("ISO3", "Year", "Pathogen", "ATC.Class", "Percent.Resistant.Isolates", "Total.Isolates")]

  # Merge datasets
  merged_data <- rbind(JOE, ATLAS)
  # Coerce to numeric in case any source file reads them as character
  merged_data$Total.Isolates <- as.numeric(merged_data$Total.Isolates)
  merged_data$Percent.Resistant.Isolates <- as.numeric(merged_data$Percent.Resistant.Isolates)
  message("[data_processing] merged ", nrow(merged_data), " observations: ",
          length(unique(merged_data$Pathogen)), " pathogens, ",
          length(unique(merged_data$ISO3)), " countries, years ",
          min(merged_data$Year, na.rm = TRUE), "--", max(merged_data$Year, na.rm = TRUE))

  # Remove rows with Total, Other, or J01X ATC.Class
  merged_data <- merged_data[!merged_data$ATC.Class %in% c("Total", "Other", "J01X"), ]
  # No years after year_cutoff
  merged_data <- merged_data[merged_data$Year <= year_cutoff, ]
  # Remove HKG
  merged_data <- merged_data[merged_data$ISO3 != "HKG", ]

  # Load GRAM consumption data
  consumption <- read.csv(consumption_path)
  # Extract the first 4 characters of ATC.level.3.class to get the ATC class
  consumption$ATC.level.3.class <- substr(consumption$ATC.level.3.class, 1, 4)
  # rename columns ("Location","Year","ATC level 3 class","Antibiotic consumption (DDD/1,000/day)") in consumption data and map country names onto ISO3
  consumption <- consumption %>%
    rename(ISO3 = Location, Year = Year, Antibiotic = ATC.level.3.class, Consumption = Antibiotic.consumption..DDD.1.000.day.) %>%
    mutate(ISO3 = iso3_ihme_mapping$iso3[match(ISO3, iso3_ihme_mapping$country_name)])

  consumption <- consumption %>%
    group_by(ISO3, Year, Antibiotic) %>%
    summarise(Consumption = sum(Consumption, na.rm = TRUE),
              .groups = "drop")

  # Safely merge the aggregated data
  consumption$Year <- as.character(consumption$Year)
  merged_data <- merged_data %>%
    left_join(consumption,
              by = c("ISO3" = "ISO3", "Year" = "Year",
                     "ATC.Class" = "Antibiotic")) %>%
    rename(Antibiotic.Consumption = Consumption)

  # 3. Create a wide format of consumption for ALL ATC classes
  wide_consumption <- consumption %>%
    pivot_wider(
      names_from = Antibiotic,
      values_from = Consumption,
      names_glue = "{Antibiotic}.Consumption"
    )

  # 4. Join the wide data and replace the redundant ATC class values with NA
  merged_data <- merged_data %>%
    left_join(wide_consumption, by = c("ISO3", "Year")) %>%
    mutate(
      across(
        # Target only the newly created *.Consumption columns (ignoring the primary one)
        ends_with(".Consumption") & !matches("^Antibiotic\\.Consumption$"),
        # If the column name matches the row's ATC.Class, replace with NA; otherwise, keep the value
        ~ ifelse(paste0(ATC.Class, ".Consumption") == cur_column(), NA_real_, .x)
      )
    )
  
  write.csv(merged_data, no_covariates_path, row.names = FALSE)

  # Load PCA covariates and merge using vectorized join
  df.pc <- read.csv(pca_path)
  # Create composite key for PCA: ISO3|Year
  merged_data$key_pca <- paste(merged_data$ISO3, merged_data$Year, sep = "|")
  df.pc$key_pca <- paste(df.pc$ISO3, df.pc$Year, sep = "|")

  idx_pca <- match(merged_data$key_pca, df.pc$key_pca)
  merged_data$PC1 <- df.pc$PC1[idx_pca]
  merged_data$PC2 <- df.pc$PC2[idx_pca]
  merged_data$PC3 <- df.pc$PC3[idx_pca]
  merged_data$GDP <- df.pc$GDP[idx_pca]
  merged_data$key_pca <- NULL  # Clean up temporary key
  df.pc$key_pca <- NULL

  # Save the merged data with PCA covariates
  write.csv(merged_data, output_path, row.names = FALSE)
  message("[data_processing] wrote ", output_path, " (rows=", nrow(merged_data), ", ",
          sum(!is.na(merged_data$Antibiotic.Consumption)), " rows with consumption, ",
          sum(complete.cases(merged_data)), " complete cases)")

  merged_data_sums <- merged_data %>%
    group_by(Pathogen, ATC.Class) %>%
    summarise(
      Total.Isolates = sum(Total.Isolates, na.rm = TRUE),
      Percent.Resistant.Isolates = mean(Percent.Resistant.Isolates, na.rm = TRUE),
      Antibiotic.Consumption = mean(Antibiotic.Consumption, na.rm = TRUE),
      PC1 = mean(PC1, na.rm = TRUE),
      PC2 = mean(PC2, na.rm = TRUE),
      PC3 = mean(PC3, na.rm = TRUE),
      GDP = mean(GDP, na.rm = TRUE)
    ) %>%
    ungroup()
  # save the merged data with sums
  write.csv(merged_data_sums, sums_output_path, row.names = FALSE)
  message("[data_processing] wrote ", sums_output_path, " (",
          nrow(merged_data_sums), " pathogen-class combinations)")
}
