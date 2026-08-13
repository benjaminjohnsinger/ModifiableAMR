## BJS May 2025
## Code to merge datasets on AMR

library(dplyr)
library(tidyr)
source("utils.R")

prepare_nagorsen_hospital_regression_data <- function(
    nagorsen_path = "Nagorsen_clean.csv",
    pca_path = "Chungman/pcato10.csv",
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
    data$PC4 <- df_pc$PC4[idx]
    data$PC5 <- df_pc$PC5[idx]
    data$PC6 <- df_pc$PC6[idx]
    data$PC7 <- df_pc$PC7[idx]
    data$PC8 <- df_pc$PC8[idx]
    data$PC9 <- df_pc$PC9[idx]
    data$PC10 <- df_pc$PC10[idx]
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
            PC4 = PC4,
            PC5 = PC5,
            PC6 = PC6,
            PC7 = PC7,
            PC8 = PC8,
            PC9 = PC9,
            PC10 = PC10,
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
    pca_path = "Chungman/pcato10.csv",
    no_covariates_path = "summed_data_new_no_covariates.csv",
    output_path = "summed_data_new.csv",
    sums_output_path = "summed_data_sums_new.csv",
    year_cutoff = 2018,
    class_mode = c("atc_class", "specific_abs")
) {
  class_mode <- match.arg(class_mode)

  load_consumption_data <- function(path) {
    consumption_raw <- read.csv(path, stringsAsFactors = FALSE)

    # Standard global-use format (DDD per 1,000/day by country name).
    if (all(c("Location", "Year", "ATC.level.3.class", "Antibiotic.consumption..DDD.1.000.day.") %in% names(consumption_raw))) {
      consumption_raw$ATC.level.3.class <- substr(consumption_raw$ATC.level.3.class, 1, 4)
      consumption <- consumption_raw %>%
        rename(
          ISO3 = Location,
          Antibiotic = ATC.level.3.class,
          Consumption = Antibiotic.consumption..DDD.1.000.day.
        ) %>%
        mutate(
          ISO3 = iso3_ihme_mapping$iso3[match(ISO3, iso3_ihme_mapping$country_name)],
          Year = as.character(Year)
        )
    } else if (all(c("ISO3", "Year", "Antibiotic", "Consumption") %in% names(consumption_raw))) {
      # Already harmonized format.
      consumption <- consumption_raw %>%
        mutate(
          Antibiotic = substr(Antibiotic, 1, 4),
          Year = as.character(Year)
        )
    } else if (all(c("ISO3", "Year", "Antimicrobial", "DDD") %in% names(consumption_raw))) {
      # Raw IQVIA format: map antimicrobial names to ATC classes then aggregate.
      consumption_raw$Antibiotic <- vapply(
        consumption_raw$Antimicrobial,
        function(x) {
          class <- get_atc_class(x)
          if (is.na(class)) NA_character_ else class
        },
        character(1)
      )
      unmapped <- sort(unique(consumption_raw$Antimicrobial[is.na(consumption_raw$Antibiotic)]))
      if (length(unmapped) > 0) {
        message("[data_processing] raw IQVIA antimicrobials not mapped to ATC and excluded: ",
                paste(unmapped, collapse = ", "))
      }

      consumption <- consumption_raw %>%
        filter(!is.na(Antibiotic)) %>%
        mutate(Year = as.character(Year)) %>%
        group_by(ISO3, Year, Antibiotic) %>%
        summarise(Consumption = sum(DDD, na.rm = TRUE), .groups = "drop")
    } else {
      stop(
        "[data_processing] unsupported consumption schema in ", path,
        ". Expected one of: ",
        "(Location, Year, ATC.level.3.class, Antibiotic.consumption..DDD.1.000.day.), ",
        "(ISO3, Year, Antibiotic, Consumption), or ",
        "(ISO3, Year, Antimicrobial, DDD).",
        call. = FALSE
      )
    }

    consumption %>%
      group_by(ISO3, Year, Antibiotic) %>%
      summarise(Consumption = sum(Consumption, na.rm = TRUE), .groups = "drop")
  }

  ## Merge Joe's data with ATLAS and GASP data
  JOE <- read.csv(joe_path)
  ATLAS <- read.csv(atlas_path)
  ATLAS2 <- read.csv(atlas2_path)
  ATLASE <- read.csv(atlase_path)
  # concatenate the three ATLAS datasets
  ATLAS <- rbind(ATLAS, ATLAS2, ATLASE)
  # remove duplicated rows
  ATLAS <- ATLAS[!duplicated(ATLAS), ]

  if (class_mode == "specific_abs") {
    map_antibiotic_group <- function(x) {
      x <- trimws(as.character(x))
      out <- rep(NA_character_, length(x))
      mapping_names <- names(antibiotic_mapping)

      # Joe data typically already has standardized values at mapping-name level.
      idx_name <- match(x, mapping_names)
      out[!is.na(idx_name)] <- x[!is.na(idx_name)]

      # ATLAS data contains specific drugs; map these to standardized groups.
      for (group_name in mapping_names) {
        mapped_values <- antibiotic_mapping[[group_name]]
        if (length(mapped_values) == 0) {
          next
        }
        idx <- is.na(out) & x %in% mapped_values
        out[idx] <- group_name
      }

      out
    }

    JOE <- JOE[, c("ISO3", "Year", "Pathogen", "ATC.Class", "Antibiotic_standardized", "Percent.Resistant.Isolates", "Total.Isolates")]
    ATLAS <- ATLAS[, c("ISO3", "Year", "Pathogen", "ATC.Class", "Antibiotic", "Percent.Resistant.Isolates", "Total.Isolates")]

    JOE$Antibiotic <- map_antibiotic_group(JOE$Antibiotic_standardized)
    ATLAS$Antibiotic <- map_antibiotic_group(ATLAS$Antibiotic)

    joe_unmapped <- sum(is.na(JOE$Antibiotic))
    atlas_unmapped <- sum(is.na(ATLAS$Antibiotic))
    if (joe_unmapped > 0 || atlas_unmapped > 0) {
      message("[data_processing] specific_abs unmapped antibiotics removed: Joe=", joe_unmapped,
              ", ATLAS=", atlas_unmapped)
    }

    JOE <- JOE[!is.na(JOE$Antibiotic), c("ISO3", "Year", "Pathogen", "ATC.Class", "Antibiotic", "Percent.Resistant.Isolates", "Total.Isolates")]
    ATLAS <- ATLAS[!is.na(ATLAS$Antibiotic), c("ISO3", "Year", "Pathogen", "ATC.Class", "Antibiotic", "Percent.Resistant.Isolates", "Total.Isolates")]
  } else {
    JOE <- JOE[, c("ISO3", "Year", "Pathogen", "ATC.Class", "Percent.Resistant.Isolates", "Total.Isolates")]
    ATLAS <- ATLAS[, c("ISO3", "Year", "Pathogen", "ATC.Class", "Percent.Resistant.Isolates", "Total.Isolates")]
  }

  # Coerce to numeric in case any source file reads them as character
  JOE$Percent.Resistant.Isolates <- as.numeric(JOE$Percent.Resistant.Isolates)
  JOE$Total.Isolates <- as.numeric(JOE$Total.Isolates)
  ATLAS$Percent.Resistant.Isolates <- as.numeric(ATLAS$Percent.Resistant.Isolates)
  ATLAS$Total.Isolates <- as.numeric(ATLAS$Total.Isolates)

  # within each dataset, if there are rows with the same ISO3, Year, Pathogen, class and Total.Isolates,
  # combine Percent.Resistant.Isolates with independent-probabilities formula.
  if (class_mode == "specific_abs") {
    JOE <- JOE %>%
      group_by(ISO3, Year, Pathogen, ATC.Class, Antibiotic, Total.Isolates) %>%
      summarise(
        Percent.Resistant.Isolates = (1 - prod(1 - Percent.Resistant.Isolates / 100)) * 100,
        .groups = "drop"
      )
    ATLAS <- ATLAS %>%
      group_by(ISO3, Year, Pathogen, ATC.Class, Antibiotic, Total.Isolates) %>%
      summarise(
        Percent.Resistant.Isolates = (1 - prod(1 - Percent.Resistant.Isolates / 100)) * 100,
        .groups = "drop"
      )
  } else {
    JOE <- JOE %>%
      group_by(ISO3, Year, Pathogen, ATC.Class, Total.Isolates) %>%
      summarise(
        Percent.Resistant.Isolates = (1 - prod(1 - Percent.Resistant.Isolates / 100)) * 100,
        .groups = "drop"
      )
    ATLAS <- ATLAS %>%
      group_by(ISO3, Year, Pathogen, ATC.Class, Total.Isolates) %>%
      summarise(
        Percent.Resistant.Isolates = (1 - prod(1 - Percent.Resistant.Isolates / 100)) * 100,
        .groups = "drop"
      )
  }

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

  # Remove rows with N. gonorrheae - not subject to same bystander exposures as other pathogens
  merged_data <- merged_data[!grepl("gonorrhoeae", tolower(merged_data$Pathogen)), ]

  # Load consumption data from standard ATC3 format, harmonized format, or raw IQVIA DDD format.
  consumption <- load_consumption_data(consumption_path)

  # Safely merge the aggregated data
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
  merged_data$PC4 <- df.pc$PC4[idx_pca]
  merged_data$PC5 <- df.pc$PC5[idx_pca]
  merged_data$PC6 <- df.pc$PC6[idx_pca]
  merged_data$PC7 <- df.pc$PC7[idx_pca]
  merged_data$PC8 <- df.pc$PC8[idx_pca]
  merged_data$PC9 <- df.pc$PC9[idx_pca]
  merged_data$PC10 <- df.pc$PC10[idx_pca]
  merged_data$GDP <- df.pc$GDP[idx_pca]
  merged_data$key_pca <- NULL  # Clean up temporary key
  df.pc$key_pca <- NULL

  # Save the merged data with PCA covariates
  write.csv(merged_data, output_path, row.names = FALSE)
  message("[data_processing] wrote ", output_path, " (rows=", nrow(merged_data), ", ",
          sum(!is.na(merged_data$Antibiotic.Consumption)), " rows with consumption, ",
          sum(complete.cases(merged_data)), " complete cases)")

  merged_data_sums <- merged_data %>%
    {
      if (class_mode == "specific_abs") {
        group_by(., Pathogen, Antibiotic) %>%
          summarise(
            ATC.Class = paste(sort(unique(ATC.Class)), collapse = "|"),
            Total.Isolates = sum(Total.Isolates, na.rm = TRUE),
            Percent.Resistant.Isolates = mean(Percent.Resistant.Isolates, na.rm = TRUE),
            Antibiotic.Consumption = mean(Antibiotic.Consumption, na.rm = TRUE),
            PC1 = mean(PC1, na.rm = TRUE),
            PC2 = mean(PC2, na.rm = TRUE),
            PC3 = mean(PC3, na.rm = TRUE),
            PC4 = mean(PC4, na.rm = TRUE),
            PC5 = mean(PC5, na.rm = TRUE),
            PC6 = mean(PC6, na.rm = TRUE),
            PC7 = mean(PC7, na.rm = TRUE),
            PC8 = mean(PC8, na.rm = TRUE),
            PC9 = mean(PC9, na.rm = TRUE),
            PC10 = mean(PC10, na.rm = TRUE),
            GDP = mean(GDP, na.rm = TRUE),
            .groups = "drop"
          )
      } else {
        group_by(., Pathogen, ATC.Class) %>%
          summarise(
            Total.Isolates = sum(Total.Isolates, na.rm = TRUE),
            Percent.Resistant.Isolates = mean(Percent.Resistant.Isolates, na.rm = TRUE),
            Antibiotic.Consumption = mean(Antibiotic.Consumption, na.rm = TRUE),
            PC1 = mean(PC1, na.rm = TRUE),
            PC2 = mean(PC2, na.rm = TRUE),
            PC3 = mean(PC3, na.rm = TRUE),
            PC4 = mean(PC4, na.rm = TRUE),
            PC5 = mean(PC5, na.rm = TRUE),
            PC6 = mean(PC6, na.rm = TRUE),
            PC7 = mean(PC7, na.rm = TRUE),
            PC8 = mean(PC8, na.rm = TRUE),
            PC9 = mean(PC9, na.rm = TRUE),
            PC10 = mean(PC10, na.rm = TRUE),
            GDP = mean(GDP, na.rm = TRUE),
            .groups = "drop"
          )
      }
    } %>%
    ungroup()
  # save the merged data with sums
  write.csv(merged_data_sums, sums_output_path, row.names = FALSE)
  message("[data_processing] wrote ", sums_output_path, " (",
          nrow(merged_data_sums), " pathogen-class combinations)")
}

prepare_main_finer_regression_data <- function(
    finer_data_path = "finer_data_new.csv",
    consumption_path = "antibiotic_consumption_by_ATC3.csv",
    pca_path = "Chungman/Chungman_pca_renamed.csv",
    no_covariates_path = "finer_data_new_no_covariates.csv",
    output_path = "finer_data_new.csv",
    sums_output_path = "finer_data_sums_new.csv",
    year_cutoff = 2018
) {
  merged_data <- read.csv(finer_data_path)

  required_cols <- c("ISO3", "Year", "Pathogen", "ATC.Class", "Antibiotic", "Total.Isolates", "Percent.Resistant.Isolates")
  missing_cols <- setdiff(required_cols, names(merged_data))
  if (length(missing_cols) > 0) {
    stop(
      "[data_processing] finer data is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  merged_data <- merged_data %>%
    mutate(
      Year = as.character(Year),
      Total.Isolates = as.numeric(Total.Isolates),
      Percent.Resistant.Isolates = as.numeric(Percent.Resistant.Isolates)
    )

  # Harmonize finer aminoglycoside labels before any filtering or summarization.
  merged_data <- merged_data %>%
    mutate(
      Antibiotic = case_when(
        Antibiotic %in% c("Aminoglycosides", "Aminoglycosides (high-level)") ~ "Other aminoglycosides",
        TRUE ~ Antibiotic
      )
    )

  merged_data <- merged_data[!merged_data$ATC.Class %in% c("Total", "Other", "J01X"), ]
  merged_data <- merged_data[as.numeric(merged_data$Year) <= year_cutoff, ]
  merged_data <- merged_data[merged_data$ISO3 != "HKG", ]
  merged_data <- merged_data[!grepl("gonorrhoeae", tolower(merged_data$Pathogen)), ]

  consumption <- read.csv(consumption_path)
  consumption$ATC.level.3.class <- substr(consumption$ATC.level.3.class, 1, 4)
  consumption <- consumption %>%
    rename(
      ISO3 = Location,
      Year = Year,
      ATC.Class = ATC.level.3.class,
      Consumption = Antibiotic.consumption..DDD.1.000.day.
    ) %>%
    mutate(ISO3 = iso3_ihme_mapping$iso3[match(ISO3, iso3_ihme_mapping$country_name)])

  consumption <- consumption %>%
    group_by(ISO3, Year, ATC.Class) %>%
    summarise(Consumption = sum(Consumption, na.rm = TRUE), .groups = "drop")

  consumption$Year <- as.character(consumption$Year)

  merged_data <- merged_data %>%
    left_join(consumption, by = c("ISO3", "Year", "ATC.Class")) %>%
    rename(Antibiotic.Consumption = Consumption)

  wide_consumption <- consumption %>%
    pivot_wider(
      names_from = ATC.Class,
      values_from = Consumption,
      names_glue = "{ATC.Class}.Consumption"
    )

  merged_data <- merged_data %>%
    left_join(wide_consumption, by = c("ISO3", "Year")) %>%
    mutate(
      across(
        ends_with(".Consumption") & !matches("^Antibiotic\\.Consumption$"),
        ~ ifelse(paste0(ATC.Class, ".Consumption") == cur_column(), NA_real_, .x)
      )
    )

  write.csv(merged_data, no_covariates_path, row.names = FALSE)

  df.pc <- read.csv(pca_path)
  merged_data$key_pca <- paste(merged_data$ISO3, merged_data$Year, sep = "|")
  df.pc$key_pca <- paste(df.pc$ISO3, df.pc$Year, sep = "|")

  idx_pca <- match(merged_data$key_pca, df.pc$key_pca)
  merged_data$PC1 <- df.pc$PC1[idx_pca]
  merged_data$PC2 <- df.pc$PC2[idx_pca]
  merged_data$PC3 <- df.pc$PC3[idx_pca]
  merged_data$PC4 <- df.pc$PC4[idx_pca]
  merged_data$PC5 <- df.pc$PC5[idx_pca]
  merged_data$PC6 <- df.pc$PC6[idx_pca]
  merged_data$PC7 <- df.pc$PC7[idx_pca]
  merged_data$PC8 <- df.pc$PC8[idx_pca]
  merged_data$PC9 <- df.pc$PC9[idx_pca]
  merged_data$PC10 <- df.pc$PC10[idx_pca]
  merged_data$GDP <- df.pc$GDP[idx_pca]
  merged_data$key_pca <- NULL
  df.pc$key_pca <- NULL

  write.csv(merged_data, output_path, row.names = FALSE)
  message("[data_processing] wrote ", output_path, " (rows=", nrow(merged_data), ")")

  merged_data_sums <- merged_data %>%
    group_by(Pathogen, Antibiotic, ATC.Class) %>%
    summarise(
      Total.Isolates = sum(Total.Isolates, na.rm = TRUE),
      Percent.Resistant.Isolates = mean(Percent.Resistant.Isolates, na.rm = TRUE),
      Antibiotic.Consumption = mean(Antibiotic.Consumption, na.rm = TRUE),
      PC1 = mean(PC1, na.rm = TRUE),
      PC2 = mean(PC2, na.rm = TRUE),
      PC3 = mean(PC3, na.rm = TRUE),
      PC4 = mean(PC4, na.rm = TRUE),
      PC5 = mean(PC5, na.rm = TRUE),
      PC6 = mean(PC6, na.rm = TRUE),
      PC7 = mean(PC7, na.rm = TRUE),
      PC8 = mean(PC8, na.rm = TRUE),
      PC9 = mean(PC9, na.rm = TRUE),
      PC10 = mean(PC10, na.rm = TRUE),
      GDP = mean(GDP, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    ungroup()

  write.csv(merged_data_sums, sums_output_path, row.names = FALSE)
  message("[data_processing] wrote ", sums_output_path, " (", nrow(merged_data_sums), " pathogen-antibiotic combinations)")
}
