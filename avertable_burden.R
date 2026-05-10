## Extrapolating linear regression results to avertable burden
## BJS May 2025
library(tidyverse)
library(dplyr)
source("utils.R")

load_and_process_ihme_data <- function(
    ihme_fitted_path = "IHME_AMR/IHME_AMR_fitted_gammas_v2.csv",
    ihme_pathogen_path = "IHME_AMR/IHME_AMR_PATHOGEN_2019_DATA_COUNTED_AB.CSV",
    pop_path = "population_by_country_and_year.csv",
    consumption_path = "antibiotic_consumption_by_ATC3.csv",
    results_path =  "Outputs/database_gradients_pathogen_ATC3_PCA_canonical_weighted_main.csv",
    results_bootstrap_path = "Outputs/database_gradients_bootstraps_pathogen_ATC3_PCA_canonical_weighted_main.csv",
    gradients_path = "Outputs/database_gradients_ATC3_PCA_canonical_weighted_all.csv",
    gradients_bootstrap_path = "Outputs/database_gradients_bootstraps_ATC3_PCA_canonical_weighted_all.csv"
) {
  ## Load the results from the linear regression
  results <- read.csv(results_path)
  results_bootstrap <- read.csv(results_bootstrap_path)

  gradients_df <- read.csv(gradients_path)
  gradients_bootstrap <- read.csv(gradients_bootstrap_path)
  class_gradients <- as.vector(gradients_df[,2])
  classes <- as.vector(gradients_df[,1])
  names(class_gradients) <- classes

  ## IHME data reformatting
  IHME <- read.csv(ihme_fitted_path)
  lower_region_names <- unique(iso3_ihme_mapping$lower_ihme_region)

  # Disambiguate lower regions that share names with upper regions
  # (e.g. South Asia, North Africa and Middle East) using location_id.
  lower_region_id_map <- unique(IHME[, c("location_name", "location_id")]) %>%
      filter(location_name %in% lower_region_names) %>%
      group_by(location_name) %>%
      summarise(location_id = min(location_id), .groups = "drop")

  lower_region_ids <- lower_region_id_map$location_id

  # Upper-region (7 super-regions) disambiguation — used for Figure 3 scenarios.
  upper_region_names <- unique(iso3_ihme_mapping$ihme_region)
  upper_region_id_map <- unique(IHME[, c("location_name", "location_id")]) %>%
      filter(location_name %in% upper_region_names) %>%
      group_by(location_name) %>%
      summarise(location_id = min(location_id), .groups = "drop")
  upper_region_ids <- upper_region_id_map$location_id

  # Apply the mapping
  for(atc_code in names(atc_mapping)) {
    IHME[IHME$antibiotic_class %in% atc_mapping[[atc_code]], "antibiotic_class"] <- atc_code
  }
  # Save full dataset (all locations, post-atc-recode) for upper-region filtering.
  IHME_all <- IHME
  # keep only lower-IHME-region rows by ID (not name) to avoid collisions with
  # identically named upper regions.
  IHME <- IHME[IHME$location_id %in% lower_region_ids,]
  n_regions <- length(unique(IHME$location_name))
  # apply bacteria mapping
  for (i in 1:nrow(IHME)) {
      ihme_name <- IHME$pathogen[i]
      match <- bacteria_mapping[bacteria_mapping$in_names == ihme_name, "canonical_names"]
      if (length(match) > 0 && !is.na(match)) {
          IHME$pathogen[i] <- match
      }
  }
  # remove duplicate rows
  IHME <- IHME[!duplicated(IHME),]
  # remove duplicates based on "location_name","age_group_name","infectious_syndrome","pathogen","antibiotic_class"
  IHME <- IHME %>%
      group_by(location_name, age_group_name, infectious_syndrome, pathogen, antibiotic_class) %>%
      slice(1) %>%
      ungroup()

  # Build upper-region IHME (Figure 3 scenarios) — same atc recode already applied.
  IHME_upper <- IHME_all[IHME_all$location_id %in% upper_region_ids, ]
  for (i in seq_len(nrow(IHME_upper))) {
      ihme_name <- IHME_upper$pathogen[i]
      match_val <- bacteria_mapping[bacteria_mapping$in_names == ihme_name, "canonical_names"]
      if (length(match_val) > 0 && !is.na(match_val)) {
          IHME_upper$pathogen[i] <- match_val
      }
  }
  IHME_upper <- IHME_upper[!duplicated(IHME_upper), ]
  IHME_upper <- IHME_upper %>%
      group_by(location_name, age_group_name, infectious_syndrome, pathogen, antibiotic_class) %>%
      slice(1) %>%
      ungroup()

  IHME_totals_raw <- read.csv(ihme_pathogen_path)
  IHME_totals <- IHME_totals_raw
  pop_by_country_year <- read.csv(pop_path)
  # aggregate population by lower_ihme_region for 2018
  pop_by_lower_ihme_region <- data.frame(
      lower_ihme_region = unique(iso3_ihme_mapping$lower_ihme_region),
      population_2018 = numeric(length(unique(iso3_ihme_mapping$lower_ihme_region)))
  )
  lower_ihme_regions_for_pop <- unique(iso3_ihme_mapping$lower_ihme_region)
  for (i in seq_along(lower_ihme_regions_for_pop)) {
      region <- lower_ihme_regions_for_pop[i]
      countries_in_region <- unique(iso3_ihme_mapping$iso3[
          iso3_ihme_mapping$lower_ihme_region == region])
      total_population <- sum(pop_by_country_year[
          pop_by_country_year$Country.Code %in% countries_in_region, "X2018"])
      pop_by_lower_ihme_region$population_2018[i] <- total_population
  }
  write.csv(pop_by_lower_ihme_region,
      "Outputs/population_by_lower_ihme_region_2018.csv",
      row.names = FALSE)

  locations <- lower_region_names
  IHME_totals <- IHME_totals %>%
    filter(location_id %in% lower_region_ids) %>%
    select(-location_id) %>%
    distinct()
  total_burden <- sum(IHME_totals$val)
  print(total_burden)

  # total burden by lower-IHME-region
  total_burden_by_region <- data.frame(
      region = unique(IHME_totals$location_name),
      total_burden = numeric(length(unique(IHME_totals$location_name))),
      population = numeric(length(unique(IHME_totals$location_name)))
  )
  for (loc in unique(IHME_totals$location_name)){
      region_burden <- sum(IHME_totals[IHME_totals$location_name == loc, "val"])
      total_burden_by_region$total_burden[total_burden_by_region$region == loc] <- region_burden
      population <- pop_by_lower_ihme_region$population_2018[
          pop_by_lower_ihme_region$lower_ihme_region == loc]
      total_burden_by_region$population[total_burden_by_region$region == loc] <- population
      print(paste0("Burden in ",loc,": ",region_burden))
  }
  # save total burden by lower-IHME-region
  write.csv(total_burden_by_region, "Outputs/total_bacterial_disease_burden_by_lower_ihme_region_v2.csv", row.names = FALSE)

  # Upper-region burden totals (for Figure 3 scenario weighting)
  IHME_totals_upper <- IHME_totals_raw %>%
      filter(location_id %in% upper_region_ids) %>%
      select(-location_id) %>%
      distinct()
  # Population by upper IHME region (7 super-regions)
  pop_by_ihme_region <- data.frame(
      ihme_region = upper_region_names,
      population_2018 = numeric(length(upper_region_names))
  )
  for (i in seq_along(upper_region_names)) {
      region <- upper_region_names[i]
      countries_in_region <- unique(iso3_ihme_mapping$iso3[
          iso3_ihme_mapping$ihme_region == region])
      pop_by_ihme_region$population_2018[i] <- sum(pop_by_country_year[
          pop_by_country_year$Country.Code %in% countries_in_region, "X2018"])
  }
  total_burden_by_region_upper <- data.frame(
      region       = unique(IHME_totals_upper$location_name),
      total_burden = numeric(length(unique(IHME_totals_upper$location_name))),
      population   = numeric(length(unique(IHME_totals_upper$location_name)))
  )
  for (loc in unique(IHME_totals_upper$location_name)) {
      region_burden <- sum(IHME_totals_upper[IHME_totals_upper$location_name == loc, "val"])
      total_burden_by_region_upper$total_burden[
          total_burden_by_region_upper$region == loc] <- region_burden
      population <- pop_by_ihme_region$population_2018[
          pop_by_ihme_region$ihme_region == loc]
      total_burden_by_region_upper$population[
          total_burden_by_region_upper$region == loc] <- population
  }

  proportion_attributable_total <- sum(IHME$true_val_att)/(sum(IHME_totals$val))
  print(paste0("Total bacterial disease burden: ", total_burden))
  print(paste0("Proportion attributable to resistant infections: ", proportion_attributable_total))

  consumption <- read.csv(consumption_path)
  # filter consumption to location rows that are in iso3_ihme_mapping$ihme_region
  consumption <- consumption[consumption$Year == 2018,]
  consumption <- consumption %>%
      rename(
          Antibiotic = ATC.level.3.class,
          Consumption = Antibiotic.consumption..DDD.1.000.day.
      )
  # crop antibiotic name before hyphen
  consumption$Antibiotic <- sub("-.*", "", consumption$Antibiotic)
  global_consumption <- consumption[consumption$Location == "Global",]
  consumption <- consumption[consumption$Location %in% iso3_ihme_mapping$ihme_region,]

  return(list(
    IHME = IHME,
    IHME_upper = IHME_upper,
    IHME_totals = IHME_totals,
    results = results,
    results_bootstrap = results_bootstrap,
    class_gradients = class_gradients,
    classes = classes,
    gradients_bootstrap = gradients_bootstrap,
    pop_by_lower_ihme_region = pop_by_lower_ihme_region,
    pop_by_country_year = pop_by_country_year,
    total_burden_by_region = total_burden_by_region,
    total_burden_by_region_upper = total_burden_by_region_upper,
    consumption = consumption,
    global_consumption = global_consumption,
    proportion_attributable_total = proportion_attributable_total
  ))
}

## Load the results from the linear regression
data_loaded <- load_and_process_ihme_data()
IHME <- data_loaded$IHME
IHME_totals <- data_loaded$IHME_totals
results <- data_loaded$results
results_bootstrap <- data_loaded$results_bootstrap
class_gradients <- data_loaded$class_gradients
classes <- data_loaded$classes
gradients_bootstrap <- data_loaded$gradients_bootstrap
pop_by_lower_ihme_region <- data_loaded$pop_by_lower_ihme_region
pop_by_country_year <- data_loaded$pop_by_country_year
total_burden_by_region <- data_loaded$total_burden_by_region
IHME_upper <- data_loaded$IHME_upper
total_burden_by_region_upper <- data_loaded$total_burden_by_region_upper
consumption <- data_loaded$consumption
global_consumption <- data_loaded$global_consumption
proportion_attributable_total <- data_loaded$proportion_attributable_total


# Initialize with all combinations at full consumption
optimistic_df <- expand.grid(
  Pathogen = c(unique(results$Pathogen), "Overall"),
  Location = unique(consumption$Location),
  Antibiotic = unique(consumption$Antibiotic),
  stringsAsFactors = FALSE
) %>%
  as.data.frame()
optimistic_df$ProportionateConsumption <- 1

proportions <- c()
print("Optimistic scenario")
for (location in unique(consumption$Location)) {
  for (pathogen in unique(results$Pathogen)) {
    # print(paste("Processing location:", location, "and pathogen:", pathogen))
    quota <- 0.1 * sum(consumption$Consumption[consumption$Location == location])
    # print(paste("Initial quota:", quota))
    # sort drugs by Response for this pathogen
    pathogen_results <- results[results$Pathogen == pathogen,]
    pathogen_results <- pathogen_results[order(pathogen_results$Response, decreasing = TRUE),]
    for (antibiotic in pathogen_results$Antibiotic) {
      if (quota <= 0) {
        break
      }
      drug_consumption <- consumption$Consumption[consumption$Location == location & consumption$Antibiotic == antibiotic]
      if (length(drug_consumption) == 0) {
        next
      }
      new_consumption <- max(0, drug_consumption - quota)
      optimistic_df[optimistic_df$Pathogen == pathogen & optimistic_df$Location == location & optimistic_df$Antibiotic == antibiotic, "ProportionateConsumption"] <-
        new_consumption / drug_consumption
      quota <- quota - drug_consumption
      # print(paste("Processed antibiotic:", antibiotic, "Remaining quota:", quota))
    }
  }
  # then include a non-pathogen specific step where you just look at the class_gradients
  quota <- 0.1 * sum(consumption$Consumption[consumption$Location == location])
  # iterate from largest class_gradient to smallest
  for (antibiotic in names(class_gradients[order(class_gradients, decreasing = TRUE)])) {
    if (quota <= 0) {
      break
    }
    drug_consumption <- consumption$Consumption[consumption$Location == location & consumption$Antibiotic == antibiotic]
    if (length(drug_consumption) == 0) {
      next
    }
    new_consumption <- max(0, drug_consumption - quota)
    print(paste("Antibiotic:", antibiotic, "reduced to", (drug_consumption - new_consumption)/drug_consumption, "% in location:", location))
    optimistic_df[optimistic_df$Pathogen == "Overall" & optimistic_df$Location == location & optimistic_df$Antibiotic == antibiotic, "ProportionateConsumption"] <-
      new_consumption / drug_consumption
    quota <- quota - drug_consumption
  }
}
# write csv
write.csv(optimistic_df, "Outputs/optimistic_proportionate_consumption_by_pathogen_location_antibiotic.csv", row.names = FALSE)

# Initialize with all combinations at full consumption
pessimistic_df <- expand.grid(
  Pathogen = c(unique(results$Pathogen), "Overall"),
  Location = unique(consumption$Location),
  Antibiotic = unique(consumption$Antibiotic),
  stringsAsFactors = FALSE
) %>%
  as.data.frame()
pessimistic_df$ProportionateConsumption <- 1

proportions <- c()
print("Pessimistic scenario")
for (location in unique(consumption$Location)) {
  for (pathogen in unique(results$Pathogen)) {
    # print(paste("Processing location:", location, "and pathogen:", pathogen))
    quota <- 0.1 * sum(consumption$Consumption[consumption$Location == location])
    # print(paste("Initial quota:", quota))
    # sort drugs by Response for this pathogen
    pathogen_results <- results[results$Pathogen == pathogen,]
    pathogen_results <- pathogen_results[order(pathogen_results$Response, decreasing = FALSE),]
    for (antibiotic in pathogen_results$Antibiotic) {
      if (quota <= 0) {
        break
      }
      drug_consumption <- consumption$Consumption[consumption$Location == location & consumption$Antibiotic == antibiotic]
      if (length(drug_consumption) == 0) {
        next
      }
      new_consumption <- max(0, drug_consumption - quota)
      pessimistic_df[pessimistic_df$Pathogen == pathogen & pessimistic_df$Location == location & pessimistic_df$Antibiotic == antibiotic, "ProportionateConsumption"] <-
        new_consumption / drug_consumption
      quota <- quota - drug_consumption
      # print(paste("Processed antibiotic:", antibiotic, "Remaining quota:", quota))
    }
  }
  # then include a non-pathogen specific step where you just look at the class_gradients
  quota <- 0.1 * sum(consumption$Consumption[consumption$Location == location])
  # iterate from largest class_gradient to smallest
  for (antibiotic in names(class_gradients[order(class_gradients, decreasing = FALSE)])) {
    if (quota <= 0) {
      break
    }
    drug_consumption <- consumption$Consumption[consumption$Location == location & consumption$Antibiotic == antibiotic]
    if (length(drug_consumption) == 0) {
      next
    }
    new_consumption <- max(0, drug_consumption - quota)
    print(paste("Antibiotic:", antibiotic, "reduced by", (drug_consumption - new_consumption)/drug_consumption, "in location:", location))
    pessimistic_df[pessimistic_df$Pathogen == "Overall" & pessimistic_df$Location == location & pessimistic_df$Antibiotic == antibiotic, "ProportionateConsumption"] <-
      new_consumption / drug_consumption
    quota <- quota - drug_consumption
  }
}
# write csv
write.csv(pessimistic_df, "Outputs/pessimistic_proportionate_consumption_by_pathogen_location_antibiotic.csv", row.names = FALSE)

##### Calculate GDP and use by lower IHME region for 2018, to use in meta-regression

gdp_by_country_year <- read.csv("Chungman/Chungman_pca_renamed.csv")
gdp_by_lower_ihme_region <- data.frame(
    lower_ihme_region = unique(iso3_ihme_mapping$lower_ihme_region),
    gdp_2018 = numeric(length(unique(iso3_ihme_mapping$lower_ihme_region)))
)
use_by_country_year <- read.csv("DDD_country_year_class.csv")
# what's the relative use in each class in 2018?
use_by_class <- use_by_country_year[use_by_country_year$Year == 2018, ] %>%
    group_by(Antimicrobial) %>%
    summarise(total_DDD = sum(DDD))
print(use_by_class)

use_by_lower_ihme_region <- data.frame(
    lower_ihme_region = unique(iso3_ihme_mapping$lower_ihme_region),
    use_2018 = numeric(length(unique(iso3_ihme_mapping$lower_ihme_region)))
)
lower_ihme_regions <- unique(iso3_ihme_mapping$lower_ihme_region)
for (i in seq_along(lower_ihme_regions)) {
    region <- lower_ihme_regions[i]
    countries_in_region <- unique(iso3_ihme_mapping$iso3[
        iso3_ihme_mapping$lower_ihme_region == region])
    # make sure there's gdp data for all
    countries_in_region <- countries_in_region[
        countries_in_region %in% gdp_by_country_year[
            gdp_by_country_year$Year == "2018", "ISO3"]]
    # weight by population
    populations <- numeric(length(countries_in_region))
    gdp <- numeric(length(countries_in_region))
    use <- numeric(length(countries_in_region))
    # ensure same order as GDP data
    for (j in seq_along(countries_in_region)) {
        iso3 <- countries_in_region[j]
        populations[j] <- pop_by_country_year[
            pop_by_country_year$Country.Code == iso3, "X2018"
        ]
        gdp[j] <- gdp_by_country_year[
            (gdp_by_country_year$Year == "2018") &
                (gdp_by_country_year$ISO3 == iso3), "GDP"
        ]
        use[j] <- sum(use_by_country_year[
            (use_by_country_year$Year == 2018) &
                (use_by_country_year$ISO3 == iso3), "DDD"
        ])
    }
    total_population <- sum(populations)
    weighted_gdp <- sum(gdp * (populations / total_population),
        na.rm = TRUE)
    gdp_by_lower_ihme_region$gdp_2018[i] <- weighted_gdp
    if (region == "Western Sub-Saharan Africa") {
        total_use <- sum(use_by_country_year[
            (use_by_country_year$Year == 2018) &
                (use_by_country_year$ISO3 == "FWA"), "DDD"
        ])
        use_by_lower_ihme_region$use_2018[i] <- total_use / total_population /
            365 * 1000
        next
    }
    use_per_1000_person_days <- sum(use) / total_population / 365 * 1000
    use_by_lower_ihme_region$use_2018[i] <- use_per_1000_person_days
}
write.csv(gdp_by_lower_ihme_region,
    "Outputs/gdp_by_lower_ihme_region_2018_test.csv",
    row.names = FALSE)
write.csv(use_by_lower_ihme_region,
    "Outputs/use_by_lower_ihme_region_2018_test.csv",
    row.names = FALSE)

# =============================================================================
# Top-level call: run burden bootstrap computation and Figure 4.
# Deferred below the function definitions so they are available.
# These are executed after all function definitions have been parsed.
# =============================================================================
# (Executed at end of file via .amr_run_burden_toplevel)
.amr_run_burden_toplevel <- TRUE

# =============================================================================
# compute_avertable_burden(): Rescued from archive/burden_bootstrap_loop_archived.r
# Runs the per-row gamma bootstrap loop over IHME data, then aggregates by
# region, pathogen, drug, region x pathogen, and region x drug.
# Writes canonical output CSVs under Outputs/.
#
# Parameters:
#   IHME               : data frame with gamma params (shape_all, scale_all,
#                        a_frac, b_frac) and location_name, pathogen,
#                        antibiotic_class columns
#   results_bootstrap  : bootstrap gradient samples per pathogen x antibiotic
#                        (columns: Pathogen, Antibiotic, Gradient)
#   gradients_bootstrap: bootstrap class-level gradient samples
#                        (columns: Antibiotic, Consumption)
#   optimistic_df      : proportionate consumption under optimistic reduction
#   pessimistic_df     : proportionate consumption under pessimistic reduction
#   total_burden_by_region : data frame with region, total_burden, population
#   n_bootstraps       : number of bootstrap iterations (default 1000)
#   output_tag         : string appended to all canonical output file names
# =============================================================================
compute_avertable_burden <- function(
    IHME,
    results_bootstrap,
    gradients_bootstrap,
    optimistic_df,
    pessimistic_df,
    total_burden_by_region,
    n_bootstraps = 1000,
    output_tag = "canonical_weighted_lower_region_v2"
) {
    set.seed(260116)
    lower_to_upper <- unique(iso3_ihme_mapping[, c("lower_ihme_region", "ihme_region")])
    n_rows <- nrow(IHME)
    optimistic_burden_bootstraps  <- matrix(0, nrow = n_rows, ncol = n_bootstraps)
    pessimistic_burden_bootstraps <- matrix(0, nrow = n_rows, ncol = n_bootstraps)

    for (i in seq_len(n_rows)) {
        if (i %% 5000 == 0) message("[burden] Row ", i, " of ", n_rows)
        pathogen   <- IHME$pathogen[i]
        antibiotic <- IHME$antibiotic_class[i]
        location   <- IHME$location_name[i]

        if (is.na(pathogen) || is.na(antibiotic) || is.na(location)) next
        if (antibiotic == "Other") next

        shape_all  <- IHME$shape_all[i]
        scale_all  <- IHME$scale_all[i]
        a_frac     <- IHME$a_frac[i]
        b_frac     <- IHME$b_frac[i]
        value_gamma <- rgamma(n_bootstraps, shape = shape_all, scale = scale_all) * a_frac * b_frac

        # Sample gradients: pathogen-specific first, fall back to class-level
        path_boots <- results_bootstrap[
            results_bootstrap$Pathogen   == pathogen &
            results_bootstrap$Antibiotic == antibiotic, "Gradient"]
        if (length(path_boots) >= n_bootstraps) {
            gradients <- sample(path_boots, n_bootstraps, replace = TRUE)
        } else {
            class_boots <- gradients_bootstrap[
                gradients_bootstrap$Antibiotic == antibiotic, "Gradient"]
            if (length(class_boots) == 0) next
            gradients <- sample(class_boots, n_bootstraps, replace = TRUE)
        }

        # Optimistic burden
        opt_row <- optimistic_df[
            optimistic_df$Pathogen == "Overall" &
            optimistic_df$Location == location &
            optimistic_df$Antibiotic == antibiotic, "ProportionateConsumption"]
        if (length(opt_row) == 0 || is.na(opt_row)) {
          upper_loc <- lower_to_upper$ihme_region[
            match(location, lower_to_upper$lower_ihme_region)]
          if (length(upper_loc) > 0 && !is.na(upper_loc)) {
            opt_row <- optimistic_df[
              optimistic_df$Pathogen == "Overall" &
              optimistic_df$Location == upper_loc &
              optimistic_df$Antibiotic == antibiotic, "ProportionateConsumption"]
          }
        }
        if (length(opt_row) == 0 || is.na(opt_row)) {
            opt_row <- 0.9  # default: uniform 10% reduction
        }
        if (opt_row == 0) {
            optimistic_burden_bootstraps[i, ] <- value_gamma
        } else {
            optimistic_burden_bootstraps[i, ] <- value_gamma * (1 - exp(gradients * log(opt_row)))
        }

        # Pessimistic burden
        pes_row <- pessimistic_df[
            pessimistic_df$Pathogen == "Overall" &
            pessimistic_df$Location == location &
            pessimistic_df$Antibiotic == antibiotic, "ProportionateConsumption"]
        if (length(pes_row) == 0 || is.na(pes_row)) {
          upper_loc <- lower_to_upper$ihme_region[
            match(location, lower_to_upper$lower_ihme_region)]
          if (length(upper_loc) > 0 && !is.na(upper_loc)) {
            pes_row <- pessimistic_df[
              pessimistic_df$Pathogen == "Overall" &
              pessimistic_df$Location == upper_loc &
              pessimistic_df$Antibiotic == antibiotic, "ProportionateConsumption"]
          }
        }
        if (length(pes_row) == 0 || is.na(pes_row)) {
            pes_row <- 0.9
        }
        if (pes_row == 0) {
            pessimistic_burden_bootstraps[i, ] <- value_gamma
        } else {
            pessimistic_burden_bootstraps[i, ] <- value_gamma * (1 - exp(gradients * log(pes_row)))
        }
    }

    # Helper: aggregate bootstrap matrix by a grouping vector
    aggregate_bootstraps <- function(boots, groups, group_vals, pop_df) {
        n_groups <- length(group_vals)
        result <- data.frame(
            region              = group_vals,
            avertable_burden    = numeric(n_groups),
            lower_bound         = numeric(n_groups),
            upper_bound         = numeric(n_groups),
            variance            = numeric(n_groups),
            total_burden        = numeric(n_groups),
            proportion_avertable              = numeric(n_groups),
            proportion_avertable_lower_bound  = numeric(n_groups),
            proportion_avertable_upper_bound  = numeric(n_groups),
            proportion_avertable_variance     = numeric(n_groups),
            avertable_burden_per_100k         = numeric(n_groups),
            lower_bound_per_100k              = numeric(n_groups),
            upper_bound_per_100k              = numeric(n_groups),
            variance_per_100k                 = numeric(n_groups),
            stringsAsFactors = FALSE
        )
        for (i in seq_len(n_groups)) {
            g <- group_vals[i]
            row_idx <- which(groups == g)
            if (length(row_idx) == 0) next
            sums <- colSums(boots[row_idx, , drop = FALSE])
            pop_row <- pop_df[pop_df$region == g, ]
            pop <- if (nrow(pop_row) > 0) pop_row$population[1] else NA
            tb  <- if (nrow(pop_row) > 0) pop_row$total_burden[1] else NA
            result$avertable_burden[i]   <- mean(sums)
            result$lower_bound[i]        <- quantile(sums, 0.025)
            result$upper_bound[i]        <- quantile(sums, 0.975)
            result$variance[i]           <- var(sums)
            result$total_burden[i]       <- tb
            if (!is.na(pop) && pop > 0) {
                result$avertable_burden_per_100k[i] <- mean(sums)  / (pop / 1e5)
                result$lower_bound_per_100k[i]      <- quantile(sums, 0.025) / (pop / 1e5)
                result$upper_bound_per_100k[i]      <- quantile(sums, 0.975) / (pop / 1e5)
                result$variance_per_100k[i]         <- var(sums)   / ((pop / 1e5)^2)
            }
        }
        result$proportion_avertable              <- result$avertable_burden   / result$total_burden
        result$proportion_avertable_lower_bound  <- result$lower_bound        / result$total_burden
        result$proportion_avertable_upper_bound  <- result$upper_bound        / result$total_burden
        result$proportion_avertable_variance     <- result$variance           / result$total_burden^2
        result
    }

    regions   <- unique(IHME$location_name)
    pathogens <- unique(IHME$pathogen)
    drugs     <- unique(IHME$antibiotic_class)

    # Aggregate by region
    avertable_by_region <- aggregate_bootstraps(
        optimistic_burden_bootstraps, IHME$location_name, regions, total_burden_by_region)
    write.csv(avertable_by_region,
        paste0("Outputs/10pc_avertable_burden_by_region_", output_tag, ".csv"),
        row.names = FALSE)

    # Aggregate by pathogen (no population normalisation needed here)
    avertable_by_pathogen <- do.call(rbind, lapply(pathogens, function(p) {
        idx  <- which(IHME$pathogen == p)
        sums <- if (length(idx) > 0) colSums(optimistic_burden_bootstraps[idx, , drop = FALSE]) else rep(0, n_bootstraps)
        data.frame(pathogen = p,
                   avertable_burden = mean(sums),
                   lower_bound      = quantile(sums, 0.025),
                   upper_bound      = quantile(sums, 0.975),
                   stringsAsFactors = FALSE)
    }))
    write.csv(avertable_by_pathogen,
        paste0("Outputs/10pc_avertable_burden_by_pathogen_", output_tag, ".csv"),
        row.names = FALSE)

    # Aggregate by region x pathogen
    pairs_rp <- expand.grid(region = regions, pathogen = pathogens, stringsAsFactors = FALSE)
    avertable_by_pathogen_and_region <- do.call(rbind, lapply(seq_len(nrow(pairs_rp)), function(k) {
        r  <- pairs_rp$region[k]
        p  <- pairs_rp$pathogen[k]
        idx <- which(IHME$location_name == r & IHME$pathogen == p)
        sums <- if (length(idx) > 0) colSums(optimistic_burden_bootstraps[idx, , drop = FALSE]) else rep(0, n_bootstraps)
        pop_row <- total_burden_by_region[total_burden_by_region$region == r, ]
        pop <- if (nrow(pop_row) > 0) pop_row$population[1] else NA
        data.frame(region = r, pathogen = p,
                   avertable_burden = mean(sums),
                   lower_bound      = quantile(sums, 0.025),
                   upper_bound      = quantile(sums, 0.975),
                   avertable_burden_per_100k = if (!is.na(pop) && pop > 0) mean(sums) / (pop / 1e5) else NA,
                   lower_bound_per_100k      = if (!is.na(pop) && pop > 0) quantile(sums, 0.025) / (pop / 1e5) else NA,
                   upper_bound_per_100k      = if (!is.na(pop) && pop > 0) quantile(sums, 0.975) / (pop / 1e5) else NA,
                   stringsAsFactors = FALSE)
    }))
    write.csv(avertable_by_pathogen_and_region,
        paste0("Outputs/10pc_avertable_burden_by_pathogen_and_region_", output_tag, ".csv"),
        row.names = FALSE)

    # Aggregate by region x drug
    pairs_rd <- expand.grid(region = regions, drug = drugs, stringsAsFactors = FALSE)
    avertable_by_drug_and_region <- do.call(rbind, lapply(seq_len(nrow(pairs_rd)), function(k) {
        r   <- pairs_rd$region[k]
        d   <- pairs_rd$drug[k]
        idx <- which(IHME$location_name == r & IHME$antibiotic_class == d)
        sums <- if (length(idx) > 0) colSums(optimistic_burden_bootstraps[idx, , drop = FALSE]) else rep(0, n_bootstraps)
        pop_row <- total_burden_by_region[total_burden_by_region$region == r, ]
        pop <- if (nrow(pop_row) > 0) pop_row$population[1] else NA
        data.frame(region = r, drug = d,
                   avertable_burden = mean(sums),
                   lower_bound      = quantile(sums, 0.025),
                   upper_bound      = quantile(sums, 0.975),
                   avertable_burden_per_100k = if (!is.na(pop) && pop > 0) mean(sums) / (pop / 1e5) else NA,
                   lower_bound_per_100k      = if (!is.na(pop) && pop > 0) quantile(sums, 0.025) / (pop / 1e5) else NA,
                   upper_bound_per_100k      = if (!is.na(pop) && pop > 0) quantile(sums, 0.975) / (pop / 1e5) else NA,
                   stringsAsFactors = FALSE)
    }))
    write.csv(avertable_by_drug_and_region,
        paste0("Outputs/10pc_avertable_burden_by_drug_and_region_", output_tag, ".csv"),
        row.names = FALSE)

    message("[burden] compute_avertable_burden() complete. Output tag: ", output_tag)
    invisible(list(
        by_region             = avertable_by_region,
        by_pathogen           = avertable_by_pathogen,
        by_pathogen_and_region = avertable_by_pathogen_and_region,
        by_drug_and_region    = avertable_by_drug_and_region
    ))
}

# =============================================================================
# Deferred top-level execution: compute_avertable_burden()
# Called after all function definitions above have been parsed.
# Figure 4 is generated by the figures stage (generate_figure4() in plotting.R).
# =============================================================================
if (isTRUE(.amr_run_burden_toplevel)) {
    .burden_scenario  <- getOption("amr_scenario", "main")
    .n_boot <- if (isTRUE(getOption("amr_smoke_mode", FALSE))) 10L else 1000L

  # Main scenario: uniform 10% reduction across antibiotic classes.
  main_uniform_df <- expand.grid(
    Pathogen = "Overall",
    Location = unique(IHME$location_name),
    Antibiotic = unique(IHME$antibiotic_class),
    stringsAsFactors = FALSE
  )
  main_uniform_df$ProportionateConsumption <- 0.9

    if (.burden_scenario %in% c("main", "burden_lower_region")) {
        compute_avertable_burden(
            IHME                   = IHME,
            results_bootstrap      = results_bootstrap,
            gradients_bootstrap    = gradients_bootstrap,
      optimistic_df          = main_uniform_df,
            pessimistic_df         = pessimistic_df,
            total_burden_by_region = total_burden_by_region,
            n_bootstraps           = .n_boot,
            output_tag             = "canonical_weighted_lower_region_v2"
        )
        # Optimistic scenario (reduce highest-gradient drugs first) — upper-region IHME
        compute_avertable_burden(
            IHME                   = IHME_upper,
            results_bootstrap      = results_bootstrap,
            gradients_bootstrap    = gradients_bootstrap,
            optimistic_df          = optimistic_df,
            pessimistic_df         = pessimistic_df,
            total_burden_by_region = total_burden_by_region_upper,
            n_bootstraps           = .n_boot,
            output_tag             = "canonical_weighted_upper_region_optimistic_overall"
        )
        # Pessimistic scenario (reduce lowest-gradient drugs first) — upper-region IHME
        compute_avertable_burden(
            IHME                   = IHME_upper,
            results_bootstrap      = results_bootstrap,
            gradients_bootstrap    = gradients_bootstrap,
            optimistic_df          = pessimistic_df,
            pessimistic_df         = pessimistic_df,
            total_burden_by_region = total_burden_by_region_upper,
            n_bootstraps           = .n_boot,
            output_tag             = "canonical_weighted_upper_region_pessimistic_overall"
        )
    }

    if (.burden_scenario %in% c("burden_upper_region", "burden_optimistic",
                                "burden_drug_region", "burden_pathogen_region")) {
        compute_avertable_burden(
            IHME                   = IHME_upper,
            results_bootstrap      = results_bootstrap,
            gradients_bootstrap    = gradients_bootstrap,
            optimistic_df          = optimistic_df,
            pessimistic_df         = pessimistic_df,
            total_burden_by_region = total_burden_by_region_upper,
            n_bootstraps           = .n_boot,
            output_tag             = "canonical_weighted_upper_region_optimistic_overall"
        )
    }

    if (.burden_scenario == "burden_pessimistic") {
        compute_avertable_burden(
            IHME                   = IHME_upper,
            results_bootstrap      = results_bootstrap,
            gradients_bootstrap    = gradients_bootstrap,
            optimistic_df          = pessimistic_df,
            pessimistic_df         = pessimistic_df,
            total_burden_by_region = total_burden_by_region_upper,
            n_bootstraps           = .n_boot,
            output_tag             = "canonical_weighted_upper_region_pessimistic_overall"
        )
    }

    message("[burden] Burden estimation complete. Intermediate CSVs written to Outputs/.")
    message("[burden] Figure 4 is generated by 'make figures' (generate_figure4() in plotting.R).")
}
