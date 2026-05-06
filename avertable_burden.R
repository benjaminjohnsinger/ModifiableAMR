## Extrapolating linear regression results to avertable burden
## BJS May 2025
library(tidyverse)
library(dplyr)
source("utils.R")

## Load the results from the linear regression
results <- read.csv("Outputs/database_gradients_pathogen_ATC3_PCA_joelike_weighted.csv")
results_bootstrap <- read.csv("Outputs/database_gradients_bootstraps_pathogen_ATC3_PCA_joelike_weighted.csv")

gradients_df <- read.csv("Outputs/database_gradients_ATC3_PCA_joelike_weighted.csv")
gradients_bootstrap <- read.csv("Outputs/database_gradients_bootstraps_ATC3_PCA_joelike_weighted.csv")
class_gradients <- as.vector(gradients_df[,2])
classes <- as.vector(gradients_df[,1])
names(class_gradients) <- classes

## IHME data reformatting
IHME <- read.csv("IHME_AMR/IHME_AMR_fitted_gammas_v2.csv")
# Apply the mapping
for(atc_code in names(atc_mapping)) {
  IHME[IHME$antibiotic_class %in% atc_mapping[[atc_code]], "antibiotic_class"] <- atc_code
}
# only keep rows where location is in iso3_ihme_mapping$ihme_region
IHME <- IHME[IHME$location_name %in% iso3_ihme_mapping$ihme_region,]
n_regions <- length(unique(IHME$location_name))
# apply bacteria mapping
for (i in 1:nrow(IHME)) {
    ihme_name <- IHME$pathogen[i]
    match <- bacteria_mapping[bacteria_mapping$in_names == ihme_name, "joe_names"]
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
# IHME_totals <- IHME %>%
#     group_by(location_name, age_group_name, infectious_syndrome, pathogen) %>%
#     slice(1) %>%
#     ungroup()
# total_burden <- sum(IHME_totals$true_val_all)

IHME_totals <- read.csv("IHME_AMR/IHME_AMR_PATHOGEN_2019_DATA_COUNTED_AB.CSV")
pop_by_country_year <- read.csv("population_by_country_and_year.csv")
# aggregate population by ihme_region for 2018. years are columns
pop_by_ihme_region <- data.frame(
    ihme_region = unique(iso3_ihme_mapping$ihme_region),
    population_2018 = numeric(length(unique(iso3_ihme_mapping$ihme_region)))
)
ihme_regions <- unique(iso3_ihme_mapping$ihme_region)
for (i in seq_along(ihme_regions)) {
    region <- ihme_regions[i]
    countries_in_region <- unique(iso3_ihme_mapping$iso3[
        iso3_ihme_mapping$ihme_region == region])
    total_population <- sum(pop_by_country_year[
        pop_by_country_year$Country.Code %in% countries_in_region, "X2018"])
    pop_by_ihme_region$population_2018[i] <- total_population
}
write.csv(pop_by_ihme_region,
    "Outputs/population_by_ihme_region_2018.csv",
    row.names = FALSE)
# IHME_totals <- IHME_totals[IHME_totals$location_name %in% iso3_ihme_mapping$ihme_region,]
# IHME_totals <- IHME_totals[!duplicated(IHME_totals),]
locations <- unique(iso3_ihme_mapping$ihme_region)
IHME_totals <- IHME_totals %>%
  filter(location_name %in% locations) %>%
  select(-location_id) %>%
  distinct()
total_burden <- sum(IHME_totals$val)
print(total_burden)
# total burden by location
total_burden_by_region <- data.frame(
    region = unique(IHME_totals$location_name),
    total_burden = numeric(length(unique(IHME_totals$location_name))),
    population = numeric(length(unique(IHME_totals$location_name)))
)
for (loc in unique(IHME_totals$location_name)){
    region_burden <- sum(IHME_totals[IHME_totals$location_name == loc, "val"])
    total_burden_by_region$total_burden[total_burden_by_region$region == loc] <- region_burden
    population <- pop_by_ihme_region$population_2018[
        pop_by_ihme_region$ihme_region == loc]
    total_burden_by_region$population[total_burden_by_region$region == loc] <- population
    print(paste0("Burden in ",loc,": ",region_burden))
}
# save total burden by region
write.csv(total_burden_by_region, "Outputs/total_bacterial_disease_burden_by_ihme_region_v2.csv", row.names = FALSE)
proportion_attributable_total <- sum(IHME$true_val_att)/(sum(IHME_totals$val))
print(paste0("Total bacterial disease burden: ", total_burden))
print(paste0("Proportion attributable to resistant infections: ", proportion_attributable_total))

consumption <- read.csv("antibiotic_consumption_by_ATC3.csv")
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


# for each row in IHME, take true_val_att (the total attributable burden), and calculate the avertable burden with the formula burden averted = burden x (1 - gradient x 0.1 x use), using the gradient from the results table for the pathogen and antibiotic, and the use from the consumption table for the location and antibiotic
# create new column in IHME called avertable_burden
n_bootstraps <- 1000
burden_bootstraps <- matrix(0, nrow = nrow(IHME), ncol = n_bootstraps)
optimistic_burden_bootstraps <- matrix(0, nrow = nrow(IHME), ncol = n_bootstraps)
pessimistic_burden_bootstraps <- matrix(0, nrow = nrow(IHME), ncol = n_bootstraps)
total_burden_bootstraps <- matrix(0, nrow = nrow(IHME), ncol = n_bootstraps)
# gamma_means <- numeric(nrow(IHME))
# gamma_vals <- numeric(nrow(IHME))
# gradient_means <- numeric(nrow(IHME))
# gradient_vals <- numeric(nrow(IHME))

# # random seed for reproducibility
# set.seed(260116)

# for (i in 1:nrow(IHME)) {
#     if (i %% 10000 == 0) {
#         print(paste("Processing row", i, "of", nrow(IHME)))
#     }
#     pathogen <- IHME$pathogen[i]
#     antibiotic <- IHME$antibiotic_class[i]
#     location <- IHME$location_name[i]
#     # print(paste(pathogen, antibiotic, location))
#     if (is.na(pathogen) || is.na(antibiotic) || is.na(location)) {
#         next
#     }
#     # # check that location is in consumption
#     # if (!location %in% consumption$Location) {
#     #     next
#     # }
#     if (antibiotic == "Other"){
#         next
#     }
#     # gamma_means[i] <- mean(value_gamma)
#     # gamma_vals[i] <- IHME$true_val_att[i]
#     if (!any(results$Pathogen == pathogen & results$Antibiotic == antibiotic)) {
#         # check if antibiotic has overall result
#         if (antibiotic %in% classes){
#             gradients <- sample(gradients_bootstrap[gradients_bootstrap$Antibiotic == antibiotic, "Consumption"], n_bootstraps, replace = TRUE)
#             pathogen <- "Overall"
#             # gradient_means[i] <- mean(gradients)
#             # gradient_vals[i] <- class_gradients[antibiotic]
#         } else {
#             next
#         }
#     } else {
#         gradients <- sample(results_bootstrap[results_bootstrap$Pathogen == pathogen & results_bootstrap$Antibiotic == antibiotic, "Gradient.Consumption"], n_bootstraps, replace = TRUE)
#         # gradient_means[i] <- mean(gradients)
#         # gradient_vals[i] <- results$Gradient[results$Pathogen == pathogen & results$Antibiotic == antibiotic]
#     }
#     shape_all <- IHME$shape_all[i]
#     scale_all <- IHME$scale_all[i]
#     a_frac <- IHME$a_frac[i]
#     b_frac <- IHME$b_frac[i]
#     value_gamma_all <- rgamma(n_bootstraps, shape = shape_all, scale = scale_all)
#     value_gamma <- value_gamma_all * a_frac * b_frac
#     # use <- consumption$Consumption[consumption$Location == location & consumption$Antibiotic == antibiotic]
#     # global_use <- global_consumption$Consumption[global_consumption$Antibiotic == antibiotic]
#     avertable_burden <- value_gamma * (1-exp(gradients * log(0.9)))
#     burden_bootstraps[i,] <- avertable_burden
#     pathogen <- "Overall"
#     optimistic_change <- optimistic_df$ProportionateConsumption[optimistic_df$Pathogen == pathogen & optimistic_df$Location == location & optimistic_df$Antibiotic == antibiotic]
#     if (optimistic_change == 0){
#       optimistic_avertable_burden <- value_gamma
#     } else {
#       optimistic_avertable_burden <- value_gamma * (1-exp(gradients * log(optimistic_change)))
#     }
#     optimistic_burden_bootstraps[i,] <- optimistic_avertable_burden
#     pessimistic_change <- pessimistic_df$ProportionateConsumption[pessimistic_df$Pathogen == pathogen & pessimistic_df$Location == location & pessimistic_df$Antibiotic == antibiotic]
#     if (pessimistic_change == 0){
#       pessimistic_avertable_burden <- value_gamma
#     } else {
#       pessimistic_avertable_burden <- value_gamma * (1-exp(gradients * log(pessimistic_change)))
#     }
#     pessimistic_burden_bootstraps[i,] <- pessimistic_avertable_burden

#     # for (j in 1:n_bootstraps) {
#     #     use <- consumption$Consumption[consumption$Location == location & consumption$Antibiotic == antibiotic]
#     #     avertable_burden <- attributable_value * max(0, gradient) * 0.01 * use
#     #     # print(gradient)
#     #     burden_bootstraps[i,j] <- avertable_burden
#     # }
# }
# # # scatter plots of means and vals
# # colors = c("#648FFF","#DC267F","#FFB000","#785EF0","#FF832B","#000000","#648FFF","#DC267F","#FFB000","#785EF0","#FF832B", "#000000","#648FFF","#DC267F","#FFB000")
# # shapes = c(18,18,18,18,18,18,16,16,16,16,16,16,15,15,15)

# # # Create a mapping of pathogens to colors and shapes
# # pathogen_list = c("Acinetobacter spp.", "E. faecalis", "E. faecium", "E. coli", "H. influenzae", "K. pneumoniae", "N. gonorrhoeae", "P. aeruginosa", "Proteus spp.", "Salmonella spp.", "Shigella spp.", "S. aureus", "S. agalactiae", "S. pneumoniae", "S. pyogenes")
# # pathogen_list = sort(pathogen_list)
# # color_mapping <- setNames(colors[1:length(pathogen_list)], pathogen_list)
# # shape_mapping <- setNames(shapes[1:length(pathogen_list)], pathogen_list)

# # ggplot(data.frame(gamma_means, gamma_vals, pathogen_class = paste(IHME$pathogen, IHME$antibiotic_class)), aes(x = gamma_means, y = gamma_vals, color = pathogen_class)) +
# #     geom_point() +
# #     labs(x = "Mean Value", y = "True Value") +
# #     ggtitle("Gamma Means vs True Values by Pathogen-ATC3 Class") +
# #     theme_minimal() +
# #     scale_x_log10() +
# #     scale_y_log10() +
# #     geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # Add x=y line
# #     theme(legend.position = "none")  # Remove legend
# # ggsave("gamma_mean_vs_true_values_by_pathogen_class.png", width = 8, height = 6)

# # ggplot(data.frame(gradient_means, gradient_vals, drug_class = IHME$antibiotic_class), aes(x = gradient_means, y = gradient_vals, color = drug_class)) +
# #     geom_point() +
# #     labs(x = "Mean Gradient Value", y = "Gradient Value") +
# #     ggtitle("Gradient Means vs Gradient Values") +
# #     theme_minimal() +
# #     scale_x_log10() +
# #     scale_y_log10() +
# #     geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")  # Add x=y line
# # ggsave("gradient_means_vs_gradient_values.png", width = 8, height = 6)
# # burden_bootstraps <- as.data.frame(burden_bootstraps)
# colnames(optimistic_burden_bootstraps) <- paste0("bootstrap_", 1:n_bootstraps)
# write.csv(optimistic_burden_bootstraps, "Outputs/10pc_avertable_burden_bootstraps_joelike_weighted_upper_region_optimistic_overall.csv", row.names = FALSE)
# colnames(pessimistic_burden_bootstraps) <- paste0("bootstrap_", 1:n_bootstraps)
# write.csv(pessimistic_burden_bootstraps, "Outputs/10pc_avertable_burden_bootstraps_joelike_weighted_upper_region_pessimistic_overall.csv", row.names = FALSE)

# # total avertable burden and 95% CI
# total_avertable_burden <- colSums(optimistic_burden_bootstraps)
# total_avertable_burden_mean <- mean(total_avertable_burden)
# total_avertable_burden_lower <- quantile(total_avertable_burden, 0.025)
# total_avertable_burden_upper <- quantile(total_avertable_burden, 0.975)
# # print total avertable burden and 95% CI
# print(paste("Total avertable burden (mean):", total_avertable_burden_mean))
# print(paste("Total avertable burden (95% CI):", total_avertable_burden_lower, "-", total_avertable_burden_upper))

# print(paste0("Proportion of total burden avertable: ", total_avertable_burden_mean / total_burden))
# print(paste0("Proportion of attributable burden avertable: ", (total_avertable_burden_mean / total_burden)/proportion_attributable_total))

# # sum bootstrap columns by region, drug, or pathogen
# n_regions <- length(unique(IHME$location_name))
# optimistic_burden_bootstraps_by_region <- matrix(0, nrow = n_regions, ncol = n_bootstraps)
# for (i in 1:n_regions) {
#     region <- unique(IHME$location_name)[i]
#     optimistic_burden_bootstraps_by_region[i,] <- colSums(optimistic_burden_bootstraps[IHME$location_name == region,])
# }
# n_drugs <- length(unique(IHME$antibiotic_class))
# optimistic_burden_bootstraps_by_drug <- matrix(0, nrow = n_drugs, ncol = n_bootstraps)
# for (i in 1:n_drugs) {
#     drug <- unique(IHME$antibiotic_class)[i]
#     optimistic_burden_bootstraps_by_drug[i,] <- colSums(optimistic_burden_bootstraps[IHME$antibiotic_class == drug,])
# }
# n_pathogens <- length(unique(IHME$pathogen))
# optimistic_burden_bootstraps_by_pathogen <- matrix(0, nrow = n_pathogens, ncol = n_bootstraps)
# for (i in 1:n_pathogens) {
#     pathogen <- unique(IHME$pathogen)[i]
#     optimistic_burden_bootstraps_by_pathogen[i,] <- colSums(optimistic_burden_bootstraps[IHME$pathogen == pathogen,])
# }
# rownames(optimistic_burden_bootstraps_by_pathogen) <- unique(IHME$pathogen)

# # for each pathogen get the mean and 95% CI of the avertable burden
# avertable_by_pathogen <- data.frame(
#     pathogen = unique(IHME$pathogen),
#     avertable_burden = numeric(length(unique(IHME$pathogen))),
#     lower_bound = numeric(length(unique(IHME$pathogen))),
#     upper_bound = numeric(length(unique(IHME$pathogen)))
# )
# for (i in 1:n_pathogens) {
#     avertable_by_pathogen$pathogen[i] <- unique(IHME$pathogen)[i]
#     avertable_by_pathogen$avertable_burden[i] <- mean(optimistic_burden_bootstraps_by_pathogen[i,])
#     avertable_by_pathogen$lower_bound[i] <- quantile(optimistic_burden_bootstraps_by_pathogen[i,], 0.025)
#     avertable_by_pathogen$upper_bound[i] <- quantile(optimistic_burden_bootstraps_by_pathogen[i,], 0.975)
# }
# write.csv(avertable_by_pathogen, "Outputs/10pc_avertable_burden_by_pathogen_joelike_weighted_upper_region_optimistic_overall.csv", row.names = FALSE)

# # for each region get the mean and 95% CI of the avertable burden
# avertable_by_region <- data.frame(
#     region = unique(IHME$location_name),
#     avertable_burden = numeric(length(unique(IHME$location_name))),
#     lower_bound = numeric(length(unique(IHME$location_name))),
#     upper_bound = numeric(length(unique(IHME$location_name))),
#     variance = numeric(length(unique(IHME$location_name))),
#     total_burden = numeric(length(unique(IHME$location_name))),
#     proportion_avertable = numeric(length(unique(IHME$location_name))),
#     proportion_avertable_lower_bound = numeric(length(unique(IHME$location_name))),
#     proportion_avertable_upper_bound = numeric(length(unique(IHME$location_name))),
#     proportion_avertable_variance = numeric(length(unique(IHME$location_name))),
#     avertable_burden_per_100k = numeric(length(unique(IHME$location_name))),
#     lower_bound_per_100k = numeric(length(unique(IHME$location_name))),
#     upper_bound_per_100k = numeric(length(unique(IHME$location_name))),
#     variance_per_100k = numeric(length(unique(IHME$location_name)))
# )
# for (i in 1:n_regions) {
#     region <- unique(IHME$location_name)[i]
#     avertable_by_region$region[i] <- region
#     avertable_by_region$avertable_burden[i] <- mean(optimistic_burden_bootstraps_by_region[i, ])
#     avertable_by_region$lower_bound[i] <- quantile(optimistic_burden_bootstraps_by_region[i, ], 0.025)
#     avertable_by_region$upper_bound[i] <- quantile(optimistic_burden_bootstraps_by_region[i, ], 0.975)
#     avertable_by_region$variance[i] <- var(optimistic_burden_bootstraps_by_region[i, ])
#     avertable_by_region$total_burden[i] <- total_burden_by_region[total_burden_by_region$region == region, "total_burden"]
#     avertable_by_region$avertable_burden_per_100k[i] <- mean(optimistic_burden_bootstraps_by_region[i, ]) / (total_burden_by_region[total_burden_by_region$region == region, "population"] / 100000)
#     avertable_by_region$lower_bound_per_100k[i] <- quantile(optimistic_burden_bootstraps_by_region[i, ], 0.025) / (total_burden_by_region[total_burden_by_region$region == region, "population"] / 100000)
#     avertable_by_region$upper_bound_per_100k[i] <- quantile(optimistic_burden_bootstraps_by_region[i, ], 0.975) / (total_burden_by_region[total_burden_by_region$region == region, "population"] / 100000)
#     avertable_by_region$variance_per_100k[i] <- var(optimistic_burden_bootstraps_by_region[i, ]) / ((total_burden_by_region[total_burden_by_region$region == region, "population"] / 100000)^2)
# }
# avertable_by_region$proportion_avertable <- avertable_by_region$avertable_burden / avertable_by_region$total_burden
# avertable_by_region$proportion_avertable_lower_bound <- avertable_by_region$lower_bound / avertable_by_region$total_burden
# avertable_by_region$proportion_avertable_upper_bound <- avertable_by_region$upper_bound / avertable_by_region$total_burden
# avertable_by_region$proportion_avertable_variance <- avertable_by_region$variance / (avertable_by_region$total_burden^2)
# write.csv(avertable_by_region, "Outputs/10pc_avertable_burden_by_region_joelike_weighted_upper_region_optimistic_overall.csv", row.names = FALSE)

# # aggregate optimistic_burden_bootstraps by pathogen and region, getting mean of and
# # 95% CI of the avertable burden
# optimistic_burden_bootstraps_by_pathogen_and_region <- matrix(0,
#   nrow = n_regions * n_pathogens, ncol = n_bootstraps)
# unique_pathogens <- unique(IHME$pathogen)
# unique_regions <- unique(IHME$location_name)

# row_index <- 1

# for (pathogen_name in unique_pathogens) {
#   for (region_name in unique_regions) {
#     summed_burden <- colSums(optimistic_burden_bootstraps[
#       IHME$pathogen == pathogen_name &
#         IHME$location_name == region_name, , drop = FALSE])
#     optimistic_burden_bootstraps_by_pathogen_and_region[row_index, ] <- summed_burden
#     row_index <- row_index + 1
#   }
# }
# # for each region and pathogen get the mean and 95% CI of the avertable burden
# avertable_by_pathogen_and_region <- data.frame(
#   pathogen = rep(unique_pathogens, each = n_regions),
#   region = rep(unique_regions, times = n_pathogens),
#   avertable_burden = numeric(n_regions * n_pathogens),
#   lower_bound = numeric(n_regions * n_pathogens),
#   upper_bound = numeric(n_regions * n_pathogens),
#   avertable_burden_per_100k = numeric(n_regions * n_pathogens),
#   lower_bound_per_100k = numeric(n_regions * n_pathogens),
#   upper_bound_per_100k = numeric(n_regions * n_pathogens)
# )
# # calculate the mean and 95% CI of the avertable burden for each region and
# # pathogen
# for (i in seq_len(nrow(avertable_by_pathogen_and_region))) {
#   avertable_by_pathogen_and_region$avertable_burden[i] <-
#     mean(optimistic_burden_bootstraps_by_pathogen_and_region[i, ])
#   avertable_by_pathogen_and_region$lower_bound[i] <-
#     quantile(optimistic_burden_bootstraps_by_pathogen_and_region[i, ], 0.025)
#   avertable_by_pathogen_and_region$upper_bound[i] <-
#     quantile(optimistic_burden_bootstraps_by_pathogen_and_region[i, ], 0.975)
  
#   # Get population for this region
#   region_name <- avertable_by_pathogen_and_region$region[i]
#   population <- total_burden_by_region[
#     total_burden_by_region$region == region_name, "population"]
  
#   # Calculate per 100k values
#   avertable_by_pathogen_and_region$avertable_burden_per_100k[i] <-
#     avertable_by_pathogen_and_region$avertable_burden[i] / 
#     (population / 100000)
#   avertable_by_pathogen_and_region$lower_bound_per_100k[i] <-
#     avertable_by_pathogen_and_region$lower_bound[i] / 
#     (population / 100000)
#   avertable_by_pathogen_and_region$upper_bound_per_100k[i] <-
#     avertable_by_pathogen_and_region$upper_bound[i] / 
#     (population / 100000)
# }
# write.csv(avertable_by_pathogen_and_region,
#   "Outputs/10pc_avertable_burden_by_pathogen_and_region_joelike_weighted_upper_region_optimistic_overall.csv",
#   row.names = FALSE)

# # aggregate optimistic_burden_bootstraps by drug and region, getting mean of and
# # 95% CI of the avertable burden
# optimistic_burden_bootstraps_by_drug_and_region <- matrix(0,
#   nrow = n_regions * n_drugs, ncol = n_bootstraps)
# unique_drugs <- unique(IHME$antibiotic_class)
# unique_regions <- unique(IHME$location_name)

# row_index <- 1

# for (drug_name in unique_drugs) {
#   for (region_name in unique_regions) {
#     summed_burden <- colSums(optimistic_burden_bootstraps[
#       IHME$antibiotic_class == drug_name &
#         IHME$location_name == region_name, , drop = FALSE])
#     optimistic_burden_bootstraps_by_drug_and_region[row_index, ] <- summed_burden
#     row_index <- row_index + 1
#   }
# }
# # for each region and pathogen get the mean and 95% CI of the avertable burden
# # print(unique_drugs)
# avertable_by_drug_and_region <- data.frame(
#   drug = rep(unique_drugs, each = n_regions),
#   region = rep(unique_regions, times = n_drugs),
#   avertable_burden = numeric(n_regions * n_drugs),
#   lower_bound = numeric(n_regions * n_drugs),
#   upper_bound = numeric(n_regions * n_drugs),
#   avertable_burden_per_100k = numeric(n_regions * n_drugs),
#   lower_bound_per_100k = numeric(n_regions * n_drugs),
#   upper_bound_per_100k = numeric(n_regions * n_drugs)
# )
# # print(avertable_by_drug_and_region)
# # calculate the mean and 95% CI of the avertable burden for each region and
# # drug
# for (i in seq_len(nrow(avertable_by_drug_and_region))) {
#   avertable_by_drug_and_region$avertable_burden[i] <-
#     mean(optimistic_burden_bootstraps_by_drug_and_region[i, ])
#   avertable_by_drug_and_region$lower_bound[i] <-
#     quantile(optimistic_burden_bootstraps_by_drug_and_region[i, ], 0.025)
#   avertable_by_drug_and_region$upper_bound[i] <-
#     quantile(optimistic_burden_bootstraps_by_drug_and_region[i, ], 0.975)
  
#   # Get population for this region
#   region_name <- avertable_by_drug_and_region$region[i]
#   population <- total_burden_by_region[
#     total_burden_by_region$region == region_name, "population"]
  
#   # Calculate per 100k values
#   avertable_by_drug_and_region$avertable_burden_per_100k[i] <-
#     avertable_by_drug_and_region$avertable_burden[i] / 
#     (population / 100000)
#   avertable_by_drug_and_region$lower_bound_per_100k[i] <-
#     avertable_by_drug_and_region$lower_bound[i] / 
#     (population / 100000)
#   avertable_by_drug_and_region$upper_bound_per_100k[i] <-
#     avertable_by_drug_and_region$upper_bound[i] / 
#     (population / 100000)
# }
# write.csv(avertable_by_drug_and_region,
#   "Outputs/10pc_avertable_burden_by_drug_and_region_joelike_weighted_upper_region_optimistic_overall.csv",
#   row.names = FALSE)

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

# Load metafor package for meta-regression
library(metafor)

file_path <- paste0("Outputs/10pc_avertable_burden_by_region_joelike_",
                    "weighted_lower_region_v2.csv")
avertable_by_region <- read.csv(file_path)

gdp_by_lower_ihme_region <- read.csv(
  "Outputs/gdp_by_lower_ihme_region_2018_test.csv")
use_by_lower_ihme_region <- read.csv(
  "Outputs/use_by_lower_ihme_region_2018_test.csv")

# plot proportion avertable vs gdp per capita
avertable_by_region <- merge(avertable_by_region,
                             gdp_by_lower_ihme_region,
                             by.x = "region",
                             by.y = "lower_ihme_region")

# Order regions by proportion avertable (descending)
avertable_by_region <- avertable_by_region[
                       order(avertable_by_region$proportion_avertable, 
                             decreasing = TRUE), ]

# Merge with antibiotic use data
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

# Predictions for totals vs GDP
pred_totals_gdp <- predict(meta_totals_gdp,
                newmods = cbind(gdp_range, gdp_range^2))
pred_df_totals_gdp <- data.frame(
  gdp = gdp_range,
  pred = pred_totals_gdp$pred,
  ci_lower = pred_totals_gdp$ci.lb,
  ci_upper = pred_totals_gdp$ci.ub
)

# Predictions for totals vs use
pred_totals_use <- predict(meta_totals_use,
                newmods = cbind(use_range, use_range^2))
pred_df_totals_use <- data.frame(
  use = use_range,
  pred = pred_totals_use$pred,
  ci_lower = pred_totals_use$ci.lb,
  ci_upper = pred_totals_use$ci.ub
)

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

# where do the predicted proportion curves peak?
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

# print population-weighted average of proportion avertable for regions South Asia, North Africa and Middle East, Eastern Europe, and Southern Latin America
total_burden_by_region <- read.csv("Outputs/total_bacterial_disease_burden_by_lower_ihme_region_v2.csv")
print(total_burden_by_region)
selected_regions <- c("Central Europe",
            "Eastern Europe", "Southern Latin America")
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
# print the gdp per capita for each of these regions
for (i in seq_along(selected_regions)) {
  region <- selected_regions[i]
  gdp <- avertable_by_region[
    avertable_by_region$region == region, "gdp_2018"]
  print(paste0("GDP per capita for ", region, ": ", round(gdp, 2)))
}
# Figure 1: Total avertable mortality
totals_by_region <- ggplot(avertable_by_region,
              aes(x = avertable_burden, 
                y = reorder(region, avertable_burden))) +
  geom_bar(stat = "identity", fill = "grey50") +
  geom_errorbar(aes(xmin = lower_bound, xmax = upper_bound), width = 0.2) +
  geom_text(aes(label = paste0(format(round(avertable_burden, 2), nsmall = 2), 
                 " (", format(round(lower_bound, 2), nsmall = 2), "–", 
                 format(round(upper_bound, 2), nsmall = 2), ")")),
        hjust = 0, size = 2.5, family = "Helvetica",
        x = max(avertable_by_region$upper_bound) + 
          0.02 * (max(avertable_by_region$upper_bound) - 
              min(avertable_by_region$lower_bound))) +
  labs(x = "Avertible AMR mortality (deaths)", y = "Region") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.title.x = element_text(size = 10, family = "Helvetica"),
      axis.title.y = element_text(size = 10, family = "Helvetica"),
      axis.text.x = element_text(size = 8, family = "Helvetica"),
      axis.text.y = element_text(size = 8, family = "Helvetica"),
      plot.title = element_text(size = 12, family = "Helvetica"),
      axis.ticks = element_line(color = "black"),
      text = element_text(family = "Helvetica")) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5.5, 100, 5.5, 5.5, "points"))

totals_vs_gdp <- ggplot(avertable_by_region,
             aes(x = gdp_2018, y = avertable_burden)) +
  geom_ribbon(data = pred_df_totals_gdp, aes(x = gdp, ymin = ci_lower,
                         ymax = ci_upper), 
        inherit.aes = FALSE, alpha = 0.2, fill = "black") +
  geom_line(data = pred_df_totals_gdp, aes(x = gdp, y = pred),
        inherit.aes = FALSE, color = "black", linewidth = 0.7, alpha = 0.7) +
  geom_point(size = 2, color = "black") +
  geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), width = 0.2,
          color = "black") +
  labs(x = "GDP per capita (USD PPP)", y = "Avertible mortality") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.title.x = element_text(size = 10, family = "Helvetica"),
      axis.title.y = element_text(size = 10, family = "Helvetica"),
      axis.text.x = element_text(size = 8, family = "Helvetica"),
      axis.text.y = element_text(size = 8, family = "Helvetica"),
      plot.title = element_text(size = 12, family = "Helvetica"),
      axis.ticks = element_line(color = "black"),
      text = element_text(family = "Helvetica"))

totals_vs_use <- ggplot(avertable_by_region_with_use,
             aes(x = use_2018, y = avertable_burden)) +
  geom_ribbon(data = pred_df_totals_use, aes(x = use, ymin = ci_lower,
                         ymax = ci_upper), 
        inherit.aes = FALSE, alpha = 0.2, fill = "black") +
  geom_line(data = pred_df_totals_use, aes(x = use, y = pred),
        inherit.aes = FALSE, color = "black", linewidth = 0.7, alpha = 0.7) +
  geom_point(size = 2, color = "black") +
  geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), width = 0.2,
          color = "black") +
  labs(x = "Antibiotic use (DDD/1000 person-days)",
     y = "Avertible mortality") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.title.x = element_text(size = 10, family = "Helvetica"),
      axis.title.y = element_text(size = 10, family = "Helvetica"),
      axis.text.x = element_text(size = 8, family = "Helvetica"),
      axis.text.y = element_text(size = 8, family = "Helvetica"),
      plot.title = element_text(size = 12, family = "Helvetica"),
      axis.ticks = element_line(color = "black"),
      text = element_text(family = "Helvetica"))

# Figure 2: Proportion avertible (converted to percentages)
proportions_by_region <- ggplot(avertable_by_region,
                 aes(x = proportion_avertable * 100,
                   y = reorder(region, 
                         proportion_avertable))) +
  geom_bar(stat = "identity", fill = "grey50") +
  geom_errorbar(aes(xmin = proportion_avertable_lower_bound * 100,
            xmax = proportion_avertable_upper_bound * 100), 
          width = 0.2) +
  geom_text(aes(label = paste0(format(round(proportion_avertable * 100, 2), 
                    nsmall = 2), 
                 "% (", 
                 format(round(proportion_avertable_lower_bound * 
                        100, 2), nsmall = 2), 
                 "-", 
                 format(round(proportion_avertable_upper_bound * 
                        100, 2), nsmall = 2), 
                 "%)")),
        # hjust = 0, size = 2.5, family = "Helvetica",
        hjust = 0, size = 5, family = "Helvetica",
        x = max(avertable_by_region$proportion_avertable_upper_bound * 
            100) + 
          0.1 * (max(avertable_by_region$proportion_avertable_upper_bound * 
               100) - 
             min(avertable_by_region$proportion_avertable_lower_bound * 
               100))) +
  labs(x = "Percentage of bacterial mortality avertible (%)", y = "Region", 
     tag = "A") +
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
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5.5, 150, 5.5, 5.5, "points"))

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

# Figure for per 100k population vs GDP
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

# Figure for per 100k population vs use
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

# Create Figure 1: Total avertable mortality with top panel spanning width
figure1 <- gridExtra::grid.arrange(
  totals_by_region,
  gridExtra::arrangeGrob(totals_vs_gdp, totals_vs_use, ncol = 2),
  heights = c(1, 1)
)

# Create Panel A figure (proportions by region)
panel_A <- proportions_by_region +
  labs(tag = "") +  # Remove tag since it will be the only panel
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(size = 24))

# Create Panels B & C figure (proportions vs GDP and use)
panels_B_C <- gridExtra::arrangeGrob(
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

panels_D_E <- gridExtra::arrangeGrob(
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

# Create Panels B-E figure
panels_B_E <- gridExtra::arrangeGrob(
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

# Save Panel A as PowerPoint slide
ggsave("Figure4_Panel_A_narrow.pdf", panel_A,
       width = 9, height = 7.5, units = "in")

# Save Panels B & C as PowerPoint slide
ggsave("Figure4_Panels_B_C_narrow.pdf", panels_B_C,
       width = 9, height = 7.5, units = "in")

# Save Panels D & E as PowerPoint slide
ggsave("Figure4_Panels_D_E_narrow.pdf", panels_D_E,
       width = 9, height = 7.5, units = "in")

# Save Panels B-E as PowerPoint slide
ggsave("Figure4_Panels_B_E_narrow.pdf", panels_B_E,
       width = 9, height = 7.5, units = "in")
    

# Save the plots
ggsave("Figure4.pdf", figure2,
     width = 6.5, height = 7, units = "in")
ggsave("percentage_avertable_mortality_figure.png", figure2,
     width = 6.5, height = 7, units = "in")

# get the avertable burden per 100,000 population by drug and region
avertable_by_drug_and_region <- merge(avertable_by_drug_and_region, pop_by_lower_ihme_region, by.x = "region", by.y = "lower_ihme_region")
avertable_by_drug_and_region$avertable_burden_per_100k <- (avertable_by_drug_and_region$avertable_burden / avertable_by_drug_and_region$population_2018) * 100000
avertable_by_drug_and_region$lower_bound_per_100k <- (avertable_by_drug_and_region$lower_bound / avertable_by_drug_and_region$population_2018) * 100000
avertable_by_drug_and_region$upper_bound_per_100k <- (avertable_by_drug_and_region$upper_bound / avertable_by_drug_and_region$population_2018) * 100000
write.csv(avertable_by_drug_and_region, "Outputs/10pc_avertable_burden_per_100k_by_drug_and_region_joelike_weighted_w_GASP.csv", row.names = FALSE)

# # remove J01X and any other antibiotics that are not in the ATC3 list
# avertable_by_antibiotic_class <- avertable_by_antibiotic_class[!avertable_by_antibiotic_class$antibiotic_class %in% c("J01X", "Other", "Multi-drug resistance in Salmonella Typhi and Paratyphi"),]

# # sum avertable burden for rows with same "pathogen","antibiotic_class","location_name"
# avertable <- IHME[c("location_name","pathogen","antibiotic_class","true_val_att","avertable_burden")] %>%
#     group_by(location_name, pathogen, antibiotic_class) %>%
#     summarise(attributable_burden = sum(true_val_att, na.rm = TRUE), avertable_burden = sum(avertable_burden, na.rm = TRUE), .groups = 'drop')
# avertable_by_country <- IHME[c("location_name","pathogen","antibiotic_class","true_val_att","avertable_burden")] %>%
#     group_by(location_name) %>%
#     summarise(attributable_burden = sum(true_val_att, na.rm = TRUE), avertable_burden = sum(avertable_burden, na.rm = TRUE), .groups = 'drop')
# avertable_by_drug <- IHME[c("pathogen","antibiotic_class","true_val_att","avertable_burden")] %>%
#     group_by(antibiotic_class) %>%
#     summarise(attributable_burden = sum(true_val_att, na.rm = TRUE), avertable_burden = sum(avertable_burden, na.rm = TRUE), .groups = 'drop')
# avertable_by_pathogen <- IHME[c("pathogen","antibiotic_class","true_val_att","avertable_burden")] %>%
#     group_by(pathogen) %>%
#     summarise(attributable_burden = sum(true_val_att, na.rm = TRUE), avertable_burden = sum(avertable_burden, na.rm = TRUE), .groups = 'drop')
# horizontal bar chart of deaths avertable by pathogen
plot <- ggplot(avertable_by_pathogen, aes(x = avertable_burden, y = fct_reorder(pathogen, avertable_burden))) +
    geom_bar(stat = "identity", fill = "grey50") +
    geom_errorbar(aes(xmin = lower_bound, xmax = upper_bound), width = 0.2, position = position_dodge(0.9)) +
    labs(x = "Deaths averted", y = "") +
    ggtitle("Deaths avertable by a 10% reduction in antibiotic use, by pathogen") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid = element_blank()) +
    theme(axis.title.x = element_text(size = 20)) +
    theme(axis.title.y = element_text(size = 20)) +
    theme(axis.text.x = element_text(size = 14)) +
    theme(axis.text.y = element_text(size = 14)) +
    theme(plot.title = element_text(size = 24))  # Set title size
# save the plot
ggsave("10pc_avertable_attributable_burden_by_pathogen_joelike_slide.png", plot, width = 13.3, height = 7.5, units = "in")

# # how many rows in IHME have NA avertable_burden
# print(paste("Proportion of rows with NA avertable_burden:", na_rows / nrow(IHME)))
# save IHME as csv
# write.csv(IHME, "IHME_AMR/IHME_AMR_with_10pc_avertable_burden.csv", row.names = FALSE)
# write.csv(avertable, "IHME_AMR/10pc_avertable_burden.csv", row.names = FALSE)
# write.csv(avertable_by_country, "IHME_AMR/10pc_avertable_burden_by_country.csv", row.names = FALSE)
# write.csv(avertable_by_drug, "IHME_AMR/10pc_avertable_burden_by_drug.csv", row.names = FALSE)
# write.csv(avertable_by_pathogen, "IHME_AMR/10pc_avertable_burden_by_pathogen.csv", row.names = FALSE)
