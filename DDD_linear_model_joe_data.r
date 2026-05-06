# ## Linear regression on database DDD data
# ## BJS March 2025
library(tidyverse)
library(dplyr)
source("utils.R")

# # # ### Nagorsen data version
# data <- read.csv("Nagorsen_clean.csv", colClasses = c("units" = "character"), na.strings = c("NA"))

# # # exclude rows where "amt_consumed" is NA
# data <- data[!is.na(data$amt_consumed),]
# # exclude rows where "units" is NA
# data <- data[!is.na(data$units),]
# # exclude amt_consumed above 10000
# data <- data[data$amt_consumed < 10000,]

# data <- data[!is.na(data$class_for_resistance),]
# # map pathogens with get_bacteria_name
# for (i in seq_len(nrow(data))) {
#     data$pathogen[i] <- get_bacteria_name(data$pathogen[i])
# }

# # # Create a mapping dictionary for antibiotic classes to ATC codes

# # # Apply the mapping
# for(atc_code in names(atc_mapping)) {
#   data[data$class_for_resistance %in% atc_mapping[[atc_code]], "class_for_resistance"] <- atc_code
# }
# # print(length(unique(data$class_for_resistance)))
# # exclude J01X and Other
# data <- data[!data$class_for_resistance %in% c("J01X", "Other"),]

# df.pc <- read.csv("Chungman/Chungman_pca_renamed.csv")
# data$ISO3 <- iso3_ihme_mapping$iso3[match(data$country, iso3_ihme_mapping$country_name)]
# # chungman data only present in some of the years represented in Mikkel data
# for (i in seq_len(nrow(data))) {
#     year <- data$end_year[i]
#     iso3 <- data$ISO3[i]
#     if (iso3 %in% df.pc$ISO3) {
#         row_index <- which(df.pc$ISO3 == iso3 & df.pc$Year == year)
#         if (length(row_index) > 0) {
#             data$PC1[i] <- df.pc$PC1[row_index[1]]
#             data$PC2[i] <- df.pc$PC2[row_index[1]]
#             data$PC3[i] <- df.pc$PC3[row_index[1]]
#             data$GDP[i] <- df.pc$GDP[row_index[1]]
#         }
#     }
# }
# # This fold-change calculation should not be used - wait until after filtering
# # # # # fold-change calculation
# # # data <- data %>% 
# # #     group_by(study, class_for_resistance, pathogen) %>%
# # #     mutate(relative_use = amt_consumed / mean(amt_consumed, na.rm = TRUE), 
# # #            relative_resistance = percent_isolates_resistant / mean(percent_isolates_resistant, na.rm = TRUE)) %>%
# # #     ungroup()

# # Replacing units (based on reading original studies to confirm equivalence)
# # Where "units" is "DDD/100 bed days", divide "amt_consumed" by 10 and replace unit with "DDD/1000 bed days"
# data %>% filter(units == "DDD/100 bed days") %>% mutate(amt_consumed = amt_consumed/10, units = "DDD/1000 bed days")
# # Where "units" is "DDD/1000 women/year" divide "amt_consumed" by 365 and replace units with "DDD/1000 women/day"
# data %>% filter(units == "DDD/1000 women/year") %>% mutate(amt_consumed = amt_consumed/365, units = "DDD/1000 women/day")
# # Where "units" is "DDD/inhabitants/day", replace with "DDD/1000 inhabitants/day"
# data[data$units == "DDD/inhabitants/day", "units"] <- "DDD/1000 inhabitants/day"
# # Where "units" is "DDD/1000 inhabitants", replace with "DDD/1000 inhabitants/day"
# data[data$units == "DDD/1000 inhabitants", "units"] <- "DDD/1000 inhabitants/day"

# # Keep only rows in which "units" contains "DDD", "1000", and "day"
# data <- data[grepl("DDD", data$units) & grepl("1000", data$units) & grepl("day", data$units),]

# # # print(unique(data$class_for_resistance))
# # print(unique(data$pathogen))
# # # # Keep only rows in which neither ab_setting or pathogen_setting contain the strings "hospital" or "Hospital"
# # print(table(data$ab_setting))
# # print(table(data$pathogen_setting))

# data <- data[!grepl("community", data$ab_setting),]
# # data <- data[!((grepl("community", data$pathogen_setting) | grepl("Community", data$pathogen_setting))|grepl("communty", data$pathogen_setting)),]

# # # print(unique(data$class_for_resistance))
# # # print(unique(data$pathogen))
# data <- data %>%
#     select(Consumption = amt_consumed,
#            Resistance = percent_isolates_resistant,
#            Pathogen = pathogen,
#            DOI = doi,
#            Antibiotic = class_for_resistance,
#            PC1 = PC1,
#            PC2 = PC2,
#            PC3 = PC3,
#            GDP = GDP,
#            Year = end_year,
#            Weight = 1,
#     )

# # remove drug-bug combos with 20 or fewer entries
# pathogen_drug_counts <- table(paste(data$Pathogen, data$Antibiotic))
# pathogen_drug_to_remove <- names(pathogen_drug_counts[pathogen_drug_counts <= 20])
# data <- data[!paste(data$Pathogen, data$Antibiotic) %in% pathogen_drug_to_remove, ]
# data <- na.omit(data)

# print(unique(data$DOI))

# stop("Stopping before combining with Joe data")

# # # ### Data from Joe version
# # data <- read.csv("pathogen_abx_analysis_all_variables_(class-specific).csv")
# # # data <- read.csv("merged_data_test.csv")

# # # relabel "Antibiotic_standardized" to "Antibiotic"
# # data <- data %>%
# #     rename(
# #         Antibiotic = ATC.Class,
# #         Consumption = Antibiotic.consumption..DDD.1.000.day.,
# #         Resistance = Percent.Resistant.Isolates,
# #         Pathogen = Pathogen
# #         )
# # data[data$Antibiotic == "Total", "Antibiotic"] <- "Other"

# ## Combining data
# # remove columns other than Consumption, Resistance, Pathogen, and Antibiotic
# # data <- data %>%
# #     select(Consumption, Resistance, Pathogen, Antibiotic)
# # data2 <- data2 %>%
# #     select(Consumption, Resistance, Pathogen, Antibiotic)

# # data <- rbind(data, data2)

# # add column for World Bank lending group from iso3_ihme_mapping
# # data <- data %>%
# #     mutate(World.Bank.Lending.Group = case_when(
# #         National.Income.per.capita..log. < 7.04 ~ "Low income",
# #         National.Income.per.capita..log. >= 7.04 & National.Income.per.capita..log. < 8.42 ~ "Lower middle income",
# #         National.Income.per.capita..log. >= 8.42 & National.Income.per.capita..log. < 9.54 ~ "Upper middle income",
# #         National.Income.per.capita..log. >= 9.54 ~ "High income"
# #     ))

# # # filter to only High income countries
# # data <- data %>%
# #     filter(World.Bank.Lending.Group == "High income")

# # ## add column for GBD super-region - High Income, Central and Eastern Europe and Central Asia, North Africa and Middle East, Southeast and East Asia and Oceania, Latin America and Carribean, South Asia, and Sub-Saharan Africa
# # data <- data %>%
# #     mutate(GBD.Super.Region = case_when(
# #         ISO3 %in% c("BGD", "BTN", "IND", "NPL", "PAK") ~ "South Asia",
# #         ISO3 %in% c("ALB", "ARM", "AZE", "BIH", "BLR", "BGR", "HRV", "CZE", "EST", "GEO", "HUN", "KAZ", "KGZ", "LTU", "LVA", "MDA", "MKD", "MNE", "MNG", "POL", "ROU", "RUS", "SRB", "SVK", "SVN", "TJK", "TKM", "UZB", "UKR") ~ "Central and Eastern Europe and Central Asia",
# #         ISO3 %in% c("AFG", "BHR", "DZA", "EGY", "IRN", "IRQ", "JOR", "KWT", "LBN", "LBY", "MAR", "PSE", "OMN", "QAT", "SAU", "SDN", "SYR", "TUN", "TUR", "ARE", "YEM") ~ "North Africa and Middle East",
# #         ISO3 %in% c("CHN", "FJI", "HKG", "IDN", "LAO", "MYS", "MNG", "PNG", "PHL", "SGP", "THA", "TLS", "VNM") ~ "Southeast and East Asia and Oceania",
# #         ISO3 %in% c("BOL", "BRA", "COL", "CRI", "CUB", "DOM", "ECU", "GUY", "HND", "MEX", "NIC", "PRY", "SLV", "URY", "VEN") ~ "Latin America and Carribean",
# #         ISO3 %in% c("BWA", "CAF", "CIV", "CMR", "COD", "DJI", "DZA", "EGY", "ETH", "GHA", "GIN", "KEN", "LSO", "MDG", "MWI", "MYS", "NAM", "NER", "NGA", "RWA", "SEN", "SLE", "SOM", "SSD", "TCD", "TGO", "TZA", "UGA", "ZAF") ~ "Sub-Saharan Africa",
# #         ISO3 %in% c("AUS", "AND", "JPN", "KOR", "NZL", "ISR", "USA", "CHL", "ARG", "CAN", "GBR", "IRL", "FRA", "DEU", "NLD", "BEL", "CHE", "AUT", "DNK", "NOR", "SWE", "FIN", "ISL") ~ "High Income")
# #     )

# # data <- data %>%
# #     filter(GBD.Super.Region == "Central and Eastern Europe and Central Asia")



# # country_covariates = c(
# # # "National.Income.per.capita..log."
# # # ,"Global.Value.National.Income.per.capita..log.",
# # "National.Income.per.capita..log..Relative.to.Global.Value"
# # # ,"Mean.Years.Schooling"
# # # ,"Global.Value.Years.Schooling"
# # ,"Mean.Years.Schooling.Relative.to.Global.Value"
# # ,"Mortality.rate.attributed.to.unsafe.water..unsafe.sanitation.and.lack.of.hygiene.from.diarrhoea..intestinal.nematode.infections..malnutrition.and.acute.respiratory.infections..deaths.per.100.000.population..Relative.to.Global.Value"
# # ,"Proportion.of.population.practicing.open.defecation..across.all.locations.....Relative.to.Global.Value"
# # ,"Proportion.of.population.using.basic.drinking.water.services..across.all.locations.....Relative.to.Global.Value"
# # ,"Proportion.of.population.using.basic.sanitation.services..across.all.locations.....Relative.to.Global.Value"
# # ,"Proportion.of.population.using.safely.managed.drinking.water.services..across.all.locations.....Relative.to.Global.Value"
# # ,"Proportion.of.population.using.safely.managed.sanitation.services..across.all.locations.....Relative.to.Global.Value"
# # ,"Proportion.of.population.with.basic.handwashing.facilities.on.premises..across.all.locations.....Relative.to.Global.Value"
# # ,"Proportion.of.schools.with.access.to.basic.drinking.water..primary.schools.only.....Relative.to.Global.Value"
# # ,"Proportion.of.schools.with.access.to.single.sex.basic.sanitation..primary.schools.only.....Relative.to.Global.Value"
# # ,"Proportion.of.schools.with.basic.handwashing.facilities..primary.schools.only.....Relative.to.Global.Value"
# # ,"Proportion.of.wastewater.treated..across.all.locations.and.activities.....Relative.to.Global.Value"
# # ,"Population"
# # # ,"Life.Expectancy"
# # # ,"Global.Value.Life.Expectancy"
# # ,"Life.Expectancy.Relative.to.Global.Value")


# # # ## Merged data from Joe and ATLAS
data <- read.csv("merged_data_N_PC3_GDP.csv")
print(table(data$ISO3))
print(length(table(data$ISO3)))
print(max(table(data$ISO3)))
print(table(data$Year))
data <- data %>%
    rename(
        Antibiotic = ATC.Class,
        Consumption = Antibiotic.Consumption,
        Resistance = Percent.Resistant.Isolates,
        Pathogen = Pathogen,
        Weight = Total.Isolates
        )
country_covariates = c("PC1", "PC2", "PC3", "GDP", "Year")
# scale GDP
data$GDP <- data$GDP / mean(data$GDP, na.rm = TRUE)
# remove rows with NA in Consumption, Resistance, Pathogen, Antibiotic, or Weight
data <- data[!is.na(data$Consumption) & !is.na(data$Resistance) & !is.na(data$Pathogen) & !is.na(data$Antibiotic) & !is.na(data$Weight),]

# remove drug-pathogen combinations with 10 or fewer entries
pathogen_drug_counts <- table(paste(data$Pathogen, data$Antibiotic))
pathogen_drug_to_remove <- names(pathogen_drug_counts[pathogen_drug_counts <= 10])
data <- data[!paste(data$Pathogen, data$Antibiotic) %in% pathogen_drug_to_remove, ]


merged_data_sums <- read.csv("merged_data_sums_N.csv")
merged_data_sums <- merged_data_sums %>%
    rename(
        Antibiotic = ATC.Class,
    )
# remove Antibiotic-Pathogen combinations with fewer than 100 Total.Isolates in merged_data_sums
pathogen_drug_counts <- merged_data_sums %>%
    group_by(Pathogen, Antibiotic) %>%
    summarise(Total.Isolates = sum(Total.Isolates, na.rm = TRUE)) %>%
    filter(Total.Isolates < 100) %>%
    select(Pathogen, Antibiotic)
pathogen_drug_to_remove <- paste(pathogen_drug_counts$Pathogen, pathogen_drug_counts$Antibiotic)
data <- data[!paste(data$Pathogen, data$Antibiotic) %in% pathogen_drug_to_remove, ]

# remove drug-pathogen combinations with all zero resistance
pathogen_drug_counts <- data %>%
    group_by(Pathogen, Antibiotic) %>%
    summarise(Resistance = sum(Resistance, na.rm = TRUE)) %>%
    filter(Resistance == 0) %>%
    select(Pathogen, Antibiotic)
pathogen_drug_to_remove <- paste(pathogen_drug_counts$Pathogen, pathogen_drug_counts$Antibiotic)
data <- data[!paste(data$Pathogen, data$Antibiotic) %in% pathogen_drug_to_remove, ]
print(table(data$Antibiotic))
print(table(data$Pathogen))
print(sum(table(data$Pathogen)))

# list of countries included for each pathogen
countries_per_pathogen <- data %>%
    group_by(Pathogen) %>%
    summarise(Countries = paste(unique(ISO3), collapse = ", ")) %>%
    ungroup()
for (i in seq_len(nrow(countries_per_pathogen))) {
    print(paste0(countries_per_pathogen$Pathogen[i], ": ", countries_per_pathogen$Countries[i]))
}


data$lending_group <- iso3_ihme_mapping$lending_group[match(data$ISO3, iso3_ihme_mapping$iso3)]
# find the total number of entries in high income countries
high_income_entries <- sum(data$lending_group == "High income")
# find the total number of isolates in high income countries
high_income_isolates <- sum(data$Weight[data$lending_group == "High income"], na.rm = TRUE)
# print the number of entries and isolates in high income countries
print(paste("HIC entries:", high_income_entries))
print(paste("HIC isolates:", high_income_isolates))
# find the total number of entris in lmics
lmics_entries <- sum(data$lending_group != "High income")
# find the total number of isolates in lmics
lmics_isolates <- sum(data$Weight[data$lending_group != "High income"], na.rm = TRUE)
# print the number of entries and isolates in lmics
print(paste("LMIC entries:", lmics_entries))
print(paste("LMIC isolates:", lmics_isolates))


# # select lending group
data_LMIC <- data[data$lending_group != "High income",]
data_HIC <- data[data$lending_group == "High income",]
# drop lending group
data <- data %>%
    select(Consumption, Resistance, Pathogen, Antibiotic, Weight, ISO3, all_of(country_covariates))

# # # print(mean(data$Resistance, na.rm = TRUE))
# # # print(mean(data$Consumption, na.rm = TRUE))

# bar chart of number of observations per antibiotic and pathogen
# library(ggplot2)
# library(ggpattern)
# # instead of default colours, use 5 IBM colors + black, then repeat with cross-hatch, then finally cross-hatch the other way
# colors <- c("#648FFF", "#DC267F", "#FFB000", "#785EF0", "#FF832B", "black","#648FFF", "#DC267F", "#FFB000", "#785EF0", "#FF832B", "black","#648FFF", "#DC267F", "#FFB000")
# patterns <- c(NA, NA, NA, NA, NA, NA,"stripe", "stripe", "stripe", "stripe", "stripe", "stripe", "circle", "circle", "circle")
# ggplot(data, aes(x = Pathogen, fill = Antibiotic)) +
#   geom_bar_pattern(aes(fill = Antibiotic, pattern=Antibiotic),
#     pattern_density = 0.3, pattern_color = "#00000000", pattern_spacing = 0.005) +
#   scale_fill_manual(values = colors) +
#   scale_pattern_manual(values = patterns) +
#   labs(title = "Total Number of Entries by ATC Class and Pathogen", x = "Pathogen", y = "Number of Entries") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# # save plot
# ggsave("data_entries_by_pathogen_and_ATC3_full.png")

# counts <- data %>%
#     group_by(Pathogen, Antibiotic) %>%
#     summarise(Counts = sum(Weight, na.rm = TRUE)) %>%
#     arrange(desc(Counts))
# for (pathogen in sort(unique(counts$Antibiotic))) {
#     print(paste("Antibiotic:", pathogen))
#     # print in alphabetical order of Antibiotic
#     counts_subset <- counts[counts$Antibiotic == pathogen,]
#     counts_subset <- counts_subset[order(counts_subset$Pathogen),]
#     print(counts_subset)
#     print(sum(counts$Counts[counts$Antibiotic == pathogen]))
# }

# load global consumption by year, country, and ATC3
consumption <- read.csv("antibiotic_consumption_by_ATC3.csv")
global_consumption <- consumption[consumption$Location == "Global",]
global_consumption <- global_consumption[global_consumption$Year == "2018",]
global_consumption <- select(global_consumption, ATC.level.3.class, Antibiotic.consumption..DDD.1.000.day.)
global_consumption <- global_consumption %>%
    rename(
        Antibiotic = ATC.level.3.class,
        Global.Consumption = Antibiotic.consumption..DDD.1.000.day.
    )
global_consumption$Antibiotic <- sub("-.*", "", global_consumption$Antibiotic)

# # scale Consumption by global mean for year and antibiotic
data <- data %>%
    left_join(global_consumption, by = c("Antibiotic")) %>%
    mutate(Consumption = Consumption / Global.Consumption) %>%
    select(-Global.Consumption)
data <- na.omit(data)

data_HIC <- data_HIC %>%
    left_join(global_consumption, by = c("Antibiotic")) %>%
    mutate(Consumption = Consumption / Global.Consumption) %>%
    select(-Global.Consumption)
data_HIC <- na.omit(data_HIC)

data_LMIC <- data_LMIC %>%
    left_join(global_consumption, by = c("Antibiotic")) %>%
    mutate(Consumption = Consumption / Global.Consumption) %>%
    select(-Global.Consumption)
data_LMIC <- na.omit(data_LMIC)

# log transform
data$Consumption <- log(data$Consumption + 1)
data$Resistance <- log(data$Resistance + 1)
data$Weight <- data$Weight / max(data$Weight, na.rm = TRUE)

data_HIC$Consumption <- log(data_HIC$Consumption + 1)
data_HIC$Resistance <- log(data_HIC$Resistance + 1)
data_HIC$Weight <- data_HIC$Weight / max(data_HIC$Weight, na.rm = TRUE)

data_LMIC$Consumption <- log(data_LMIC$Consumption + 1)
data_LMIC$Resistance <- log(data_LMIC$Resistance + 1)
data_LMIC$Weight <- data_LMIC$Weight / max(data_LMIC$Weight, na.rm = TRUE)

# compute year-on-year changes in resistance in each country for each antibiotic and pathogen
data <- data %>%
    group_by(ISO3, Antibiotic, Pathogen) %>%
    arrange(Year) %>%
    mutate(Resistance = (Resistance - lag(Resistance)) / lag(Resistance)) %>%
    ungroup()


# require(boot)
# require(car)
# # table of gradient with confidence intervals for separate linear model for each antibiotic - pathogen combination
# gradients <- c()
# # calculate confidence intervals
# conf_intervals <- c()
# pathogens <- c()
# abs <- c()
# r_squareds <- c()
# results_variation <- c()
# bootstraps <- numeric()
# for (antibiotic in sort(unique(data$Antibiotic))) {
#     print(antibiotic)
#     for (pathogen in sort(unique(data$Pathogen))) {
#         data_subset <- data[data$Pathogen == pathogen & data$Antibiotic == antibiotic,]
#         # filter subset to just Resistance, Consumption, and country covariates
#         # data_subset <- data_subset %>%
#         #     select(Consumption, Resistance, Pathogen, Antibiotic, all_of(country_covariates))
#         # remove NAs
#         data_subset <- data_subset[complete.cases(data_subset),]
#         # remove infs
#         data_subset <- data_subset[!is.infinite(data_subset$Consumption),]
#         data_subset <- data_subset[!is.infinite(data_subset$Resistance),]
#         if (nrow(data_subset) <= 1) {
#             next
#         }
#         model <- lm(Resistance ~ Consumption
#         + PC1 + PC2 + PC3 + GDP + Year
#         #  + National.Income.per.capita..log..Relative.to.Global.Value + Mean.Years.Schooling.Relative.to.Global.Value + Mortality.rate.attributed.to.unsafe.water..unsafe.sanitation.and.lack.of.hygiene.from.diarrhoea..intestinal.nematode.infections..malnutrition.and.acute.respiratory.infections..deaths.per.100.000.population..Relative.to.Global.Value + Proportion.of.population.practicing.open.defecation..across.all.locations.....Relative.to.Global.Value + Proportion.of.population.using.basic.drinking.water.services..across.all.locations.....Relative.to.Global.Value + Proportion.of.population.using.basic.sanitation.services..across.all.locations.....Relative.to.Global.Value + Proportion.of.population.using.safely.managed.drinking.water.services..across.all.locations.....Relative.to.Global.Value + Proportion.of.population.using.safely.managed.sanitation.services..across.all.locations.....Relative.to.Global.Value + Proportion.of.population.with.basic.handwashing.facilities.on.premises..across.all.locations.....Relative.to.Global.Value + Proportion.of.schools.with.access.to.basic.drinking.water..primary.schools.only.....Relative.to.Global.Value + Proportion.of.schools.with.access.to.single.sex.basic.sanitation..primary.schools.only.....Relative.to.Global.Value + Proportion.of.schools.with.basic.handwashing.facilities..primary.schools.only.....Relative.to.Global.Value + Proportion.of.wastewater.treated..across.all.locations.and.activities.....Relative.to.Global.Value + Population + Life.Expectancy.Relative.to.Global.Value
#         , data = data_subset,
#         weights = Weight)
#         r_squared <- summary(model)$r.squared
#         # Calculate variation explained by each variable individually
#         variables <- c("Consumption", "PC1", "PC2", "PC3", "GDP", "Year")
#         variation_explained <- c()
        
#         for (var in variables) {
#             # Create formula excluding the current variable
#             other_vars <- variables[variables != var]
#             formula_str <- paste("Resistance ~", paste(other_vars, collapse = " + "))
            
#             model_null <- lm(as.formula(formula_str), data = data_subset, weights = Weight)
#             r_squared_null <- summary(model_null)$r.squared
#             variation_explained <- c(variation_explained, r_squared - r_squared_null)
#         }
        
#         # Store all variation explained values
#         variation_explained_df <- data.frame(t(setNames(variation_explained, variables)))
#         results_variation <- rbind(results_variation, variation_explained_df)
#         r_squareds <- c(r_squareds, r_squared)
#         # print(paste("Antibiotic:", antibiotic, "Pathogen:", pathogen, "Coefficient:", coef(model)["Consumption"]))
#         if (is.na(coef(model)["Consumption"])) {
#             next
#         }
#         gradients <- c(gradients, coef(model)["Consumption"])
#         conf_intervals <- c(conf_intervals, confint(model)["Consumption",])
#         pathogens <- c(pathogens, pathogen)
#         abs <- c(abs, antibiotic)
#         bs <- Boot(model, R=1000)
#         # data frame with columns Pathogen, Antibiotic, Gradient (for Consumption), 1000 rows
#         outdf <- data.frame(Pathogen = pathogen, Antibiotic = antibiotic, Gradient = bs$t)
#         bootstraps <- rbind(bootstraps, outdf)
#     }
# }
# conf_intervals <- matrix(conf_intervals, nrow = length(gradients), ncol = 2, byrow = TRUE)
# # matrix with columns: Pathogen, Antibiotic, Gradient, Lower CI, Upper CI
# results <- data.frame(Antibiotic = abs,
#                     Pathogen = pathogens,
#                     Response = gradients,
#                     Lower_CI = conf_intervals[,1],
#                     Upper_CI = conf_intervals[,2],
#                     R_squared = r_squareds,
#                     Variation_Explained = results_variation
#                     )
# # print(results)
# # # save results in csv
# write.csv(results, paste0("Outputs/database_gradients_pathogen_ATC3_PCA_joelike_weighted_lagged.csv"), row.names = FALSE)
# write.csv(bootstraps, paste0("Outputs/database_gradients_bootstraps_pathogen_ATC3_PCA_joelike_weighted_lagged.csv"), row.names = FALSE)

# # # # # remove S. pneumoniae, S. aureus, E. faecalis, E. faecium, and N. gonorrhoeae
# # # # results <- results %>%
# # # #     filter(Pathogen != "S. pneumoniae" & Pathogen != "S. aureus" & Pathogen != "E. faecalis" & Pathogen != "E. faecium" & Pathogen != "N. gonorrhoeae")
# # # # filter to only S. pneumoniae, S. aureus, E. faecalis, E. faecium, and N. gonorrhoeae
# # # # results <- results %>%
# # # #     filter(Pathogen == "S. pneumoniae" | Pathogen == "S. aureus" | Pathogen == "E. faecalis" | Pathogen == "E. faecium" | Pathogen == "N. gonorrhoeae")
# # # # transform Response into character with confidence intervals in parentheses
# # # results <- results %>%
# # #     mutate(across(Response, ~paste(signif(., 2), " (", signif(Lower_CI, 2), ", ", signif(Upper_CI, 2), ")", sep = ""))
# # #            )
# # # # remove Lower_CI and Upper_CI columns
# # # results <- results %>%
# # #     select(-Lower_CI, -Upper_CI)
# # # # remove index column
# # # rownames(results) <- NULL
# # # # print in latex format
# # # results <- results %>%
# # #     pivot_wider(names_from = Pathogen, values_from = Response)
# # # # print(results)
# # # results <- results %>%
# # #     select(Antibiotic, everything())

# # # print(xtable::xtable(results), type = "latex")

# # # require(lme4)
# # # require("lmerTest")
# # # # crossed random effects since antibiotics can be used for multiple pathogens
# # # model_crossed <- lmer(Resistance ~ Consumption + (Consumption||Pathogen) + (Consumption||Antibiotic), data = data)
# # # print(summary(model_crossed))
# # # require(nlme)
# # # model_nested <- lme(Resistance ~ Consumption,random=~1|Pathogen/Antibiotic,data=data,na.action = na.exclude)
# # # print(summary(model_nested)$tTable)

# get gradient and p-value for each pathogen

require(lme4)
income <- "all"
for (income in c("all")) {
    print(income)
    data_ <- if (income == "HIC") {
        data_HIC
    } else if (income == "LMIC") {
        data_LMIC
    } else {
        data
    }
    class_gradients <- c()
    class_intercepts <- c()
    class_lowerCI <- c()
    class_upperCI <- c()
    classes <- sort(unique(data_$Antibiotic))
    bootstraps <- numeric()
    # print(classes)÷
    for (class in classes) {
        print(class)
        class_data <- data_[data_$Antibiotic == class,]
        class_data <- class_data[!is.infinite(class_data$Consumption),]
        class_data <- class_data[!is.infinite(class_data$Resistance),]
        # make sure there is more than one Antibiotic in data
        if (length(unique(class_data$Pathogen)) > 1) {
            model <- lmer(Resistance ~ Consumption + (Consumption||Pathogen)
            + PC1 + PC2 + PC3 + GDP + Year
            #  + National.Income.per.capita..log..Relative.to.Global.Value + Mean.Years.Schooling.Relative.to.Global.Value + Mortality.rate.attributed.to.unsafe.water..unsafe.sanitation.and.lack.of.hygiene.from.diarrhoea..intestinal.nematode.infections..malnutrition.and.acute.respiratory.infections..deaths.per.100.000.population..Relative.to.Global.Value + Proportion.of.population.practicing.open.defecation..across.all.locations.....Relative.to.Global.Value + Proportion.of.population.using.basic.drinking.water.services..across.all.locations.....Relative.to.Global.Value + Proportion.of.population.using.basic.sanitation.services..across.all.locations.....Relative.to.Global.Value + Proportion.of.population.using.safely.managed.drinking.water.services..across.all.locations.....Relative.to.Global.Value + Proportion.of.population.using.safely.managed.sanitation.services..across.all.locations.....Relative.to.Global.Value + Proportion.of.population.with.basic.handwashing.facilities.on.premises..across.all.locations.....Relative.to.Global.Value + Proportion.of.schools.with.access.to.basic.drinking.water..primary.schools.only.....Relative.to.Global.Value + Proportion.of.schools.with.access.to.single.sex.basic.sanitation..primary.schools.only.....Relative.to.Global.Value + Proportion.of.schools.with.basic.handwashing.facilities..primary.schools.only.....Relative.to.Global.Value + Proportion.of.wastewater.treated..across.all.locations.and.activities.....Relative.to.Global.Value + Population + Life.Expectancy.Relative.to.Global.Value
            , data = class_data,
            weights = Weight)
            print(model)
            # if model is singular, exclude class
            if (isSingular(model)) {
                print(paste("Model is singular for class", class, "- excluding"))
                classes <- classes[classes != class]
                next
            }
            # print(summary(model))
            class_gradients <- c(class_gradients, summary(model)$coefficients["Consumption",1])
            # get intercepts
            class_intercepts <- c(class_intercepts, summary(model)$coefficients["(Intercept)",1])
            # intervals <- confint(model, method="boot")
            # bs <- Boot(model, R=1000)
            bs <- bootMer(model, FUN = function(x) fixef(x)["Consumption"], nsim = 1000)
            intervals <- confint(bs)
            outdf <- data.frame(Antibiotic = class, Gradient = bs$t)
            bootstraps <- rbind(bootstraps, outdf)
            class_lowerCI <- c(class_lowerCI, intervals["Consumption","2.5 %"])
            class_upperCI <- c(class_upperCI, intervals["Consumption","97.5 %"])
        } else {
            model <- lm(Resistance ~ Consumption
            + PC1 + PC2 + PC3 + GDP + Year
            #  + National.Income.per.capita..log..Relative.to.Global.Value + Mean.Years.Schooling.Relative.to.Global.Value + Mortality.rate.attributed.to.unsafe.water..unsafe.sanitation.and.lack.of.hygiene.from.diarrhoea..intestinal.nematode.infections..malnutrition.and.acute.respiratory.infections..deaths.per.100.000.population..Relative.to.Global.Value + Proportion.of.population.practicing.open.defecation..across.all.locations.....Relative.to.Global.Value + Proportion.of.population.using.basic.drinking.water.services..across.all.locations.....Relative.to.Global.Value + Proportion.of.population.using.basic.sanitation.services..across.all.locations.....Relative.to.Global.Value + Proportion.of.population.using.safely.managed.drinking.water.services..across.all.locations.....Relative.to.Global.Value + Proportion.of.population.using.safely.managed.sanitation.services..across.all.locations.....Relative.to.Global.Value + Proportion.of.population.with.basic.handwashing.facilities.on.premises..across.all.locations.....Relative.to.Global.Value + Proportion.of.schools.with.access.to.basic.drinking.water..primary.schools.only.....Relative.to.Global.Value + Proportion.of.schools.with.access.to.single.sex.basic.sanitation..primary.schools.only.....Relative.to.Global.Value + Proportion.of.schools.with.basic.handwashing.facilities..primary.schools.only.....Relative.to.Global.Value + Proportion.of.wastewater.treated..across.all.locations.and.activities.....Relative.to.Global.Value + Population + Life.Expectancy.Relative.to.Global.Value
            , data = class_data,
            weights = Weight)
            class_gradients <- c(class_gradients, summary(model)$coefficients["Consumption",1])
            # get intercepts
            class_intercepts <- c(class_intercepts, summary(model)$coefficients["(Intercept)",1])
            bs <- Boot(model, R=1000)
            outdf <- data.frame(Antibiotic = class, Consumption = bs$t[,"Consumption"])
            bootstraps <- rbind(bootstraps, outdf)
            intervals <- confint(model)
            class_lowerCI <- c(class_lowerCI, intervals["Consumption",1])
            class_upperCI <- c(class_upperCI, intervals["Consumption",2])
        }
    }
    class_gradients <- setNames(class_gradients, classes)
    class_intercepts <- setNames(class_intercepts, classes)
    class_lowerCI <- setNames(class_lowerCI, classes)
    class_upperCI <- setNames(class_upperCI, classes)
    write.csv(class_gradients, paste0("Outputs/database_gradients_ATC3_PCA_joelike_weighted_",income,"_lagged.csv"), row.names = TRUE)
    write.csv(class_lowerCI, paste0("Outputs/database_lowerCI_ATC3_PCA_joelike_weighted_",income,"_lagged.csv"), row.names = TRUE)
    write.csv(class_upperCI, paste0("Outputs/database_upperCI_ATC3_PCA_joelike_weighted_",income,"_lagged.csv"), row.names = TRUE)
    write.csv(bootstraps, paste0("Outputs/database_gradients_bootstraps_ATC3_PCA_joelike_weighted_",income,"_lagged.csv"), row.names = TRUE)
}

# # results <- read.csv("Outputs/database_gradients_pathogen_ATC3_PCA.csv", stringsAsFactors = FALSE)
# # class_gradients <- setNames(as.vector(read.csv("Outputs/database_gradients_ATC3_PCA.csv", stringsAsFactors = FALSE)[,2]), 
# #                             read.csv("Outputs/database_gradients_ATC3_PCA.csv", stringsAsFactors = FALSE)[,1])

# library(ggallin)

# ab_classes <- unique(results$Antibiotic)

# # Create a mapping data frame
# segment_data <- data.frame(
#     Antibiotic = ab_classes,
#     Value = class_gradients[ab_classes],
#     Lower_CI = class_lowerCI[ab_classes],
#     Upper_CI = class_upperCI[ab_classes]
# )

# pathogen_list = c("Acinetobacter spp.", "E. faecalis", "E. faecium", "Enterococcus spp.", "E. coli", "H. influenzae", "K. pneumoniae", "Morganella spp.", "N. gonorrhoeae", "P. aeruginosa", "Salmonella spp.", "S. aureus", "S. agalactiae", "S. pneumoniae", "S. pyogenes")
# # in hospital data set need to expand list
# # pathogen_list = c(pathogen_list, "S. epidermidis","Enterobacteriaceae","Proteus spp.","Streptococcus spp.","Enterobacter spp.","Gram-negatives")
# pathogen_list = sort(pathogen_list)
# colors = c("#648FFF","#DC267F","#FFB000","#785EF0","#FF832B","#000000","#648FFF","#DC267F","#FFB000","#785EF0","#FF832B", "#000000","#648FFF","#DC267F","#FFB000")
# # # in hospital data set need to expand colors
# # colors = c("#648FFF","#DC267F","#FFB000","#785EF0","#FF832B","#000000","#648FFF","#DC267F","#FFB000","#785EF0","#FF832B", "#000000","#648FFF","#DC267F","#FFB000")

# # colors = c("#648FFF","#DC267F","#FFB000",
# # "#785EF0","#648FFF","#DC267F",
# # "#FF832B","#000000","#648FFF",
# # "#DC267F","#FFB000","#785EF0",
# # "#FFB000", "#FF832B","#000000",
# # "#785EF0","#648FFF","#DC267F",
# # "#FFB000","#FF832B","#000000")
# colors = setNames(colors, pathogen_list)
# color_vector = c()
# # get position of pathogen in pathogen_list and assign color
# for (pathogen in results$Pathogen) {
#     color_vector = c(color_vector, colors[pathogen])
# }
# shapes = c(18,18,18,18,18,18,16,16,16,16,16,16,15,15,15)
# # in hospital data set need to expand shapes
# # shapes = c(18,18,18,18,15,15,15,18,18,16,16,15,16,16,15,16,16,15,15,15,17)
# # shapes = c(18,18,18,
# # 18,17,17,
# # 18,18,16,
# # 16,16,16,
# # 17,16,16,
# # 17,15,15,
# # 15,17,17)
# shapes = setNames(shapes, pathogen_list)
# shape_vector = c()
# for (pathogen in results$Pathogen) {
#     shape_vector = c(shape_vector, shapes[pathogen])
# }
# # Plot elasticities with confidence intervals, grouped by Antibiotic. Random effects model plotted using errorbar
# # for (i in 1:nrow(results)){
# #     results$Antibiotic[i] <- atc_names[[results$Antibiotic[i]]]
# # }
# # for (i in 1:nrow(segment_data)){
# #     segment_data$Antibiotic[i] <- atc_names[[segment_data$Antibiotic[i]]]
# # }
# # # Define the custom order for Antibiotic
# # custom_order <- c(  "Tetracyclines",
# #   "Glycopeptides and Lipopeptides",
# #   "Penicillins",
# #   "Other Beta-Lactams",
# #   "Sulfonamides and Trimethoprim",
# #   "Macrolides, Lincosamides and Streptogramins",
# #   "Aminoglycosides",
# #   "Quniolones") # Replace with your actual antibiotic names

# # # Convert Antibiotic to a factor with the custom order
# # results$Antibiotic <- factor(results$Antibiotic, levels = custom_order)
# ggplot(results, aes(x = Antibiotic, y = Response, color = Pathogen, shape = Pathogen)) +
#     geom_errorbar(data = segment_data,
#                                 aes(x = Antibiotic, y = Value, ymin = Lower_CI, ymax = Upper_CI), 
#                                 width = 0, inherit.aes = FALSE, color = "#eaeaea", linewidth = 15) +
#         geom_errorbar(data = segment_data,
#                                                                 aes(x = Antibiotic, y = Value, ymin = Value, ymax = Value), 
#                                                                 width = 0.6, inherit.aes = FALSE, color = "gray", linewidth = 1) +
#         geom_point(position = position_dodge(width = 0.5), size = 2.5) +
#         geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0, position = position_dodge(width = 0.5)) +
#         geom_hline(yintercept = 0, color = "gray", linetype = "dashed", alpha = 0.5) +
#         scale_color_manual(values = color_vector, name = NULL) +
#         scale_shape_manual(values = shape_vector, name = NULL) +
#         scale_y_continuous(limits = c(-1, 2.4)) +  # Limit y scale between -0.5 and 1
#         theme_minimal() +
#         theme(panel.background = element_rect(fill = "white")) +
#         theme(panel.grid = element_blank()) +
#         labs(x = "Antibiotic Class", y = "Gradient") +
#         ggtitle("Elasticity by Antibiotic Class and Pathogen") +
#         theme(plot.title = element_text(hjust = 0.5)) +
#         theme(axis.title.x = element_text(size = 11)) +
#         theme(axis.title.y = element_text(size = 11)) +
#         theme(axis.text.x = element_text(size = 9)) +
#         theme(axis.text.y = element_text(size = 9)) +
#         theme(legend.position = "bottom", legend.text = element_text(size = 8))  # Reduced legend text size
#         # scale_y_continuous(trans = pseudolog10_trans, breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000))

# # save the plot
# ggsave("database_drug_PC3controlled_joelike_60cutoff_no0resistance_no0overlap_slide.png",width=13.3,height=7.5,units="in")

# # # ## combined Nagorsen plot            
# library(ggallin)
# # # # # Plot elasticities with confidence intervals, grouped by Antibiotic
# # # # # First, get the unique pathogen levels in the correct order
# results_h2h <- read.csv("Outputs/Nagorsen_gradients_pathogen_ATC3_PCA_joelike_hospital_to_all.csv", stringsAsFactors = FALSE)
# results_c2c <- read.csv("Outputs/Nagorsen_gradients_pathogen_ATC3_PCA_joelike_community_to_all.csv", stringsAsFactors = FALSE)
# class_gradients_h2h <- setNames(as.vector(read.csv("Outputs/Nagorsen_gradients_ATC3_PCA_joelike_hospital_to_all.csv", stringsAsFactors = FALSE)[,2]), 
#                             read.csv("Outputs/Nagorsen_gradients_ATC3_PCA_joelike_hospital_to_all.csv", stringsAsFactors = FALSE)[,1])
# class_gradients_c2c <- setNames(as.vector(read.csv("Outputs/Nagorsen_gradients_ATC3_PCA_joelike_community_to_all.csv", stringsAsFactors = FALSE)[,2]), 
#                             read.csv("Outputs/Nagorsen_gradients_ATC3_PCA_joelike_community_to_all.csv", stringsAsFactors = FALSE)[,1])
# class_lowerCI_h2h <- setNames(as.vector(read.csv("Outputs/Nagorsen_lowerCI_ATC3_PCA_joelike_hospital_to_all.csv", stringsAsFactors = FALSE)[,2]), 
#                             read.csv("Outputs/Nagorsen_lowerCI_ATC3_PCA_joelike_hospital_to_all.csv", stringsAsFactors = FALSE)[,1])
# class_lowerCI_c2c <- setNames(as.vector(read.csv("Outputs/Nagorsen_lowerCI_ATC3_PCA_joelike_community_to_all.csv", stringsAsFactors = FALSE)[,2]),
#                             read.csv("Outputs/Nagorsen_lowerCI_ATC3_PCA_joelike_community_to_all.csv", stringsAsFactors = FALSE)[,1])
# class_upperCI_h2h <- setNames(as.vector(read.csv("Outputs/Nagorsen_upperCI_ATC3_PCA_joelike_hospital_to_all.csv", stringsAsFactors = FALSE)[,2]),
#                             read.csv("Outputs/Nagorsen_upperCI_ATC3_PCA_joelike_hospital_to_all.csv", stringsAsFactors = FALSE)[,1])
# class_upperCI_c2c <- setNames(as.vector(read.csv("Outputs/Nagorsen_upperCI_ATC3_PCA_joelike_community_to_all.csv", stringsAsFactors = FALSE)[,2]),
#                             read.csv("Outputs/Nagorsen_upperCI_ATC3_PCA_joelike_community_to_all.csv", stringsAsFactors = FALSE)[,1])
# # # ab_classes <- unique(results$Antibiotic)
# ab_classes <- unique(c(names(class_gradients_h2h), names(class_gradients_c2c)))
# print(results_h2h)
# # # Create a mapping data frame
# # segment_data <- data.frame(
# #     Antibiotic = ab_classes,
# #     Value = class_gradients[ab_classes],
# #     Lower_CI = class_lowerCI[ab_classes],
# #     Upper_CI = class_upperCI[ab_classes]
# # )
# segment_data_h2h <- data.frame(
#     Antibiotic = ab_classes,
#     Value = class_gradients_h2h[ab_classes],
#     Lower_CI = class_lowerCI_h2h[ab_classes],
#     Upper_CI = class_upperCI_h2h[ab_classes]
# )
# segment_data_c2c <- data.frame(
#     Antibiotic = ab_classes,
#     Value = class_gradients_c2c[ab_classes],
#     Lower_CI = class_lowerCI_c2c[ab_classes],
#     Upper_CI = class_upperCI_c2c[ab_classes]
# )


# pathogen_list = c("Acinetobacter spp.", "E. coli", "E. faecium", "Enterobacter spp.", "H. influenzae", 
# "K. pneumoniae", "M. tuberculosis", "N. meningitidis", 
# "P. aeruginosa", "S. aureus", "S. pneumoniae", 
# "S. pyogenes", "Shigella spp.")
# pathogen_list = sort(pathogen_list)
# print(pathogen_list)
# colors = c("#648FFF","#DC267F","#FFB000","#648FFF","#000000",
# "#648FFF", "#DC267F", "#FFB000",
# "#785EF0","#000000","#648FFF",
# "#DC267F","#FFB000")

# colors = setNames(colors, pathogen_list)
# color_vector = c()
# # get position of pathogen in pathogen_list and assign color
# all_pathogens = unique(c(results_h2h$Pathogen, results_c2c$Pathogen))
# print(all_pathogens)
# # all_pathogens = unique(data$Pathogen)
# for (pathogen in all_pathogens) {
#     color_vector = c(color_vector, colors[pathogen])
# }

# shapes = c(18,18,18,17,18,
# 16,17,17,
# 16,16,15,
# 15,17)
# shapes = setNames(shapes, pathogen_list)
# shape_vector = c()
# for (pathogen in all_pathogens) {
#     shape_vector = c(shape_vector, shapes[pathogen])
# }
# # # Plot elasticities with confidence intervals, grouped by Antibiotic. Random effects model plotted using errorbar
# # # for (i in 1:nrow(results)){
# # #     results$Antibiotic[i] <- atc_names[[results$Antibiotic[i]]]
# # # }
# # # for (i in 1:nrow(segment_data)){
# # #     segment_data$Antibiotic[i] <- atc_names[[segment_data$Antibiotic[i]]]
# # # }
# # # # Define the custom order for Antibiotic
# # # custom_order <- c(  "Tetracyclines",
# # #   "Glycopeptides and Lipopeptides",
# # #   "Penicillins",
# # #   "Other Beta-Lactams",
# # #   "Sulfonamides and Trimethoprim",
# # #   "Macrolides, Lincosamides and Streptogramins",
# # #   "Aminoglycosides",
# # #   "Quniolones") # Replace with your actual antibiotic names

# # # Convert Antibiotic to a factor with the custom order
# # results$Antibiotic <- factor(results$Antibiotic, levels = custom_order)
# # add a row of nas to results_c2c
# # remove results with NA confidence intervals
# results_h2h <- results_h2h[!is.na(results_h2h$Lower_CI) & !is.na(results_h2h$Upper_CI),]
# results_c2c <- results_c2c[!is.na(results_c2c$Lower_CI) & !is.na(results_c2c$Upper_CI),]
# # segment_data_c2c <- segment_data_c2c[!is.na(segment_data_c2c$Lower_CI) & !is.na(segment_data_c2c$Upper_CI),]
# # segment_data_h2h <- segment_data_h2h[!is.na(segment_data_h2h$Lower_CI) & !is.na(segment_data_h2h$Upper_CI),]
# print(segment_data_h2h)
# results_c2c <- rbind(results_c2c, data.frame(Antibiotic = "J01M", Pathogen = "S. aureus", Response = NA, Lower_CI = NA, Upper_CI = NA))
# results_c2c <- rbind(results_c2c, data.frame(Antibiotic = "J01M", Pathogen = "Enterobacter spp.", Response = NA, Lower_CI = NA, Upper_CI = NA))
# results_c2c <- rbind(results_c2c, data.frame(Antibiotic = "J01M", Pathogen = "E. faecium", Response = NA, Lower_CI = NA, Upper_CI = NA))
# plot_h2h <- ggplot(results_h2h, aes(x = Antibiotic, y = Response, color = Pathogen, shape = Pathogen)) +
#     geom_errorbar(data = segment_data_h2h,
#                                 aes(x = Antibiotic, y = Value, ymin = Lower_CI, ymax = Upper_CI), 
#                                 width = 0, inherit.aes = FALSE, color = "#eaeaea", linewidth = 15) +
#         geom_errorbar(data = segment_data_h2h,
#                                                                 aes(x = Antibiotic, y = Value, ymin = Value, ymax = Value), 
#                                                                 width = 0.6, inherit.aes = FALSE, color = "gray", linewidth = 1) +
#         geom_point(position = position_dodge(width = 0.5), size = 2.5) +
#         geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0, position = position_dodge(width = 0.5)) +
#         # geom_vline(xintercept = seq(1.5, length(unique(results$Antibiotic)), by = 1), color = "gray", linetype = "dashed", alpha = 0.5) +
#         geom_hline(yintercept = 0, color = "gray", linetype = "dashed", alpha = 0.5) +
#         scale_color_manual(values = color_vector, name = NULL) +
#         scale_shape_manual(values = shape_vector, name = NULL) +
#         theme_minimal() +
#         theme(panel.background = element_rect(fill = "white")) +
#         theme(panel.grid = element_blank()) +
#         labs(x = "Antibiotic Class", y = "Gradient") +
#         theme(plot.title = element_text(hjust = 0.5)) +
#         theme(axis.title.x = element_text(size = 11)) +
#         theme(axis.title.y = element_text(size = 11)) +
#         theme(axis.text.x = element_text(size = 9)) +
#         theme(axis.text.y = element_text(size = 9)) +
#         theme(legend.position = "bottom", legend.text = element_text(size = 8)) +  # Reduced legend text size
#         scale_y_continuous(limits = c(-5.7,9.3), trans = pseudolog10_trans)
# plot_c2c <- ggplot(results_c2c, aes(x = Antibiotic, y = Response, color = Pathogen, shape = Pathogen)) +
#     geom_errorbar(data = segment_data_c2c,
#                                 aes(x = Antibiotic, y = Value, ymin = Lower_CI, ymax = Upper_CI), 
#                                 width = 0, inherit.aes = FALSE, color = "#eaeaea", linewidth = 15) +
#         geom_errorbar(data = segment_data_c2c,
#                                                                 aes(x = Antibiotic, y = Value, ymin = Value, ymax = Value), 
#                                                                 width = 0.6, inherit.aes = FALSE, color = "gray", linewidth = 1) +
#         geom_point(position = position_dodge(width = 0.5), size = 2.5) +
#         geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0, position = position_dodge(width = 0.5)) +
#         # geom_vline(xintercept = seq(1.5, length(unique(results$Antibiotic)), by = 1), color = "gray", linetype = "dashed", alpha = 0.5) +
#         geom_hline(yintercept = 0, color = "gray", linetype = "dashed", alpha = 0.5) +
#         scale_color_manual(values = color_vector, name = NULL) +
#         scale_shape_manual(values = shape_vector, name = NULL) +
#         theme_minimal() +
#         theme(panel.background = element_rect(fill = "white")) +
#         theme(panel.grid = element_blank()) +
#         labs(x = "Antibiotic Class", y = "Gradient") +
#         theme(plot.title = element_text(hjust = 0.5)) +
#         theme(axis.title.x = element_text(size = 11)) +
#         theme(axis.title.y = element_text(size = 11)) +
#         theme(axis.text.x = element_text(size = 9)) +
#         theme(axis.text.y = element_text(size = 9)) +
#         theme(legend.position = "bottom", legend.text = element_text(size = 8)) +  # Reduced legend text size
#         scale_y_continuous(limits = c(-5.7,9.3), trans = pseudolog10_trans)
# library(ggpubr)
# ggarrange(plot_c2c + theme(legend.position="none"), plot_h2h + theme(legend.position="none"),
#           ncol = 1, nrow = 2,
#           labels = c("A", "B"),
#           common.legend = TRUE, legend = "bottom",
#           align = "v", heights = c(1, 1)) +
#     theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
#     theme(plot.title = element_text(hjust = 0.5))

# # save the plot
# ggsave("Nagorsen_drug_PC3controlled_joelike_combined_to_all.png",width=6.5,height=7.5,units="in")




# for (pathogen in unique(data$Pathogen)) {
#     # filter data for pathogen
#     data_subset <- data[data$Pathogen == pathogen,]
#     ## multipanel plot with separate plot for each pathogen, with Antibiotic random effects
#     ggplot(data_subset, aes(x = Consumption, y = Resistance, color = Pathogen, shape = Pathogen)) +
#         geom_point() +
#         theme_minimal() +
#         theme(panel.background = element_rect(fill = "white")) +
#         labs(x = "Relative consumption", y = "Percent resistant") +
#         theme(plot.title = element_text(hjust = 0.5)) +
#         theme(axis.title.x = element_text(size = 11)) +
#         theme(axis.title.y = element_text(size = 11)) +
#         theme(axis.text.x = element_text(size = 9)) +
#         theme(axis.text.y = element_text(size = 9)) +
#         theme(legend.position = "bottom", legend.text = element_text(size = 8)) +
#         facet_wrap(~Antibiotic, scales = "free") +
#         scale_color_manual(values = color_vector, name = NULL) +
#         scale_shape_manual(values = shape_vector, name = NULL) +
#         scale_x_log10() +  # Logarithmic scale for x-axis
#         scale_y_log10()    # Logarithmic scale for y-axis

#     # save the multipanel plot
#     ggsave(paste0("database_scatter_",pathogen,"_loglog.png"),width=6.5,height=6.5,units="in")
# }