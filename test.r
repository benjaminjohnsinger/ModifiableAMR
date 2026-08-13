data <- read.csv("summed_data_new.csv")
head(data)


source("utils.r")
# use the iso3_ihme_mapping to give a lending group to each row based on ISO3 code
data$Lending.Group <- iso3_ihme_mapping$lending_group[match(data$ISO3, iso3_ihme_mapping$iso3)]
data$lower_ihme_region <- iso3_ihme_mapping$lower_ihme_region[match(data$ISO3, iso3_ihme_mapping$iso3)]

# Sum the total number of isolates by IHME region.
isolates_by_region <- aggregate(Total.Isolates ~ lower_ihme_region, data = data, sum, na.rm = TRUE)
isolates_by_region <- isolates_by_region[order(isolates_by_region$Total.Isolates, decreasing = TRUE), ]
print(isolates_by_region)

isolates_by_year <- aggregate(Total.Isolates ~ Year, data = data, sum, na.rm = TRUE)
isolates_by_year <- isolates_by_year[order(isolates_by_year$Year), ]
print(isolates_by_year)
isolates_by_year <- isolates_by_year[isolates_by_year$Year < 2017, ]
# plot isolates by year with a fitted line and print estimated slope
require(ggplot2)
plot_isolates_by_year <- ggplot(isolates_by_year, aes(x = Year, y = Total.Isolates)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Total Isolates by Year", x = "Year", y = "Total Isolates") +
  theme_minimal()
png("plot_isolates_by_year.png", width = 6.5, height = 3, units = "in", res = 300)
print(plot_isolates_by_year)
dev.off()


# isolates by region and year
isolates_by_region_and_year <- aggregate(Total.Isolates ~ lower_ihme_region + Year, data = data, sum, na.rm = TRUE)
isolates_by_region_and_year <- isolates_by_region_and_year[order(isolates_by_region_and_year$lower_ihme_region, isolates_by_region_and_year$Year), ]
# table with years as columns
isolates_by_region_and_year_table <- reshape(isolates_by_region_and_year, idvar = "lower_ihme_region", timevar = "Year", direction = "wide")
# sort columns by year
year_cols <- as.numeric(gsub("Total.Isolates.", "", names(isolates_by_region_and_year_table)[-1]))
year_order <- order(year_cols)
isolates_by_region_and_year_table <- isolates_by_region_and_year_table[, c(1, year_order + 1)]
print(isolates_by_region_and_year_table)



# # filter to high income countries only
# data <- data %>% filter(Lending.Group != "High income")

library(tidyr)
# for each antibiotic class and pathogen, calculate the proportion of of times when Percent.Resistant.Isolates goes up from one year to the next, and the proportion of times when Percent.Resistant.Isolates goes down from one year to the next.
# plot this as a heatmap
plotting_data <- data %>%
  group_by(ATC.Class, Pathogen) %>%
  filter(n() > 10) %>% #then filter out if sum of Total.Isolates is less than 100
  filter(sum(Total.Isolates, na.rm = TRUE) > 100) %>%
  arrange(Year) %>%
  summarise(
    prop_up = mean(diff(Percent.Resistant.Isolates) > 0, na.rm = TRUE),
    prop_down = mean(diff(Percent.Resistant.Isolates) < 0, na.rm = TRUE),
    prop_up3y = mean(diff(Percent.Resistant.Isolates, lag = 3) > 0, na.rm = TRUE),
    prop_down3y = mean(diff(Percent.Resistant.Isolates, lag = 3) < 0, na.rm = TRUE),
    prop_up5y = mean(diff(Percent.Resistant.Isolates, lag = 5) > 0, na.rm = TRUE),
    prop_down5y = mean(diff(Percent.Resistant.Isolates, lag = 5) < 0, na.rm = TRUE)
  ) %>%
  ungroup()

library(ggplot2)

plotting_data_long <- plotting_data %>%
  pivot_longer(cols = c(prop_up, prop_down, prop_up3y, prop_down3y, prop_up5y, prop_down5y), names_to = "direction", values_to = "proportion")

data_up_minus_down_by_lag <- plotting_data %>%
  mutate(
    prop_up_minus_down = (prop_up - prop_down)*100,
    prop_up_minus_down3y = (prop_up3y - prop_down3y)*100,
    prop_up_minus_down5y = (prop_up5y - prop_down5y)*100
  ) %>%
  dplyr::select(ATC.Class, Pathogen, prop_up_minus_down, prop_up_minus_down3y, prop_up_minus_down5y) %>%
  pivot_longer(cols = c(prop_up_minus_down, prop_up_minus_down3y, prop_up_minus_down5y), names_to = "lag", values_to = "proportion")

# label the lag column to be "1 year lag", "3 year lag", "5 years"
data_up_minus_down_by_lag <- data_up_minus_down_by_lag %>%
  mutate(lag = recode(lag,
                      prop_up_minus_down = "1 year",
                      prop_up_minus_down3y = "3 years",
                      prop_up_minus_down5y = "5 years"))
# label 'proportion' as "Difference of proportions", add percentage signs to colorbar
plot_up_minus_down <- ggplot(data_up_minus_down_by_lag, aes(x = ATC.Class, y = Pathogen, fill = proportion)) +
  labs(fill = "Difference in\nproportions") +
  geom_tile() +
  facet_wrap(~lag) +
  scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, labels = scales::percent_format(scale = 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


png("heatmap_prop_up_minus_down.png", width = 6.5, height = 3, units = "in", res = 300)
print(plot_up_minus_down)
dev.off()

# library(mice)
# library(dplyr)
# source("utils.r")

# # drop columns beginning with "J"
# data <- data[, !grepl("^J", names(data))]

# # for any country-year pairs not currently in the dataset, make a row with NA values for all variables except Country and Year
# # don't include nans for pahtogen or antibiotic class
# ISO3 <- unique(data$ISO3)
# Year <- unique(data$Year)
# Pathogen <- unique(data$Pathogen)
# Pathogen <- Pathogen[!is.na(Pathogen)]
# ATC.Class <- unique(data$ATC.Class)
# ATC.Class <- ATC.Class[!is.na(ATC.Class)]
# all_country_years <- expand.grid(ISO3 = ISO3, Year = Year, Pathogen = Pathogen, ATC.Class = ATC.Class)
# data <- merge(all_country_years, data, by = c("ISO3", "Year", "Pathogen", "ATC.Class"), all.x = TRUE)

# consumption_path = "antibiotic_consumption_by_ATC3.csv"
# # Load GRAM consumption data
# consumption <- read.csv(consumption_path)
# # Extract the first 4 characters of ATC.level.3.class to get the ATC class
# consumption$ATC.level.3.class <- substr(consumption$ATC.level.3.class, 1, 4)
# # rename columns ("Location","Year","ATC level 3 class","Antibiotic consumption (DDD/1,000/day)") in consumption data and map country names onto ISO3
# consumption <- consumption %>%
# rename(ISO3 = Location, Year = Year, ATC.Class = ATC.level.3.class, Antibiotic.Consumption = Antibiotic.consumption..DDD.1.000.day.) %>%
# mutate(ISO3 = iso3_ihme_mapping$iso3[match(ISO3, iso3_ihme_mapping$country_name)])
# consumption <- consumption %>%
# group_by(ISO3, Year, ATC.Class) %>%
# summarise(Antibiotic.Consumption = sum(Antibiotic.Consumption, na.rm = TRUE),
#             .groups = "drop")

# # for any country-year pairs without consumption data yet, add the values from the consumption data to the main dataset. Don't duplicate the column
# data["Antibiotic.Consumption"] <- consumption$Antibiotic.Consumption[match(paste(data$ISO3, data$Year, data$ATC.Class, sep = "|"), paste(consumption$ISO3, consumption$Year, consumption$ATC.Class, sep = "|"))]
  


# pca_path = "Chungman/pcato10.csv"
# # Load PCA covariates and merge using vectorized join
# df.pc <- read.csv(pca_path)
# # Create composite key for PCA: ISO3|Year
# data$key_pca <- paste(data$ISO3, data$Year, sep = "|")
# df.pc$key_pca <- paste(df.pc$ISO3, df.pc$Year, sep = "|")
# idx_pca <- match(data$key_pca, df.pc$key_pca)
# data$PC1 <- df.pc$PC1[idx_pca]
# data$PC2 <- df.pc$PC2[idx_pca]
# data$PC3 <- df.pc$PC3[idx_pca]
# data$PC4 <- df.pc$PC4[idx_pca]
# data$PC5 <- df.pc$PC5[idx_pca]
# data$PC6 <- df.pc$PC6[idx_pca]
# data$PC7 <- df.pc$PC7[idx_pca]
# data$PC8 <- df.pc$PC8[idx_pca]
# data$PC9 <- df.pc$PC9[idx_pca]
# data$PC10 <- df.pc$PC10[idx_pca]
# data$GDP <- df.pc$GDP[idx_pca]
# data$key_pca <- NULL  # Clean up temporary key
# df.pc$key_pca <- NULL

# # save filled-in data to csv
# write.csv(data, "filled_in_data.csv", row.names = FALSE)

# dataMI5 = mice(data, m=20, seed=260731, printFlag = TRUE, maxit = 5)
# # 20 imputed datasets, 5 iterations

# # save the imputed datasets to CSV files
# for (i in 1:20) {
#   imputed_data <- complete(dataMI5, i)
#   write.csv(imputed_data, paste0("imputed_data_", i, ".csv"), row.names = FALSE)
# }