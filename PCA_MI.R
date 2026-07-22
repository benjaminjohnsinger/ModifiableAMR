# random seed
set.seed(260602)

requiredPackages = c("dplyr", "tidyr", "data.table", "readr","tibble","MASS", "missMDA")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

df.1 <- fread("Chungman/UNIndex/1to2.csv"  , select = c("Series Name","Country Name", "2001 [YR2001]", "2002 [YR2002]"), col.names = c("Series", "Country", "2001", "2002"))
df.2 <- fread("Chungman/UNIndex/3to5.csv"  , select = c("Series Name","Country Name", "2003 [YR2003]", "2004 [YR2004]", "2005 [YR2005]"), col.names = c("Series", "Country", "2003", "2004", "2005"))
df.3 <- fread("Chungman/UNIndex/6to8.csv"  , select = c("Series Name","Country Name", "2006 [YR2006]", "2007 [YR2007]", "2008 [YR2008]"), col.names = c("Series", "Country", "2006", "2007", "2008"))
df.4 <- fread("Chungman/UNIndex/9to11.csv" , select = c("Series Name","Country Name", "2009 [YR2009]", "2010 [YR2010]", "2011 [YR2011]"), col.names = c("Series", "Country", "2009", "2010", "2011"))
df.5 <- fread("Chungman/UNIndex/12to13.csv", select = c("Series Name","Country Name", "2012 [YR2012]", "2013 [YR2013]"), col.names = c("Series", "Country", "2012", "2013"))
df.6 <- fread("Chungman/UNIndex/14to16.csv", select = c("Series Name","Country Name", "2014 [YR2014]", "2015 [YR2015]", "2016 [YR2016]"), col.names = c("Series", "Country", "2014", "2015", "2016"))
df.7 <- fread("Chungman/UNIndex/17to19.csv", select = c("Series Name","Country Name", "2017 [YR2017]", "2018 [YR2018]", "2019 [YR2019]"), col.names = c("Series", "Country", "2017", "2018", "2019"))
df.8 <- fread("Chungman/UNIndex/20to22.csv", select = c("Series Name","Country Name", "2020 [YR2020]", "2021 [YR2021]", "2022 [YR2022]"), col.names = c("Series", "Country", "2020", "2021", "2022"))

df_list <- list(df.1, df.2, df.3, df.4, df.5, df.6, df.7, df.8)
df_list <- lapply(df_list, function(df) df[1:(.N - 5)])

first_two_cols <- lapply(df_list, function(df) df[, .(Series, Country)])
sapply(2:8, function(i) all.equal(first_two_cols[[1]], first_two_cols[[i]])) # check whether the Country series are same.

dfs.years <- lapply(df_list, function(df) df[, -c("Series", "Country"), with = FALSE])
df.combined <- cbind(df_list[[1]][,1:2], do.call(cbind, dfs.years))
dim(df.combined); colnames(df.combined)
df.long <- df.combined %>% pivot_longer(cols = -c(Series, Country), names_to = "Year", values_to = "Value")
df.long <- df.long %>% filter(Series != "")
df.ready <- df.long %>% pivot_wider(names_from = Series, values_from = Value)
df.ready <- df.ready[, -((ncol(df.ready)-1):ncol(df.ready))] #(5874,330)
dim(df.ready) # (5874,328)
unique(df.ready$Country)

non_countries <- c("Africa Eastern and Southern", "Africa Western and Central", "Arab World",
                   "Caribbean small states", "Central Europe and the Baltics", "East Asia & Pacific",
                   "East Asia & Pacific (excluding high income)", "East Asia & Pacific (IDA & IBRD countries)",
                   "Early-demographic dividend", "Euro area", "Europe & Central Asia",
                   "Europe & Central Asia (excluding high income)", "Europe & Central Asia (IDA & IBRD countries)",
                   "European Union", "Fragile and conflict affected situations",
                   "Heavily indebted poor countries (HIPC)", "High income", "IBRD only",
                   "IDA & IBRD total", "IDA blend", "IDA only", "IDA total",
                   "Late-demographic dividend", "Latin America & Caribbean",
                   "Latin America & Caribbean (excluding high income)",
                   "Latin America & the Caribbean (IDA & IBRD countries)",
                   "Least developed countries: UN classification", "Low & middle income", "Low income",
                   "Lower middle income", "Middle East & North Africa",
                   "Middle East & North Africa (excluding high income)",
                   "Middle East & North Africa (IDA & IBRD countries)", "Middle income",
                   "North America", "Not classified", "OECD members", "Other small states",
                   "Pacific island small states", "Post-demographic dividend", "Pre-demographic dividend",
                   "Small states", "South Asia", "South Asia (IDA & IBRD)", "Sub-Saharan Africa",
                   "Sub-Saharan Africa (excluding high income)", "Sub-Saharan Africa (IDA & IBRD countries)",
                   "Upper middle income", "World")

# Example: Filter from your dataframe
df.ready <- df.ready %>% filter(!(Country %in% non_countries))
print(colnames(df.ready))
#GDP from world bank (2021 international $): Yemen and Venezuela from IMF data (GDP per capita current international dollar -> should be deflated)
GDP.PPP <- read.csv("Chungman/GDP_PPP.csv")
colnames(GDP.PPP)[-1] <- gsub("^X", "", colnames(GDP.PPP)[-1])
GDP.PPP <- GDP.PPP %>% pivot_longer(cols = -Country, names_to = "Year", values_to = "GDP")
GDP.PPP$Year <- as.numeric(GDP.PPP$Year)
a <- GDP.PPP$Country
b <- df.ready$Country
setdiff(b,a)
df.ready$Year <- as.integer(df.ready$Year)
df.ready.gdp <- df.ready %>% inner_join(GDP.PPP, by = c("Country" = "Country", "Year" = "Year"))
unique(df.ready.gdp[!complete.cases(df.ready.gdp), ]$Country)
df.ready.gdp <- df.ready.gdp[!is.na(df.ready.gdp$GDP),]
# df.ready.gdp <- df.ready.gdp[df.ready.gdp$GDP < 20500, ]

meta_cols <- c("Country", "Year", "GDP")
df.pca <- df.ready.gdp %>% mutate(across(-all_of(meta_cols), ~ na_if(., "..")   )) %>% # convert ".." to NA
  mutate(across(-all_of(meta_cols), ~ as.numeric(.), .names = "{.col}")) # convert character to numeric

hist(round(colSums(is.na(df.pca))/nrow(df.pca), digits = 3),)
na_prop <- colSums(is.na(df.pca)) / nrow(df.pca)
df.pca.ready <- df.pca[, na_prop <= 0.3] 
dim(df.pca.ready) #   <0.3:(5874,169)
df.numeric <- df.pca.ready[, !(names(df.pca.ready) %in% c("Country", "Year"))]

mi.pca <- MIPCA(df.numeric, ncp = 2, nboot = 1)  # 5 imputed datasets
# Access the imputed datasets
imputed.list <- mi.pca$res.MI  # This is a list of imputed datasets
# Run PCA on each imputed dataset
pca.results <- lapply(imputed.list, prcomp, scale. = TRUE)
# screeplot for first imputed set
screeplot(pca.results[[1]], main = "Screeplot of First Imputed Dataset", type = "lines")
summary(pca.results[[1]]); pca.results[[1]]$x

pc.scores <- pca.results[[1]]$x[, 1:3]  # Get PC1 to PC3
colnames(pc.scores) <- c("PC1", "PC2", "PC3")
df.pc <- bind_cols(df.pca.ready[, c("Country", "Year","GDP")], as.tibble(pc.scores))

# what proportion of variance is explained by PC1, PC2, and PC3?
variance_explained <- summary(pca.results[[1]])$importance[2, 1:3]  # Proportion of variance explained by PC1, PC2, and PC3
print(sum(variance_explained))

# Create a dataframe with PC, variance explained, and top 3 variables
# with loadings
variance_explained <- summary(pca.results[[1]])$importance[2, 1:10]
cum_variance_explained <- 
  cumsum(summary(pca.results[[1]])$importance[2, 1:10])
top_vars_data <- lapply(1:10, function(i) {
  pc_loadings <- pca.results[[1]]$rotation[, i]
  top_indices <- order(abs(pc_loadings), decreasing = TRUE)[1:3]
  top_var_names <- names(pc_loadings)[top_indices]
  top_var_loadings <- pc_loadings[top_indices]

  data.frame(
    PC = paste0("PC", i),
    Variance_Explained = round(variance_explained[i], 4),
    Cumulative_Variance = round(cum_variance_explained[i], 4),
    Var1 = paste0(top_var_names[1], " (", 
                  round(top_var_loadings[1], 3), ")"),
    Var2 = paste0(top_var_names[2], " (", 
                  round(top_var_loadings[2], 3), ")"),
    Var3 = paste0(top_var_names[3], " (", 
                  round(top_var_loadings[3], 3), ")")
  )
})
top_vars_df <- do.call(rbind, top_vars_data)
top_vars_df[] <- lapply(top_vars_df, gsub, pattern = '"', replacement = "")
write.table(top_vars_df, file = "Outputs/pca_tsv.txt", 
            row.names = FALSE, sep = "\t")

# str(df.pc)
# par(mfrow=c(2,1))
# hist(df.pc$PC1); hist(df.pc$PC2)#; hist(df.pc$PC3)
# par(mfrow=c(1,1))
# # scatter plots of PC1 vs PC2, PC1 vs PC3, and PC2 vs PC3, in panels
# library(ggplot2)
# df.pc <- df.pc[(df.pc$PC3 < 2) & (df.pc$PC3 > 0), ]
# ggplot(df.pc, aes(x = PC1, y = PC2, color = Country)) +
#   geom_point() +
#   labs(title = "PCA: PC1 vs PC2", x = "PC1", y = "PC2")
#   # scale_color_gradient(low = "blue", high = "red") +
#   # theme_minimal()
# ggsave("Chungman/pc1_vs_pc2_Country_colour_big_subset.png", width = 8, height = 6)
# ggplot(df.pc, aes(x = PC1, y = PC3, color = Country)) +
#   geom_point() +
#   labs(title = "PCA: PC1 vs PC3", x = "PC1", y = "PC3")
#   # scale_color_gradient(low = "blue", high = "red") +
#   # theme_minimal()
# ggsave("Chungman/pc1_vs_pc3_Country_colour_big_subset.png", width = 8, height = 6)
# ggplot(df.pc, aes(x = PC2, y = PC3, color = Country)) +
#   geom_point() +
#   labs(title = "PCA: PC2 vs PC3", x = "PC2", y = "PC3")
#   # scale_color_gradient(low = "blue", high = "red") +
#   # theme_minimal()
# ggsave("Chungman/pc2_vs_pc3_Country_colour_big_subset.png", width = 8, height = 6)

# saveRDS(df.pc, "Chungman/UNIndex/pc_gdp.RDS")
# write.csv(df.pc, "Chungman/pc3_gdp.csv", row.names = FALSE)


