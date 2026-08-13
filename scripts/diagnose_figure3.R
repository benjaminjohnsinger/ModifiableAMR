source("utils.R")
ihme <- read.csv("IHME_AMR/IHME_AMR_fitted_gammas_v2.csv")
for(atc_code in names(atc_mapping)) {
  ihme[ihme$antibiotic_class %in% atc_mapping[[atc_code]], "antibiotic_class"] <- atc_code
}
lower_names <- unique(iso3_ihme_mapping$lower_ihme_region)
upper_ids <- c(4, 31, 64, 103, 137, 158, 166)
IHME_upper <- ihme[ihme$location_id %in% upper_ids, ]
IHME_lower  <- ihme[ihme$location_name %in% lower_names, ]

cat("IHME_upper rows:", nrow(IHME_upper), "| unique regions:", length(unique(IHME_upper$location_name)), "\n")
cat("IHME_lower rows:", nrow(IHME_lower),  "| unique regions:", length(unique(IHME_lower$location_name)), "\n")

ef_u <- IHME_upper[IHME_upper$pathogen == "Enterococcus faecium", ]
cat("\nE. faecium true_val_att by region (UPPER):\n")
print(tapply(ef_u$true_val_att, ef_u$location_name, sum, na.rm=TRUE))

ef_l <- IHME_lower[IHME_lower$pathogen == "Enterococcus faecium", ]
cat("\nE. faecium true_val_att by region (LOWER) — summed to upper:\n")
lower_to_upper <- unique(iso3_ihme_mapping[, c("lower_ihme_region", "ihme_region")])
ef_l$upper_region <- lower_to_upper$ihme_region[match(ef_l$location_name, lower_to_upper$lower_ihme_region)]
print(tapply(ef_l$true_val_att, ef_l$upper_region, sum, na.rm=TRUE))

cat("\nE. faecium total: UPPER =", sum(ef_u$true_val_att, na.rm=TRUE),
    " LOWER =", sum(ef_l$true_val_att, na.rm=TRUE), "\n")

# S. agalactiae
sa_u <- IHME_upper[IHME_upper$pathogen == "Streptococcus agalactiae", ]
sa_l <- IHME_lower[IHME_lower$pathogen  == "Streptococcus agalactiae", ]
cat("\nS. agalactiae total: UPPER =", sum(sa_u$true_val_att, na.rm=TRUE),
    " LOWER =", sum(sa_l$true_val_att, na.rm=TRUE), "\n")

# Check response values for E. faecium antibiotics
results <- read.csv(getOption("amr_burden_results_path",
    "Outputs/database_gradients_pathogen_ATC3_PCA_canonical_weighted_main.csv"))
cat("\nE. faecium response values from regression:\n")
print(results[results$Pathogen == "E. faecium", c("Antibiotic","Response")])

# Check what antibiotics appear in E. faecium IHME rows
cat("\nE. faecium antibiotics in IHME_upper:\n")
print(unique(ef_u$antibiotic_class))
cat("\nE. faecium antibiotics in IHME_lower:\n")
print(unique(ef_l$antibiotic_class))

# Check the optimistic proportions for E. faecium
opt <- read.csv("Outputs/optimistic_proportionate_consumption_by_pathogen_location_antibiotic.csv")
ef_opt <- opt[opt$Pathogen == "E. faecium", ]
cat("\nE. faecium optimistic proportions (non-1.0):\n")
print(ef_opt[ef_opt$ProportionateConsumption != 1.0, ])
