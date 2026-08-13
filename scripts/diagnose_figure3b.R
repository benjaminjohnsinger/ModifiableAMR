source("utils.R")
ihme_raw <- read.csv("IHME_AMR/IHME_AMR_fitted_gammas_v2.csv")

# Apply bacteria mapping (same as production)
for (i in seq_len(nrow(ihme_raw))) {
  match_val <- bacteria_mapping[bacteria_mapping$in_names == ihme_raw$pathogen[i], "canonical_names"]
  if (length(match_val) > 0 && !is.na(match_val)) ihme_raw$pathogen[i] <- match_val
}

upper_ids <- c(4, 31, 64, 103, 137, 158, 166)
IHME_upper <- ihme_raw[ihme_raw$location_id %in% upper_ids, ]

cat("All pathogens in IHME_upper (canonical names):\n")
print(sort(unique(IHME_upper$pathogen)))

# Check S. agalactiae
cat("\nAny row with 'agalactiae' in pathogen name?\n")
print(unique(ihme_raw$pathogen[grepl("agalactiae|Streptococcus ag|S\\. ag", ihme_raw$pathogen, ignore.case=TRUE)]))
cat("Any row with 'group B' or 'GBS'?\n")
print(unique(ihme_raw$pathogen[grepl("group B|GBS|faecalis|faecium", ihme_raw$pathogen, ignore.case=TRUE)]))

# Show "Overall" proportions for J01M and J01X (E. faecium's antibiotics in IHME)
opt <- read.csv("Outputs/optimistic_proportionate_consumption_by_pathogen_location_antibiotic.csv")
pes <- read.csv("Outputs/pessimistic_proportionate_consumption_by_pathogen_location_antibiotic.csv")

cat("\n--- Overall opt proportions for J01M by region ---\n")
print(opt[opt$Pathogen == "Overall" & opt$Antibiotic == "J01M", c("Location","ProportionateConsumption")])
cat("\n--- Overall pes proportions for J01M by region ---\n")
print(pes[pes$Pathogen == "Overall" & pes$Antibiotic == "J01M", c("Location","ProportionateConsumption")])

cat("\n--- Overall opt proportions for J01X by region ---\n")
print(opt[opt$Pathogen == "Overall" & opt$Antibiotic == "J01X", c("Location","ProportionateConsumption")])
cat("\n--- Overall pes proportions for J01X by region ---\n")
print(pes[pes$Pathogen == "Overall" & pes$Antibiotic == "J01X", c("Location","ProportionateConsumption")])

# Show class_gradients ranking
cat("\nClass gradients (determines order in 'Overall' scenario):\n")
gradients_df <- read.csv("Outputs/database_gradients_ATC3_PCA_canonical_weighted_all.csv")
print(gradients_df[order(gradients_df[,2], decreasing=TRUE),])

# Show E. faecium rows in IHME_upper with their true_val_att and antibiotic_class
for(atc_code in names(atc_mapping)) {
  IHME_upper[IHME_upper$antibiotic_class %in% atc_mapping[[atc_code]], "antibiotic_class"] <- atc_code
}
cat("\nE. faecium rows in IHME_upper:\n")
ef <- IHME_upper[IHME_upper$pathogen == "E. faecium", c("location_name","antibiotic_class","true_val_att","shape_all","scale_all","a_frac","b_frac")]
print(ef[ef$antibiotic_class != "Other",])
