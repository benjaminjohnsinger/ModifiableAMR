# test_glmnet.R

# 1. Load your existing environment and functions safely
source("config.R")
source("utils.R")
source("regression_models.r") # (Ensure the bottom loop is commented out here!)
library(glmnet)
library(dplyr)

# 2. Force 'smoke_mode' ON for a fast test run
Sys.setenv(AMR_DEV_SMOKE = "false") 
runtime_options <- get_runtime_options()
set.seed(runtime_options$random_seed)

# 3. Manually load and prepare the data
message("Loading and preparing data...")
inputs <- load_model_inputs(
    merged_data_path = "merged_data_new.csv",
    merged_sums_path = "merged_data_sums_new.csv"
)

data <- inputs$data
global_consumption <- build_global_consumption_reference()
data <- scale_and_log_transform(data, global_consumption)
data_test <- limit_for_smoke_mode(data, runtime_options)

# 4. Run ONLY the glmnet model
message("Running GLMNET model...")
fit_combined_pathogen_drug_glmnet(
    data_ = data_test,
    output_tag = "glmnet_nopfac",
    runtime_options = runtime_options,
    output_prefix = "database"
)

message("Test complete. Check the Outputs/ directory for glmnet_test files.")