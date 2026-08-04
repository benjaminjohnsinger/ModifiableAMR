# ## Linear regression on database DDD data
# ## BJS March 2025
library(tidyverse)
library(dplyr)
library(glmnet)
library(skedastic)
library(MuMIn)
library(future)
library(future.apply)
source("utils.R")
source("config.R")

# NOTE: Lagged-response modeling is preserved as exploratory functionality only.
# Main and supplementary publication analyses should run with apply_lagged_response = FALSE.

is_truthy_env <- function(name, default = FALSE) {
    value <- Sys.getenv(name, unset = if (default) "true" else "")
    if (identical(value, "")) {
        return(default)
    }
    tolower(value) %in% c("1", "true", "t", "yes", "y", "on")
}

get_integer_env <- function(name, default) {
    value <- Sys.getenv(name, unset = as.character(default))
    suppressWarnings(parsed <- as.integer(value))
    if (is.na(parsed)) {
        return(default)
    }
    parsed
}

resolve_reference_year <- function(scenario = NULL, job = NULL) {
    env_year <- Sys.getenv("REF_YEAR", unset = "")
    if (nzchar(env_year)) {
        return(as.character(get_integer_env("REF_YEAR", default = 2018)))
    }

    if (!is.null(job) && !is.null(job$reference_year)) {
        return(as.character(job$reference_year))
    }

    if (!is.null(scenario) && identical(scenario, "main_2000")) {
        return("2000")
    }

    "2018"
}

log_info <- function(..., verbose = TRUE) {
    if (isTRUE(verbose)) {
        message(...)
    }
}

get_runtime_options <- function() {
    # Smoke mode is developer-facing only; it is intentionally not part of normal user docs.
    smoke_mode <- is_truthy_env("AMR_DEV_SMOKE", default = FALSE) ||
        is_truthy_env("AMR_SMOKE", default = FALSE)
    verbose <- is_truthy_env("AMR_VERBOSE", default = TRUE)
    list(
        smoke_mode = smoke_mode,
        verbose = verbose,
        random_seed = get_integer_env("AMR_RANDOM_SEED", default = 20260506),
        boot_nsim = if (smoke_mode) 0 else 1000,
        smoke_max_classes = if (smoke_mode) 3 else Inf,
        smoke_max_pathogens = if (smoke_mode) 3 else Inf,
        smoke_max_pairs = if (smoke_mode) 9 else Inf,
        smoke_max_rows_per_pair = if (smoke_mode) 40 else Inf
    )
}

pool_rubins_rules <- function(estimates, std_errors) {
    # Remove any NAs caused by complete model failures
    valid_idx <- !is.na(estimates) & !is.na(std_errors)
    estimates <- estimates[valid_idx]
    std_errors <- std_errors[valid_idx]
    m <- length(estimates)
    
    if (m == 0) {
        return(list(gradient = NA_real_, lower_ci = NA_real_, upper_ci = NA_real_, se = NA_real_))
    }
    
    # Pooled point estimate
    pooled_est <- mean(estimates)
    
    # Variances
    vw <- mean(std_errors^2)
    vb <- var(estimates)
    
    # Handle single successful model edge case
    if (is.na(vb)) vb <- 0 
    
    vt <- vw + (1 + (1/m)) * vb
    pooled_se <- sqrt(vt)
    
    # Calculate 95% CIs
    lower_ci <- pooled_est - (1.96 * pooled_se)
    upper_ci <- pooled_est + (1.96 * pooled_se)
    
    list(
        gradient = pooled_est,
        lower_ci = lower_ci,
        upper_ci = upper_ci,
        se = pooled_se
    )
}

limit_for_smoke_mode <- function(df, runtime_options) {
    if (!isTRUE(runtime_options$smoke_mode)) {
        return(df)
    }

    # Keep a bounded number of pathogen-antibiotic pairs and rows per pair for fast smoke checks.
    pair_counts <- df %>%
        count(Pathogen, Antibiotic, name = "n") %>%
        arrange(desc(n), Pathogen, Antibiotic)
    keep_pairs <- head(pair_counts, runtime_options$smoke_max_pairs)

    if (nrow(keep_pairs) == 0) {
        return(df[0, , drop = FALSE])
    }

    df %>%
        inner_join(keep_pairs %>% select(Pathogen, Antibiotic), by = c("Pathogen", "Antibiotic")) %>%
        group_by(Pathogen, Antibiotic) %>%
        slice_head(n = runtime_options$smoke_max_rows_per_pair) %>%
        ungroup()
}

select_income_slice <- function(inputs, income) {
    if (income == "HIC") {
        return(inputs$data_HIC)
    }
    if (income == "LMIC") {
        return(inputs$data_LMIC)
    }
    inputs$data
}

map_nagorsen_to_atc_class <- function(class_vector) {
    mapped <- class_vector
    for (atc_code in names(atc_mapping)) {
        mapped[mapped %in% atc_mapping[[atc_code]]] <- atc_code
    }
    mapped
}

load_nagorsen_model_inputs <- function(
    nagorsen_path = "Nagorsen_clean.csv",
    pca_path = "Chungman/Chungman_pca_renamed.csv",
    prepared_data_path = "summed_data_Nagorsen_hospital_to_all_filtered.csv",
    min_entries_per_combo = 20
) {
    runtime_options <- get_runtime_options()

    if (file.exists(prepared_data_path)) {
        data <- read.csv(prepared_data_path)
        data <- data[complete.cases(data), ]
        pathogen_drug_counts <- table(paste(data$Pathogen, data$Antibiotic))
        pathogen_drug_to_remove <- names(pathogen_drug_counts[pathogen_drug_counts <= min_entries_per_combo])
        data <- data[!paste(data$Pathogen, data$Antibiotic) %in% pathogen_drug_to_remove, ]

        log_info("[ddd-linear-model] Loaded prepared Nagorsen data: ", prepared_data_path, verbose = runtime_options$verbose)
        log_info("[ddd-linear-model] Nagorsen rows retained: ", nrow(data), verbose = runtime_options$verbose)
        log_info("[ddd-linear-model] Nagorsen pathogens retained: ", length(unique(data$Pathogen)), verbose = runtime_options$verbose)
        log_info("[ddd-linear-model] Nagorsen ATC classes retained: ", length(unique(data$Antibiotic)), verbose = runtime_options$verbose)

        return(list(data = data, data_HIC = data, data_LMIC = data))
    }

    data <- read.csv(
        nagorsen_path,
        colClasses = c("units" = "character"),
        na.strings = c("NA")
    )

    data <- data[
        !is.na(data$amt_consumed) &
            !is.na(data$units) &
            !is.na(data$class_for_resistance) &
            !is.na(data$pathogen),
    ]
    data <- data[data$amt_consumed < 10000, ]

    data$pathogen <- vapply(data$pathogen, get_bacteria_name, character(1))
    data$class_for_resistance <- map_nagorsen_to_atc_class(data$class_for_resistance)
    data <- data[!data$class_for_resistance %in% c("J01X", "Other"), ]

    # Align unit conventions to DDD/1000/day before filtering.
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

    log_info("[ddd-linear-model] Nagorsen rows retained: ", nrow(data), verbose = runtime_options$verbose)
    log_info("[ddd-linear-model] Nagorsen pathogens retained: ", length(unique(data$Pathogen)), verbose = runtime_options$verbose)
    log_info("[ddd-linear-model] Nagorsen ATC classes retained: ", length(unique(data$Antibiotic)), verbose = runtime_options$verbose)

    list(data = data, data_HIC = data, data_LMIC = data)
}

job_or_default <- function(job, key, default_value) {
    if (!is.null(job[[key]])) {
        return(job[[key]])
    }
    default_value
}
load_model_inputs <- function(
    summed_data_path = "summed_data_new.csv",
    merged_sums_path = "summed_data_sums_new.csv",
    antibiotic_col = "ATC.Class",
    consumption_class_col = "ATC.Class",
    min_entries_per_combo = 10,
    min_isolates_per_combo = 100
) {
    data <- read.csv(summed_data_path)

    if (!antibiotic_col %in% names(data)) {
        stop(
            "[ddd-linear-model] antibiotic_col not found in input data: ",
            antibiotic_col,
            call. = FALSE
        )
    }
    if (!consumption_class_col %in% names(data)) {
        stop(
            "[ddd-linear-model] consumption_class_col not found in input data: ",
            consumption_class_col,
            call. = FALSE
        )
    }

    runtime_options <- get_runtime_options()
    log_info("[ddd-linear-model] Input file: ", summed_data_path, verbose = runtime_options$verbose)
    log_info("[ddd-linear-model] Unique ISO3 count: ", length(table(data$ISO3)), verbose = runtime_options$verbose)
    log_info("[ddd-linear-model] Max rows for one ISO3: ", max(table(data$ISO3)), verbose = runtime_options$verbose)
    log_info("[ddd-linear-model] Years included: ", paste(sort(unique(data$Year)), collapse = ", "), verbose = runtime_options$verbose)

    data <- data %>%
        rename(
            Consumption = Antibiotic.Consumption,
            Resistance = Percent.Resistant.Isolates,
            Weight = Total.Isolates
        ) %>%
        mutate(
            Antibiotic = .data[[antibiotic_col]],
            Consumption.Class = .data[[consumption_class_col]]
        )

    country_covariates <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "GDP", "Year")
    data$GDP <- data$GDP / mean(data$GDP, na.rm = TRUE)

    data <- data[
        !is.na(data$Consumption) & !is.na(data$Resistance) &
            !is.na(data$Pathogen) & !is.na(data$Antibiotic) & !is.na(data$Weight),
    ]

    pathogen_drug_counts <- table(paste(data$Pathogen, data$Antibiotic))
    pathogen_drug_to_remove <- names(pathogen_drug_counts[pathogen_drug_counts <= min_entries_per_combo])
    data <- data[!paste(data$Pathogen, data$Antibiotic) %in% pathogen_drug_to_remove, ]

    summed_data_sums <- read.csv(merged_sums_path)

    if (!antibiotic_col %in% names(summed_data_sums)) {
        stop(
            "[ddd-linear-model] antibiotic_col not found in merged sums data: ",
            antibiotic_col,
            call. = FALSE
        )
    }

    summed_data_sums <- summed_data_sums %>%
        mutate(Antibiotic = .data[[antibiotic_col]])

    pathogen_drug_counts <- summed_data_sums %>%
        group_by(Pathogen, Antibiotic) %>%
        summarise(Total.Isolates = sum(Total.Isolates, na.rm = TRUE), .groups = "drop") %>%
        filter(Total.Isolates < min_isolates_per_combo) %>%
        select(Pathogen, Antibiotic)

    pathogen_drug_to_remove <- paste(pathogen_drug_counts$Pathogen, pathogen_drug_counts$Antibiotic)
    data <- data[!paste(data$Pathogen, data$Antibiotic) %in% pathogen_drug_to_remove, ]

    pathogen_drug_counts <- data %>%
        group_by(Pathogen, Antibiotic) %>%
        summarise(Resistance = sum(Resistance, na.rm = TRUE), .groups = "drop") %>%
        filter(Resistance == 0) %>%
        select(Pathogen, Antibiotic)

    pathogen_drug_to_remove <- paste(pathogen_drug_counts$Pathogen, pathogen_drug_counts$Antibiotic)
    data <- data[!paste(data$Pathogen, data$Antibiotic) %in% pathogen_drug_to_remove, ]

    log_info("[ddd-linear-model] Remaining rows after filtering: ", nrow(data), verbose = runtime_options$verbose)
    log_info("[ddd-linear-model] Antibiotic classes retained: ", length(unique(data$Antibiotic)), verbose = runtime_options$verbose)
    log_info("[ddd-linear-model] Pathogens retained: ", length(unique(data$Pathogen)), verbose = runtime_options$verbose)

    # countries_per_pathogen <- data %>%
    #     group_by(Pathogen) %>%
    #     summarise(Countries = paste(unique(ISO3), collapse = ", "), .groups = "drop")

    # for (i in seq_len(nrow(countries_per_pathogen))) {
    #     log_info(paste0(countries_per_pathogen$Pathogen[i], ": ", countries_per_pathogen$Countries[i]), verbose = runtime_options$verbose)
    # }

    data$lending_group <- iso3_ihme_mapping$lending_group[match(data$ISO3, iso3_ihme_mapping$iso3)]

    high_income_entries <- sum(data$lending_group == "High income")
    high_income_isolates <- sum(data$Weight[data$lending_group == "High income"], na.rm = TRUE)
    log_info(paste("HIC entries:", high_income_entries), verbose = runtime_options$verbose)
    log_info(paste("HIC isolates:", high_income_isolates), verbose = runtime_options$verbose)

    lmics_entries <- sum(data$lending_group != "High income")
    lmics_isolates <- sum(data$Weight[data$lending_group != "High income"], na.rm = TRUE)
    log_info(paste("LMIC entries:", lmics_entries), verbose = runtime_options$verbose)
    log_info(paste("LMIC isolates:", lmics_isolates), verbose = runtime_options$verbose)

    data_LMIC <- data[data$lending_group != "High income", ]
    data_HIC <- data[data$lending_group == "High income", ]

    data <- data %>%
        select(Consumption, Resistance, Pathogen, Antibiotic, Consumption.Class, Weight, ISO3, 
               all_of(country_covariates), ends_with(".Consumption"))

    list(data = data, data_HIC = data_HIC, data_LMIC = data_LMIC)
}
build_global_consumption_reference <- function(consumption_path = "antibiotic_consumption_by_ATC3.csv", year = "2018") {
    consumption <- read.csv(consumption_path)
    global_consumption <- consumption[consumption$Location == "Global", ]
    global_consumption <- global_consumption[global_consumption$Year == year, ]
    global_consumption <- select(global_consumption, ATC.level.3.class, Antibiotic.consumption..DDD.1.000.day.)
    global_consumption <- global_consumption %>%
        rename(
            Antibiotic = ATC.level.3.class,
            Global.Consumption = Antibiotic.consumption..DDD.1.000.day.
        )
    global_consumption$Antibiotic <- sub("-.*", "", global_consumption$Antibiotic)
    global_consumption
}

scale_and_log_transform <- function(df, global_consumption, model_family = "gaussian") {
    join_by <- if ("Consumption.Class" %in% names(df)) {
        c("Consumption.Class" = "Antibiotic")
    } else {
        c("Antibiotic" = "Antibiotic")
    }

    df <- df %>%
        left_join(global_consumption, by = join_by)

    if (all(is.na(df$Global.Consumption))) {
        stop(
            "[ddd-linear-model] Global consumption normalization failed: no matching ATC classes between model data and reference consumption table.",
            call. = FALSE
        )
    }

    df <- df %>%
        mutate(Consumption = Consumption / Global.Consumption) %>%
        select(-Global.Consumption)

    df <- df %>%
        filter(!is.na(Consumption) & !is.na(Resistance) & !is.na(Weight))

    df$Consumption <- log(df$Consumption + 1)
    if (model_family == "binomial") {
        df$Resistance <- df$Resistance / 100
        df$resistant_count <- round(df$Resistance * df$Weight)
    } else if (model_family == "ppml") {
        # Keep Resistance as the raw percentage, do not log it
        df$Weight <- df$Weight / max(df$Weight, na.rm = TRUE)
    } else {
        # Gaussian defaults
        df$Resistance <- log(df$Resistance + 1)
        df$Weight <- df$Weight / max(df$Weight, na.rm = TRUE)
    }

    # make GDP and Year have mean 0 and sd 1
    df$GDP <- scale(df$GDP)
    df$Year <- scale(df$Year)

    df
}

get_fixed_effects_formula <- function(extra_pcs = FALSE) {
    if (extra_pcs) {
        Resistance ~ Consumption + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + GDP + Year
    } else {
        Resistance ~ Consumption + PC1 + PC2 + PC3 + GDP + Year
    }
}

fit_weighted_lm <- function(data_subset, extra_pcs = FALSE) {
    lm(
        formula = get_fixed_effects_formula(extra_pcs = extra_pcs),
        data = data_subset,
        weights = Weight
    )
}

fit_ppml_glm <- function(data_subset, extra_pcs = FALSE) {
    glm(
        formula = get_fixed_effects_formula(extra_pcs = extra_pcs),
        data = data_subset,
        family = quasipoisson(link = "log"),
        weights = Weight
    )
}

fit_binomial_glm <- function(data_subset, extra_pcs = FALSE) {
    glm(
        formula = get_fixed_effects_formula(extra_pcs = extra_pcs) |> update(cbind(resistant_count, Weight - resistant_count) ~ .),
        data = data_subset,
        family = binomial(link = "logit")
    )
}

fit_random_lmer <- function(data_subset, random_effect_var, extra_pcs = FALSE) {
    formula_str <- paste0(
        paste0(
            "Resistance ~ Consumption + (Consumption||", random_effect_var,
            ") + PC1 + PC2 + PC3",
            if (extra_pcs) " + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10" else "",
            " + GDP + Year"
        )
    )
    lmer(
        formula = as.formula(formula_str),
        data = data_subset,
        weights = Weight
    )
}

fit_binomial_glmer <- function(data_subset, random_effect_var, extra_pcs = FALSE) {
    formula_str <- paste0(
        "cbind(resistant_count, Weight - resistant_count) ~ Consumption + (Consumption||", random_effect_var,
        paste0(
            ") + PC1 + PC2 + PC3",
            if (extra_pcs) " + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10" else "",
            " + GDP + Year"
        )
    )
    glmer(
        formula = as.formula(formula_str),
        data = data_subset,
        family = binomial(link = "logit"),
    )
}

# install.packages("glmmTMB")
library(glmmTMB)

fit_random_ppml_glmer <- function(data_subset, random_effect_var, extra_pcs = FALSE) {
    formula_str <- paste0(
        "Resistance ~ Consumption + (Consumption||", random_effect_var,
        ") + PC1 + PC2 + PC3",
        if (extra_pcs) " + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10" else "",
        " + GDP + Year"
    )

    # suppressWarnings hides the "non-integer counts" message
    suppressWarnings(
        glmmTMB::glmmTMB(
            formula = as.formula(formula_str),
            data = data_subset,
            family = poisson(link = "log"), 
            weights = Weight
        )
    )
}

build_output_path <- function(prefix, output_tag) {
    paste0("Outputs/", prefix, "_", output_tag, ".csv")
}

write_random_effects_outputs <- function(
    gradient_prefix, lower_prefix, upper_prefix, bootstrap_prefix, output_tag,
    gradients, lower_ci, upper_ci, bootstraps,
    r_squareds = NULL, white_tests = NULL, rmses = NULL, aics = NULL, bics = NULL, vif_list = NULL
) {
    # If the new metrics are provided, attach them to the gradients dataframe
    if (!is.null(r_squareds)) {
        gradients_df <- data.frame(
            Gradient = gradients,
            R_squared = r_squareds,
            White_test_p_value = white_tests,
            RMSE = rmses,
            AIC = aics,
            BIC = bics
        )
        
        # Format VIFs into a dataframe if available
        if (!is.null(vif_list) && length(vif_list) > 0) {
            vif_rows <- lapply(vif_list, function(x) {
                # Keep one row per fitted model, even when VIF fails or is partially unavailable.
                if (is.null(x) || length(x) == 0 || all(is.na(x))) {
                    return(data.frame(.vif_fallback = NA_real_))
                }

                x_names <- names(x)
                x <- as.numeric(x)
                if (is.null(x_names) || any(x_names == "")) {
                    x_names <- paste0("V", seq_along(x))
                }
                names(x) <- x_names
                as.data.frame(as.list(x), check.names = FALSE)
            })

            vif_df <- dplyr::bind_rows(vif_rows)
            if (".vif_fallback" %in% colnames(vif_df) && ncol(vif_df) > 1) {
                vif_df$.vif_fallback <- NULL
            }
            colnames(vif_df) <- paste0("VIF.", colnames(vif_df))
            gradients_df <- cbind(gradients_df, vif_df)
        }
        
        write.csv(gradients_df, build_output_path(gradient_prefix, output_tag), row.names = TRUE)
    } else {
        write.csv(gradients, build_output_path(gradient_prefix, output_tag), row.names = TRUE)
    }
    
    write.csv(lower_ci, build_output_path(lower_prefix, output_tag), row.names = TRUE)
    write.csv(upper_ci, build_output_path(upper_prefix, output_tag), row.names = TRUE)
    write.csv(bootstraps, build_output_path(bootstrap_prefix, output_tag), row.names = TRUE)
}

initialize_random_effects_accumulator <- function() {
    list(
        labels = character(),
        gradients = numeric(),
        intercepts = numeric(),
        lower_ci = numeric(),
        upper_ci = numeric(),
        bootstraps = data.frame(),
        r_squareds = numeric(),
        white_tests = numeric(),
        rmses = numeric(),
        aics = numeric(),
        bics = numeric(),
        vif_list = list()
    )
}

append_random_effects_result <- function(
    accumulator, label, gradient, intercept, lower_ci, upper_ci, bootstrap_df,
    r_sq, white_p, rmse, aic, bic, vifs
) {
    accumulator$labels <- c(accumulator$labels, label)
    accumulator$gradients <- c(accumulator$gradients, gradient)
    accumulator$intercepts <- c(accumulator$intercepts, intercept)
    accumulator$lower_ci <- c(accumulator$lower_ci, lower_ci)
    accumulator$upper_ci <- c(accumulator$upper_ci, upper_ci)
    accumulator$bootstraps <- rbind(accumulator$bootstraps, bootstrap_df)
    
    # New metrics
    accumulator$r_squareds <- c(accumulator$r_squareds, r_sq)
    accumulator$white_tests <- c(accumulator$white_tests, white_p)
    accumulator$rmses <- c(accumulator$rmses, rmse)
    accumulator$aics <- c(accumulator$aics, aic)
    accumulator$bics <- c(accumulator$bics, bic)
    accumulator$vif_list[[length(accumulator$vif_list) + 1]] <- vifs
    
    accumulator
}

format_combo_results <- function(named_vector, value_col_name) {
    # Convert the named vector into a clean data frame
    df <- data.frame(
        Combo_Label = names(named_vector),
        Value = unname(named_vector)
    )
    # Rename the value column to whatever is appropriate (Gradient, Lower_CI, etc.)
    colnames(df)[2] <- value_col_name
    
    # Split the Combo_Label into Pathogen and Antibiotic columns
    df %>%
        tidyr::separate(Combo_Label, into = c("Pathogen", "Antibiotic"), sep = "___")
}

extract_lm_consumption_summary <- function(model, boot_nsim, data_subset = NULL, extra_pcs = FALSE) {
    gradient <- summary(model)$coefficients["Consumption", 1]
    intercept <- summary(model)$coefficients["(Intercept)", 1]

    if (boot_nsim > 0 && !is.null(data_subset)) {
        boot_result <- tryCatch(
            {
                unique_countries <- unique(data_subset$ISO3)
                n_countries <- length(unique_countries)
                
                boot_values <- replicate(boot_nsim, {
                    sampled_countries <- sample(unique_countries, n_countries, replace = TRUE)
                    idx <- unlist(lapply(sampled_countries, function(iso) {
                        which(data_subset$ISO3 == iso)
                    }))
                    
                    bdf <- data_subset[idx, ]
                    mb <- lm(
                        formula = get_fixed_effects_formula(extra_pcs = extra_pcs),
                        data = bdf,
                        weights = Weight
                    )
                    summary(mb)$coefficients["Consumption", 1]
                })
                
                list(
                    lower_ci = quantile(boot_values, 0.025, na.rm = TRUE),
                    upper_ci = quantile(boot_values, 0.975, na.rm = TRUE),
                    bootstrap_values = boot_values
                )
            },
            error = function(e) {
                warning("Cluster bootstrap failed for LM Consumption term; using Wald CI fallback. Error: ", conditionMessage(e))
                intervals <- confint(model)
                list(
                    lower_ci = intervals["Consumption", 1],
                    upper_ci = intervals["Consumption", 2],
                    bootstrap_values = gradient
                )
            }
        )
        lower_ci <- boot_result$lower_ci
        upper_ci <- boot_result$upper_ci
        bootstrap_values <- boot_result$bootstrap_values
    } else {
        intervals <- confint(model)
        lower_ci <- intervals["Consumption", 1]
        upper_ci <- intervals["Consumption", 2]
        bootstrap_values <- gradient
    }

    # Inside extract_lm_consumption_summary and extract_lmer_consumption_summary...
    # Calculate SE from bootstraps, fallback to Wald SE if bootstraps failed/0
    if (boot_nsim > 0 && !is.null(data_subset) && length(bootstrap_values) > 1) {
        calc_se <- sd(bootstrap_values, na.rm = TRUE)
    } else {
        # Fallback to analytic standard error
        calc_se <- summary(model)$coefficients["Consumption", 2] 
    }

    list(
        gradient = gradient,
        intercept = intercept,
        lower_ci = lower_ci,
        upper_ci = upper_ci,
        se = calc_se, # NEW: Standard error for MI
        bootstrap_values = bootstrap_values
    )
}

extract_glm_ppml_consumption_summary <- function(model, boot_nsim, data_subset = NULL, extra_pcs = FALSE) {
    gradient <- coef(model)[["Consumption"]]
    intercept <- coef(model)[["(Intercept)"]]

    if (boot_nsim > 0 && !is.null(data_subset)) {
        boot_result <- tryCatch({
            unique_countries <- unique(data_subset$ISO3)
            n_countries <- length(unique_countries)

            # Per-replicate error handling keeps successful bootstrap draws
            # even when some cluster resamples fail to fit.
            boot_values <- replicate(boot_nsim, {
                tryCatch({
                    sampled_countries <- sample(unique_countries, n_countries, replace = TRUE)
                    idx <- unlist(lapply(sampled_countries, function(iso) {
                        which(data_subset$ISO3 == iso)
                    }))

                    bdf <- data_subset[idx, ]
                    mb <- fit_ppml_glm(bdf, extra_pcs = extra_pcs)
                    coef(mb)[["Consumption"]]
                }, error = function(e) {
                    NA_real_
                })
            })

            valid_boot_values <- boot_values[is.finite(boot_values)]
            invalid_n <- length(boot_values) - length(valid_boot_values)

            if (invalid_n > 0) {
                warning(
                    "PPML cluster bootstrap dropped ", invalid_n,
                    " of ", length(boot_values), " resamples due to fit failures."
                )
            }

            if (length(valid_boot_values) >= 2) {
                list(
                    lower_ci = quantile(valid_boot_values, 0.025, na.rm = TRUE),
                    upper_ci = quantile(valid_boot_values, 0.975, na.rm = TRUE),
                    bootstrap_values = valid_boot_values
                )
            } else {
                intervals <- suppressMessages(confint.default(model))
                warning(
                    "PPML cluster bootstrap produced fewer than 2 valid resamples; ",
                    "falling back to Wald CI for interval estimation."
                )
                list(
                    lower_ci = intervals["Consumption", 1],
                    upper_ci = intervals["Consumption", 2],
                    bootstrap_values = if (length(valid_boot_values) == 1) valid_boot_values else gradient
                )
            }
        }, error = function(e) {
            warning("Cluster bootstrap failed for PPML Consumption term. Error: ", conditionMessage(e))
            intervals <- suppressMessages(confint.default(model))
            list(
                lower_ci = intervals["Consumption", 1],
                upper_ci = intervals["Consumption", 2],
                bootstrap_values = gradient
            )
        })
        lower_ci <- boot_result$lower_ci
        upper_ci <- boot_result$upper_ci
        bootstrap_values <- boot_result$bootstrap_values
    } else {
        intervals <- suppressMessages(confint.default(model))
        lower_ci <- intervals["Consumption", 1]
        upper_ci <- intervals["Consumption", 2]
        bootstrap_values <- gradient
    }
    
    calc_se <- if (length(bootstrap_values) > 1) sd(bootstrap_values, na.rm = TRUE) else summary(model)$coefficients["Consumption", 2]

    list(
        gradient = gradient,
        intercept = intercept,
        lower_ci = lower_ci,
        upper_ci = upper_ci,
        se = calc_se,
        bootstrap_values = bootstrap_values
    )
}

extract_glm_binomial_consumption_summary <- function(model, boot_nsim, data_subset = NULL, extra_pcs = FALSE) {
    gradient <- coef(model)[["Consumption"]]
    intercept <- coef(model)[["(Intercept)"]]

    # Fast/stable CI for pipeline runs
    ci <- suppressMessages(confint.default(model))
    lower_ci <- ci["Consumption", 1]
    upper_ci <- ci["Consumption", 2]

    if (boot_nsim > 0) {
        boot_values <- tryCatch(
            {
                n <- nrow(data_subset)
                replicate(boot_nsim, {
                    idx <- sample(n, n, replace = TRUE)
                    bdf <- data_subset[idx, ]
                    mb <- glm(
                        formula = get_fixed_effects_formula(extra_pcs = extra_pcs),
                        data = bdf,
                        family = binomial(link = "logit")
                    )
                    coef(mb)[["Consumption"]]
                })
            },
            error = function(e) gradient
        )
    } else {
        boot_values <- gradient
    }

    list(
        gradient = gradient,
        intercept = intercept,
        lower_ci = lower_ci,
        upper_ci = upper_ci,
        bootstrap_values = boot_values
    )
}

extract_lmer_consumption_summary <- function(model, boot_nsim, data_subset = NULL) {
    gradient <- summary(model)$coefficients["Consumption", 1]
    intercept <- summary(model)$coefficients["(Intercept)", 1]
    
    # NEW: Determine exactly how many fixed effects the healthy main model has
    expected_coef_count <- length(lme4::fixef(model))

    if (boot_nsim > 0 && !is.null(data_subset)) {
        boot_result <- tryCatch(
            {
                unique_countries <- unique(data_subset$ISO3)
                n_countries <- length(unique_countries)
                
                boot_values <- replicate(boot_nsim, {
                    
                    # NEW: Retry mechanism (prevents throwing away too many iterations)
                    max_attempts <- 10
                    attempt <- 1
                    valid_estimate <- NA
                    
                    while(attempt <= max_attempts && is.na(valid_estimate)) {
                        
                        # 1. Resample countries
                        sampled_countries <- sample(unique_countries, n_countries, replace = TRUE)
                        idx <- unlist(lapply(sampled_countries, function(iso) {
                            which(data_subset$ISO3 == iso)
                        }))
                        
                        bdf <- data_subset[idx, ]
                        
                        valid_estimate <- tryCatch({
                            # 2. Fit the mixed model
                            mb <- lme4::lmer(
                                formula = formula(model),
                                data = bdf,
                                weights = Weight,
                                control = lme4::lmerControl(calc.derivs = FALSE)
                            )
                            
                            # NEW: Rank deficiency check for LMER
                            # If the bootstrap model has fewer coefficients than the main model, a factor was dropped.
                            if (length(lme4::fixef(mb)) != expected_coef_count) {
                                stop("Rank deficient (dropped factors) in LMER fit.")
                            }
                            
                            # 3. Handle Singularities
                            if (lme4::isSingular(mb)) {
                                mb_fallback <- lm(
                                    formula = lme4::nobars(formula(model)), 
                                    data = bdf,
                                    weights = Weight
                                )
                                
                                # NEW: Rank deficiency check for LM fallback
                                # In base R 'lm', dropped factors result in NA coefficients.
                                if (any(is.na(coef(mb_fallback)))) {
                                    stop("Rank deficient (dropped factors) in LM fallback.")
                                }
                                
                                coef(mb_fallback)["Consumption"]
                                
                            } else {
                                lme4::fixef(mb)["Consumption"]
                            }
                            
                        }, error = function(e) {
                            # If an error or stop() occurs, return NA to trigger the next while-loop attempt
                            NA 
                        })
                        
                        attempt <- attempt + 1
                    }
                    
                    # If it failed 10 times in a row, it returns NA and moves to the next bootstrap iteration
                    if (is.na(valid_estimate)) {
                        warning("Bootstrap iteration completely failed after 10 attempts due to data structure/rank deficiency.")
                    }
                    
                    return(valid_estimate)
                })
                
                # Clean out any persistent NAs before calculating quantiles
                valid_boot_values <- boot_values[!is.na(boot_values)]

                list(
                    lower_ci = quantile(valid_boot_values, 0.025, na.rm = TRUE),
                    upper_ci = quantile(valid_boot_values, 0.975, na.rm = TRUE),
                    bootstrap_values = boot_values
                )
            },
            error = function(e) {
                warning("Cluster bootstrap failed entirely; using Wald CI fallback. Error: ", conditionMessage(e))
                intervals <- suppressMessages(confint(model, parm = "Consumption", method = "Wald"))
                list(
                    lower_ci = intervals[1, 1],
                    upper_ci = intervals[1, 2],
                    bootstrap_values = gradient
                )
            }
        )
        lower_ci <- boot_result$lower_ci
        upper_ci <- boot_result$upper_ci
        bootstrap_values <- boot_result$bootstrap_values
    } else {
        intervals <- suppressMessages(confint(model, parm = "Consumption", method = "Wald"))
        lower_ci <- intervals[1, 1]
        upper_ci <- intervals[1, 2]
        bootstrap_values <- gradient
    }

    # Calculate SE from bootstraps, fallback to Wald SE if bootstraps failed/0
    if (boot_nsim > 0 && !is.null(data_subset) && length(bootstrap_values) > 1) {
        calc_se <- sd(bootstrap_values, na.rm = TRUE)
    } else {
        # Fallback to analytic standard error
        calc_se <- summary(model)$coefficients["Consumption", 2] 
    }

    list(
        gradient = gradient,
        intercept = intercept,
        lower_ci = lower_ci,
        upper_ci = upper_ci,
        se = calc_se, # NEW: Standard error for MI
        bootstrap_values = bootstrap_values
    )
}

extract_glmer_binomial_consumption_summary <- function(model, boot_nsim) {
    gradient <- lme4::fixef(model)[["Consumption"]]
    intercept <- lme4::fixef(model)[["(Intercept)"]]

    if (boot_nsim > 0) {
        boot_result <- tryCatch(
            {
                bs <- bootMer(model, FUN = function(x) fixef(x)[["Consumption"]], nsim = boot_nsim)
                ci <- confint(bs)
                list(
                    lower_ci = ci["t*", "2.5 %"],
                    upper_ci = ci["t*", "97.5 %"],
                    bootstrap_values = as.vector(bs$t)
                )
            },
            error = function(e) {
                ci <- suppressMessages(confint(model, parm = "Consumption", method = "Wald"))
                list(
                    lower_ci = ci[1, 1],
                    upper_ci = ci[1, 2],
                    bootstrap_values = gradient
                )
            }
        )
    } else {
        ci <- suppressMessages(confint(model, parm = "Consumption", method = "Wald"))
        boot_result <- list(
            lower_ci = ci[1, 1],
            upper_ci = ci[1, 2],
            bootstrap_values = gradient
        )
    }

    list(
        gradient = gradient,
        intercept = intercept,
        lower_ci = boot_result$lower_ci,
        upper_ci = boot_result$upper_ci,
        bootstrap_values = boot_result$bootstrap_values
    )
}

extract_glmer_ppml_consumption_summary <- function(model, boot_nsim, data_subset = NULL, extra_pcs = FALSE) {
    # glmmTMB fixed effects for the conditional model are in a list accessed via $cond
    gradient <- glmmTMB::fixef(model)$cond["Consumption"]
    intercept <- glmmTMB::fixef(model)$cond["(Intercept)"]

    if (boot_nsim > 0 && !is.null(data_subset)) {
        fallback_fixed_ppml_bootstrap <- function() {
            # Last-resort path: use fixed-effects PPML bootstrap distribution
            # rather than collapsing to a single point estimate.
            fallback_model <- fit_ppml_glm(data_subset, extra_pcs = extra_pcs)
            fallback_stats <- extract_glm_ppml_consumption_summary(
                fallback_model,
                boot_nsim = boot_nsim,
                data_subset = data_subset,
                extra_pcs = extra_pcs
            )
            list(
                lower_ci = fallback_stats$lower_ci,
                upper_ci = fallback_stats$upper_ci,
                bootstrap_values = fallback_stats$bootstrap_values
            )
        }

        boot_result <- tryCatch({
            unique_countries <- unique(data_subset$ISO3)
            n_countries <- length(unique_countries)

            boot_values <- replicate(boot_nsim, {
                sampled_countries <- sample(unique_countries, n_countries, replace = TRUE)
                idx <- unlist(lapply(sampled_countries, function(iso) {
                    which(data_subset$ISO3 == iso)
                }))
                bdf <- data_subset[idx, ]

                # Try random-effects PPML first; if it is non-converged or fails,
                # fall back to fixed-effects PPML for that resample.
                estimate <- tryCatch({
                    mb <- glmmTMB::glmmTMB(
                        formula = formula(model),
                        data = bdf,
                        family = poisson(link = "log"),
                        weights = Weight
                    )

                    if (isTRUE(mb$sdr$pdHess) &&
                        "Consumption" %in% names(glmmTMB::fixef(mb)$cond) &&
                        is.finite(glmmTMB::fixef(mb)$cond[["Consumption"]])) {
                        glmmTMB::fixef(mb)$cond[["Consumption"]]
                    } else {
                        mb_fallback <- fit_ppml_glm(bdf, extra_pcs = extra_pcs)
                        coef(mb_fallback)[["Consumption"]]
                    }
                }, error = function(e) {
                    NA_real_
                })

                if (!is.finite(estimate)) {
                    estimate <- tryCatch({
                        mb_fallback <- fit_ppml_glm(bdf, extra_pcs = extra_pcs)
                        coef(mb_fallback)[["Consumption"]]
                    }, error = function(e) {
                        NA_real_
                    })
                }

                estimate
            })

            valid_boot_values <- boot_values[is.finite(boot_values)]
            invalid_n <- length(boot_values) - length(valid_boot_values)

            if (invalid_n > 0) {
                message(
                    "PPML random-effects bootstrap dropped ", invalid_n,
                    " of ", length(boot_values), " resamples due to fit failures."
                )
            }

            if (length(valid_boot_values) >= 2) {
                list(
                    lower_ci = quantile(valid_boot_values, 0.025, na.rm = TRUE),
                    upper_ci = quantile(valid_boot_values, 0.975, na.rm = TRUE),
                    bootstrap_values = valid_boot_values
                )
            } else {
                warning(
                    "PPML random-effects bootstrap produced fewer than 2 valid resamples; ",
                    "attempting fixed-effects PPML bootstrap fallback."
                )
                fallback_fixed_ppml_bootstrap()
            }
        }, error = function(e) {
            warning(
                "PPML random-effects bootstrap failed entirely; attempting fixed-effects PPML bootstrap fallback. Error: ",
                conditionMessage(e)
            )
            fallback_fixed_ppml_bootstrap()
        })
        lower_ci <- boot_result$lower_ci
        upper_ci <- boot_result$upper_ci
        bootstrap_values <- boot_result$bootstrap_values
    } else {
        intervals <- suppressMessages(confint(model))
        # Dynamically find the correct row name
        rn <- intersect(c("Consumption", "cond.Consumption"), rownames(intervals))[1]
        lower_ci <- intervals[rn, 1]
        upper_ci <- intervals[rn, 2]
        bootstrap_values <- gradient
    }

    calc_se <- if (length(bootstrap_values) > 1) sd(bootstrap_values, na.rm = TRUE) else summary(model)$coefficients$cond["Consumption", 2]

    list(gradient = gradient, intercept = intercept, lower_ci = lower_ci, upper_ci = upper_ci, se = calc_se, bootstrap_values = bootstrap_values)
}

build_bootstrap_df <- function(label_name, label_value, bootstrap_values) {
    outdf <- data.frame(Gradient = bootstrap_values)
    outdf[[label_name]] <- label_value
    outdf[, c(label_name, "Gradient"), drop = FALSE]
}

fit_combined_pathogen_drug_lm <- function(data_, output_tag = "lagged", runtime_options = get_runtime_options(), output_prefix = "database", model_family = "gaussian", extra_pcs = FALSE) {
    gradients <- c()
    conf_intervals <- c()
    pathogens <- c()
    abs <- c()
    bootstraps <- data.frame()
    r_squareds <- c()
    variation_explained_list <- list()
    vif_list <- list()
    white_tests <- c()
    RMSEs <- c()
    AICs <- c()
    BICs <- c()

    antibiotics_to_fit <- sort(unique(data_$Antibiotic))
    pathogens_to_fit <- sort(unique(data_$Pathogen))

    if (isTRUE(runtime_options$smoke_mode)) {
        antibiotics_to_fit <- head(antibiotics_to_fit, runtime_options$smoke_max_classes)
        pathogens_to_fit <- head(pathogens_to_fit, runtime_options$smoke_max_pathogens)
    }

    for (antibiotic in antibiotics_to_fit) {
        for (pathogen in pathogens_to_fit) {
            data_subset <- data_[data_$Pathogen == pathogen & data_$Antibiotic == antibiotic, ]
            if (extra_pcs) {
                vars_needed <- c("Resistance", "Consumption", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "GDP", "Year", "Weight")
            } else {
                vars_needed <- c("Resistance", "Consumption", "PC1", "PC2", "PC3", "GDP", "Year", "Weight")
            }
            if (model_family == "binomial") {
                vars_needed <- c(vars_needed, "resistant_count")
            }
            data_subset <- data_subset[complete.cases(data_subset[, vars_needed]), ]
            data_subset <- data_subset[!is.infinite(data_subset$Consumption), ]
            data_subset <- data_subset[!is.infinite(data_subset$Resistance), ]
            if (nrow(data_subset) <= 1) {
                next
            }

            if (model_family == "gaussian") {
                model <- fit_weighted_lm(data_subset, extra_pcs = extra_pcs)
            } else if (model_family == "binomial") {
                model <- fit_binomial_glm(data_subset, extra_pcs = extra_pcs)
            } else if (model_family == "ppml") {
                model <- fit_ppml_glm(data_subset, extra_pcs = extra_pcs)
            } else {
                stop("Unsupported model family: ", model_family)
            }

            # calculate White test for heteroskedasticity
            white_test <- white(model)$p.value
            RMSE <- sqrt(mean(model$residuals^2))
            AIC <- AIC(model)
            BIC <- BIC(model)

            # # save a plot of residuals
            # plot_path <- paste0("residuals_", antibiotic, "_", pathogen, "_year_colors_nonppml.png")
            # png(plot_path)
            # # color by year
            # plot(model$residuals, main = paste("Residuals for", antibiotic, "and", pathogen), ylab = "Residuals", col = rainbow(length(unique(data_subset$Year)))[as.numeric(as.factor(data_subset$Year))])
            # abline(h = 0, col = "red")
            # dev.off()

            # #another version colored by country
            # plot_path_country <- paste0("residuals_", antibiotic, "_", pathogen, "_country_colors.png")
            # png(plot_path_country)
            # plot(model$residuals, main = paste("Residuals for", antibiotic, "and", pathogen), ylab = "Residuals", col = rainbow(length(unique(data_subset$ISO3)))[as.numeric(as.factor(data_subset$ISO3))])
            # abline(h = 0, col = "red")
            # dev.off()

            # # plot cooks distance
            # plot_path_cooks <- paste0("cooks_distance_", antibiotic, "_", pathogen, "_year_colors.png")
            # png(plot_path_cooks)
            # plot(cooks.distance(model), main = paste("Cook's Distance for", antibiotic, "and", pathogen), ylab = "Cook's Distance", col = rainbow(length(unique(data_subset$Year)))[as.numeric(as.factor(data_subset$Year))])
            # abline(h = 1, col = "red")
            # abline(h = 4/(nrow(data_subset)-length(coef(model))-2), col = "red")
            # dev.off()

            # Compute R-squared and per-variable variance explained
            if (extra_pcs) {
                ve_vars <- c("Consumption", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "GDP", "Year")
            } else {
                ve_vars <- c("Consumption", "PC1", "PC2", "PC3", "GDP", "Year")
            }
            ve <- setNames(rep(NA_real_, length(ve_vars)), ve_vars)
            vifs <- setNames(rep(NA_real_, length(ve_vars)), ve_vars)
            
            tryCatch({
                # 1. Check if the model dropped any variables (NA coefficients)
                na_coefs <- names(coef(model))[is.na(coef(model))]
                
                if (length(na_coefs) > 0) {
                    # Map the dropped coefficients back to the original variable names
                    # (e.g., "Year2018" matches back to "Year")
                    aliased_vars <- ve_vars[sapply(ve_vars, function(v) any(grepl(v, na_coefs)))]
                    
                    if (length(aliased_vars) > 0) {
                        # Create a safe list of variables, removing the aliased ones
                        safe_vars <- setdiff(ve_vars, aliased_vars)
                        
                        # Refit a temporary model just for VIF calculation
                        if (model_family == "gaussian") {
                            safe_formula <- as.formula(paste("Resistance ~", paste(safe_vars, collapse = " + ")))
                            safe_model <- lm(safe_formula, data = data_subset, weights = data_subset$Weight)
                        } else if (model_family == "binomial") {
                            safe_formula <- as.formula(paste("cbind(resistant_count, Weight - resistant_count) ~", 
                                                             paste(safe_vars, collapse = " + ")))
                            safe_model <- glm(safe_formula, data = data_subset, family = binomial(link = "logit"))
                        } else if (model_family == "ppml") {
                            safe_formula <- as.formula(paste("Resistance ~", paste(safe_vars, collapse = " + ")))
                            safe_model <- glm(safe_formula, data = data_subset, weights = data_subset$Weight, family = quasipoisson(link = "log"))
                        }
                        calc_vifs <- car::vif(safe_model)
                    } else {
                        calc_vifs <- car::vif(model)
                    }
                } else {
                    # No aliased coefficients, calculate normally
                    calc_vifs <- car::vif(model)
                }
                
                # Safely map calculated VIFs into our array
                matching_vars <- intersect(names(calc_vifs), names(vifs))
                vifs[matching_vars] <- calc_vifs[matching_vars]
                
            }, error = function(e) {
                # If VIF still fails, print a warning to the console but KEEP RUNNING the loop.
                # The VIFs for this specific model will just remain as NAs.
                message("Notice: VIF calculation skipped for ", antibiotic, " / ", pathogen, " - ", e$message)
            })

            if (model_family == "gaussian") {
                r_sq <- summary(model)$r.squared
            } else if (model_family %in% c("binomial", "ppml")) {
                r_sq <- 1 - (model$deviance / model$null.deviance)
            }

            # Estimate predictor contribution from Type II ANOVA statistics.
            tryCatch({
                anova_tbl <- tryCatch(
                    car::Anova(model, type = 2),
                    error = function(e) car::Anova(model, type = 2, test.statistic = "Wald")
                )
                anova_df <- as.data.frame(anova_tbl)
                stat_cols <- c("Sum Sq", "Chisq", "LR Chisq", "Wald", "F value", "F")
                stat_col <- stat_cols[stat_cols %in% colnames(anova_df)][1]

                if (!is.na(stat_col)) {
                    term_stats <- anova_df[[stat_col]]
                    names(term_stats) <- rownames(anova_df)
                    term_stats <- term_stats[names(term_stats) %in% ve_vars]

                    stat_total <- sum(term_stats, na.rm = TRUE)
                    if (is.finite(stat_total) && stat_total > 0) {
                        ve[names(term_stats)] <- term_stats / stat_total
                    }
                }
            }, error = function(e) {
                message("Notice: ANOVA contribution calculation skipped for ", antibiotic, " / ", pathogen, " - ", e$message)
            })

            if (is.na(coef(model)["Consumption"])) {
                next
            }

            if (model_family == "gaussian") {
                summary_stats <- extract_lm_consumption_summary(
                    model,
                    runtime_options$boot_nsim,
                    data_subset,
                    extra_pcs = extra_pcs
                )
            } else if (model_family == "binomial") {
                summary_stats <- extract_glm_binomial_consumption_summary(
                    model,
                    runtime_options$boot_nsim,
                    data_subset,
                    extra_pcs = extra_pcs
                )
            } else if (model_family == "ppml") {
                summary_stats <- extract_glm_ppml_consumption_summary(model, runtime_options$boot_nsim, data_subset, extra_pcs = extra_pcs)
            }

            gradients <- c(gradients, summary_stats$gradient)
            conf_intervals <- c(conf_intervals, summary_stats$lower_ci,
                                summary_stats$upper_ci)
            pathogens <- c(pathogens, pathogen)
            abs <- c(abs, antibiotic)
            r_squareds <- c(r_squareds, r_sq)
            variation_explained_list[[length(variation_explained_list) + 1]] <- ve
            vif_list[[length(vif_list) + 1]] <- vifs
            white_tests <- c(white_tests, white_test)
            RMSEs <- c(RMSEs, RMSE)
            AICs <- c(AICs, AIC)
            BICs <- c(BICs, BIC)

            outdf <- data.frame(
                Pathogen = pathogen,
                Antibiotic = antibiotic,
                Gradient = summary_stats$bootstrap_values
            )
            bootstraps <- rbind(bootstraps, outdf)
        }
    }

    if (length(gradients) == 0) {
        warning("No valid combined pathogen-drug LM models were fit for output tag: ", output_tag)
        return(invisible(NULL))
    }

    conf_intervals <- matrix(conf_intervals, nrow = length(gradients), ncol = 2, byrow = TRUE)
    variation_explained_df <- do.call(rbind, lapply(variation_explained_list, function(x) as.data.frame(t(x))))
    colnames(variation_explained_df) <- paste0("Variation_Explained.", colnames(variation_explained_df))
    vif_df <- do.call(rbind, lapply(vif_list, function(x) as.data.frame(t(x))))
    colnames(vif_df) <- paste0("VIF.", colnames(vif_df))
    results <- data.frame(
        Antibiotic = abs,
        Pathogen = pathogens,
        Response = gradients,
        Lower_CI = conf_intervals[, 1],
        Upper_CI = conf_intervals[, 2],
        R_squared = r_squareds,
        White_test_p_value = white_tests,
        RMSE = RMSEs,
        AIC = AICs,
        BIC = BICs
    )
    results <- cbind(results, variation_explained_df, vif_df)

    gradient_prefix <- if (identical(output_prefix, "Nagorsen")) {
        "Nagorsen_gradients_pathogen_ATC3_PCA_canonical"
    } else {
        "database_gradients_pathogen_ATC3_PCA_canonical_weighted"
    }
    bootstrap_prefix <- if (identical(output_prefix, "Nagorsen")) {
        "Nagorsen_gradients_bootstraps_pathogen_ATC3_PCA_canonical"
    } else {
        "database_gradients_bootstraps_pathogen_ATC3_PCA_canonical_weighted"
    }

    write.csv(results, build_output_path(gradient_prefix, output_tag), row.names = FALSE)
    write.csv(bootstraps, build_output_path(bootstrap_prefix, output_tag), row.names = FALSE)
}

fit_combined_pathogen_drug_mi_lm <- function(mi_datasets, output_tag = "lagged", runtime_options = get_runtime_options(), output_prefix = "database", model_family = "gaussian", extra_pcs = FALSE) {
    
    gradients <- c()
    conf_intervals <- c()
    pathogens <- c()
    abs <- c()
    
    # Use the first dataset to define the combinations to fit
    template_data <- mi_datasets[[1]]
    antibiotics_to_fit <- sort(unique(template_data$Antibiotic))
    pathogens_to_fit <- sort(unique(template_data$Pathogen))
    
    for (antibiotic in antibiotics_to_fit) {
        for (pathogen in pathogens_to_fit) {
            
            mi_estimates <- numeric(length(mi_datasets))
            mi_ses <- numeric(length(mi_datasets))
            
            for (i in seq_along(mi_datasets)) {
                data_subset <- mi_datasets[[i]][mi_datasets[[i]]$Pathogen == pathogen & mi_datasets[[i]]$Antibiotic == antibiotic, ]
                
                if (nrow(data_subset) <= 1) {
                    mi_estimates[i] <- NA
                    mi_ses[i] <- NA
                    next
                }
                
                # Fit the model on this specific imputation
                model <- fit_weighted_lm(data_subset, extra_pcs = extra_pcs)
                
                # Extract gradient and bootstrap standard error
                summary_stats <- extract_lm_consumption_summary(
                    model,
                    runtime_options$boot_nsim, # This is now hardcoded to 100
                    data_subset,
                    extra_pcs = extra_pcs
                )
                
                mi_estimates[i] <- summary_stats$gradient
                mi_ses[i] <- summary_stats$se
            }
            
            # Pool the 20 results using Rubin's Rules
            pooled_results <- pool_rubins_rules(mi_estimates, mi_ses)
            
            # Skip if pooling failed (e.g., all models crashed)
            if (is.na(pooled_results$gradient)) next
            
            gradients <- c(gradients, pooled_results$gradient)
            conf_intervals <- c(conf_intervals, pooled_results$lower_ci, pooled_results$upper_ci)
            pathogens <- c(pathogens, pathogen)
            abs <- c(abs, antibiotic)
        }
    }
    
    conf_intervals <- matrix(conf_intervals, nrow = length(gradients), ncol = 2, byrow = TRUE)
    results <- data.frame(
        Antibiotic = abs,
        Pathogen = pathogens,
        Response = gradients,
        Lower_CI = conf_intervals[, 1],
        Upper_CI = conf_intervals[, 2]
    )
    
    gradient_prefix <- "database_MI_pooled_gradients"
    write.csv(results, build_output_path(gradient_prefix, output_tag), row.names = FALSE)
}

fit_combined_pathogen_drug_glmnet <- function(data_, output_tag = "lagged", runtime_options = get_runtime_options(), output_prefix = "database") {
    gradients <- c()
    conf_intervals <- c()
    pathogens <- c()
    abs <- c()
    lambdas <- c()
    bootstraps <- data.frame()
    cross_class_effects <- data.frame()

    antibiotics_to_fit <- sort(unique(data_$Antibiotic))
    pathogens_to_fit <- sort(unique(data_$Pathogen))

    if (isTRUE(runtime_options$smoke_mode)) {
        antibiotics_to_fit <- head(antibiotics_to_fit, runtime_options$smoke_max_classes)
        pathogens_to_fit <- head(pathogens_to_fit, runtime_options$smoke_max_pathogens)
    }

    for (antibiotic in antibiotics_to_fit) {
        for (pathogen in pathogens_to_fit) {
            # 1. Subset and clean data
            data_subset <- data_[data_$Pathogen == pathogen & data_$Antibiotic == antibiotic, ]
            data_subset <- data_subset %>%
                select(where(~ !all(is.na(.)))) %>%
                na.omit()
            
            if (nrow(data_subset) < 10) {
                next
            }

            # 2. Build Primary Matrix
            y_vector <- data_subset$Resistance
            weights_vector <- data_subset$Weight
            
            # Drop character metadata completely so model.matrix doesn't see them
            x_data <- data_subset %>%
                select(-Pathogen, -Antibiotic, -ISO3, -Weight, -Resistance)

            # Ensure Year is strictly numeric to act as a continuous secular trend
            if ("Year" %in% colnames(x_data)) {
                x_data$Year <- as.numeric(as.character(x_data$Year))
            }

            # Build matrix purely on the numeric predictor columns
            x_matrix <- model.matrix(~ . - 1, data = x_data)

            # --- NEW: Create Selective Penalty Factor ---
            # Default all columns to a penalty of 1 (fully penalized)
            p_fac <- rep(1, ncol(x_matrix))
            
            # # Find the index of our target "Consumption" column and set its penalty to 0
            # target_idx <- which(colnames(x_matrix) == "Consumption")
            # if (length(target_idx) > 0) {
            #     p_fac[target_idx] <- 0 
            # } else {
            #     warning("Consumption column not found in x_matrix for ", pathogen, " x ", antibiotic, "; proceeding without unpenalized term.")
            # }

            # 3. Fit Primary CV Model to find optimal Lambda (alpha = 0 for Ridge, alpha = 1 for Lasso)
            cv_model <- tryCatch({
                cv.glmnet(x = x_matrix, y = y_vector, weights = weights_vector,
                alpha = 1, penalty.factor = p_fac)
            }, error = function(e) NULL)

            if (is.null(cv_model)) {
                next
            }

            # Extract point estimate at lambda.min
            best_lambda <- cv_model$lambda.min
            model_coefs <- as.matrix(coef(cv_model, s = best_lambda))
            
            # Extract focal gradient
            gradient <- if ("Consumption" %in% rownames(model_coefs)) model_coefs["Consumption", 1] else 0 

            # Extract cross-class coefficients that survived the penalty
            cross_vars <- grep("\\.Consumption$", rownames(model_coefs), value = TRUE)
            
            # --- CRITICAL FIX: Force R to keep names even if length == 1 ---
            cross_coefs <- setNames(as.numeric(model_coefs[cross_vars, 1]), cross_vars)
            # ---------------------------------------------------------------
            
            active_cross_coefs <- cross_coefs[cross_coefs != 0]

            # Define all variables we want the bootstrap to track
            vars_to_track <- c("Consumption", names(active_cross_coefs))
            
            # 5. Bootstrap Loop
            if (runtime_options$boot_nsim > 0) {
                n <- nrow(data_subset)
                
                # Run resamples using lapply to safely build a matrix
                boot_list <- lapply(1:runtime_options$boot_nsim, function(i) {
                    idx <- sample(n, n, replace = TRUE)
                    
                    x_boot_data <- data_subset[idx, ] %>%
                        select(-Pathogen, -Antibiotic, -ISO3, -Weight, -Resistance)
                    if ("Year" %in% colnames(x_boot_data)) {
                        x_boot_data$Year <- as.numeric(as.character(x_boot_data$Year))
                    }
                    
                    x_boot <- model.matrix(~ . - 1, data = x_boot_data)
                    y_boot <- data_subset$Resistance[idx]
                    w_boot <- data_subset$Weight[idx]
                    
                    p_fac_boot <- rep(1, ncol(x_boot))
                    target_idx_boot <- which(colnames(x_boot) == "Consumption")
                    if (length(target_idx_boot) > 0) p_fac_boot[target_idx_boot] <- 0 
                    
                    b_model <- tryCatch({
                        glmnet(x = x_boot, y = y_boot, weights = w_boot, 
                               alpha = 1, lambda = best_lambda, penalty.factor = p_fac_boot)
                    }, error = function(e) NULL)
                    
                    if (is.null(b_model)) return(NA)
                    
                    b_coefs <- as.matrix(coef(b_model))
                    
                    # Extract the exact tracked variables (returning 0 if Lasso dropped them)
                    sapply(vars_to_track, function(v) {
                        if (v %in% rownames(b_coefs)) b_coefs[v, 1] else 0
                    })
                })
                
                # Remove failed iterations and combine into a matrix
                valid_boots <- boot_list[!is.na(boot_list)]
                if (length(valid_boots) > 0) {
                    boot_matrix <- do.call(cbind, valid_boots)
                    
                    cis <- apply(boot_matrix, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
                    
                    focal_lower <- cis[1, "Consumption"]
                    focal_upper <- cis[2, "Consumption"]
                    cross_lowers <- cis[1, names(active_cross_coefs)]
                    cross_uppers <- cis[2, names(active_cross_coefs)]
                    
                    outdf <- data.frame(Pathogen = pathogen, Antibiotic = antibiotic, Gradient = boot_matrix["Consumption", ])
                } else {
                    focal_lower <- NA; focal_upper <- NA
                    cross_lowers <- rep(NA_real_, length(active_cross_coefs))
                    cross_uppers <- rep(NA_real_, length(active_cross_coefs))
                    outdf <- data.frame(Pathogen = pathogen, Antibiotic = antibiotic, Gradient = gradient)
                }
            } else {
                focal_lower <- NA; focal_upper <- NA
                # Use NA_real_ to strictly ensure numeric columns
                cross_lowers <- rep(NA_real_, length(active_cross_coefs)) 
                cross_uppers <- rep(NA_real_, length(active_cross_coefs))
                outdf <- data.frame(Pathogen = pathogen, Antibiotic = antibiotic, Gradient = gradient)
            }

            # --- Store Cross-Class CIs alongside the effect ---
            if (length(active_cross_coefs) > 0) {
                cross_df <- data.frame(
                    Pathogen = pathogen,
                    Target_Antibiotic = antibiotic,
                    Cross_Class_Antibiotic = names(active_cross_coefs),
                    Coefficient = as.numeric(active_cross_coefs), # Force numeric
                    Lower_CI = as.numeric(cross_lowers),
                    Upper_CI = as.numeric(cross_uppers),
                    row.names = NULL
                )
                cross_class_effects <- rbind(cross_class_effects, cross_df)
            }

            # 5. Store Results
            gradients <- c(gradients, gradient)
            conf_intervals <- c(conf_intervals, focal_lower, focal_upper)
            pathogens <- c(pathogens, pathogen)
            abs <- c(abs, antibiotic)
            lambdas <- c(lambdas, best_lambda)
            bootstraps <- rbind(bootstraps, outdf)
        }
    }

    if (length(gradients) == 0) {
        warning("No valid glmnet models were fit for output tag: ", output_tag)
        return(invisible(NULL))
    }

    # Format the confidence intervals into a matrix
    conf_intervals <- matrix(conf_intervals, nrow = length(gradients), ncol = 2, byrow = TRUE)

    results <- data.frame(
        Antibiotic = abs,
        Pathogen = pathogens,
        Response = gradients,
        Lower_CI = conf_intervals[, 1],
        Upper_CI = conf_intervals[, 2],
        Optimal_Lambda = lambdas
    )

    # Output Paths
    gradient_prefix <- if (identical(output_prefix, "Nagorsen")) {
        "Nagorsen_gradients_pathogen_ATC3_glmnet"
    } else {
        "database_gradients_pathogen_ATC3_glmnet_weighted"
    }
    cross_prefix <- if (identical(output_prefix, "Nagorsen")) {
        "Nagorsen_cross_class_effects_glmnet"
    } else {
        "database_cross_class_effects_glmnet_weighted"
    }
    bootstrap_prefix <- if (identical(output_prefix, "Nagorsen")) {
        "Nagorsen_gradients_bootstraps_pathogen_ATC3_glmnet"
    } else {
        "database_gradients_bootstraps_pathogen_ATC3_glmnet_weighted"
    }

    # Write files
    write.csv(results, build_output_path(gradient_prefix, output_tag), row.names = FALSE)
    write.csv(bootstraps, build_output_path(bootstrap_prefix, output_tag), row.names = FALSE)
    # Write cross-class effects if any were found
    if (nrow(cross_class_effects) > 0) {
        cross_class_effects$Cross_Class_Antibiotic <- gsub("\\.Consumption$", "", cross_class_effects$Cross_Class_Antibiotic)
        write.csv(cross_class_effects, build_output_path(cross_prefix, output_tag), row.names = FALSE)
    }
}

fit_class_random_effects_models <- function(data_, output_tag = "all_lagged", runtime_options = get_runtime_options(), allow_fallback = FALSE) {
    fit_random_effects_models(
        data_ = data_,
        output_tag = output_tag,
        runtime_options = runtime_options,
        mode = "class",
        allow_fallback = allow_fallback
    )
}

fit_pathogen_random_effects_models <- function(data_, output_tag = "all_lagged", runtime_options = get_runtime_options(), allow_fallback = FALSE) {
    fit_random_effects_models(
        data_ = data_,
        output_tag = output_tag,
        runtime_options = runtime_options,
        mode = "pathogen",
        allow_fallback = allow_fallback
    )
}

fit_country_random_effects_models <- function(data_, output_tag = "all_lagged", runtime_options = get_runtime_options(), allow_fallback = FALSE) {
    fit_random_effects_models(
        data_ = data_,
        output_tag = output_tag,
        runtime_options = runtime_options,
        mode = "country",
        allow_fallback = allow_fallback
    )
}

fit_random_effects_models <- function(
    data_,
    output_tag,
    runtime_options = get_runtime_options(),
    mode = c("class", "pathogen", "country"),
    output_prefix = "database",
    allow_fallback = FALSE,
    model_family = "gaussian",
    extra_pcs = FALSE
) {
    mode <- match.arg(mode)
    accumulator <- initialize_random_effects_accumulator()

    if (mode == "class") {
        label_var <- "Antibiotic"
        random_effect_var <- "Pathogen"
        singular_msg <- "class"
        smoke_max <- runtime_options$smoke_max_classes
        if (identical(output_prefix, "Nagorsen")) {
            gradient_prefix <- "Nagorsen_gradients_ATC3_PCA_canonical"
            lower_prefix <- "Nagorsen_lowerCI_ATC3_PCA_canonical"
            upper_prefix <- "Nagorsen_upperCI_ATC3_PCA_canonical"
            bootstrap_prefix <- "Nagorsen_gradients_bootstraps_ATC3_PCA_canonical"
        } else {
            gradient_prefix <- "database_gradients_ATC3_PCA_canonical_weighted"
            lower_prefix <- "database_lowerCI_ATC3_PCA_canonical_weighted"
            upper_prefix <- "database_upperCI_ATC3_PCA_canonical_weighted"
            bootstrap_prefix <- "database_gradients_bootstraps_ATC3_PCA_canonical_weighted"
        }
    } else if (mode == "pathogen") {
        label_var <- "Pathogen"
        random_effect_var <- "Antibiotic"
        singular_msg <- "pathogen"
        smoke_max <- runtime_options$smoke_max_pathogens
        if (identical(output_prefix, "Nagorsen")) {
            gradient_prefix <- "Nagorsen_gradients_pathogen_PCA_canonical"
            lower_prefix <- "Nagorsen_lowerCI_pathogen_PCA_canonical"
            upper_prefix <- "Nagorsen_upperCI_pathogen_PCA_canonical"
            bootstrap_prefix <- "Nagorsen_gradients_bootstraps_pathogen_PCA_canonical"
        } else {
            gradient_prefix <- "database_gradients_pathogen_PCA_canonical_weighted"
            lower_prefix <- "database_lowerCI_pathogen_PCA_canonical_weighted"
            upper_prefix <- "database_upperCI_pathogen_PCA_canonical_weighted"
            bootstrap_prefix <- "database_gradients_bootstraps_pathogen_PCA_canonical_weighted"
        }
    } else if (mode == "country") {# Create a composite label directly in the dataset
        data_$Combo_Label <- paste(data_$Pathogen, data_$Antibiotic, sep = "___")
        label_var <- "Combo_Label"
        random_effect_var <- "ISO3"
        singular_msg <- "country"
        smoke_max <- Inf
        
        gradient_prefix <- if (identical(output_prefix, "Nagorsen")) "Nagorsen_gradients_country" else "database_gradients_country"
        lower_prefix    <- if (identical(output_prefix, "Nagorsen")) "Nagorsen_lowerCI_country" else "database_lowerCI_country"
        upper_prefix    <- if (identical(output_prefix, "Nagorsen")) "Nagorsen_upperCI_country" else "database_upperCI_country"
        bootstrap_prefix<- if (identical(output_prefix, "Nagorsen")) "Nagorsen_bootstraps_country" else "database_bootstraps_country"
    }

    labels <- sort(unique(data_[[label_var]]))

    if (isTRUE(runtime_options$smoke_mode)) {
        labels <- head(labels, smoke_max)
    }

    for (label in labels) {
        log_info("[ddd-linear-model] ", tools::toTitleCase(mode), " model: ", label, verbose = runtime_options$verbose)
        subset_data <- data_[data_[[label_var]] == label, ]
        subset_data <- subset_data[!is.infinite(subset_data$Consumption), ]
        subset_data <- subset_data[!is.infinite(subset_data$Resistance), ]

        if (nrow(subset_data) <= 1) {
            next
        }

        n_levels <- length(unique(subset_data[[random_effect_var]]))

        # Fit lmer if there are >2 levels. If singular/too few levels, check fallback.
        if (n_levels > 2 && nrow(subset_data) > n_levels) {
            if (model_family == "gaussian") {
                model <- fit_random_lmer(subset_data, random_effect_var = random_effect_var, extra_pcs = extra_pcs)
            } else if (model_family == "binomial") {
                model <- fit_binomial_glmer(subset_data, random_effect_var = random_effect_var, extra_pcs = extra_pcs)
            } else if (model_family == "ppml") {
                model <- fit_random_ppml_glmer(subset_data, random_effect_var = random_effect_var, extra_pcs = extra_pcs)
            } else {
                stop("Unsupported model family: ", model_family)
            }
            log_info(capture.output(print(model)), verbose = runtime_options$verbose)

            # Safely detect singularity/convergence issues for both package types
            is_singular_flag <- FALSE
            if (inherits(model, "merMod")) {
                is_singular_flag <- lme4::isSingular(model)
            } else if (inherits(model, "glmmTMB")) {
                # glmmTMB indicates singular fits with a non-positive-definite Hessian
                is_singular_flag <- !isTRUE(model$sdr$pdHess)
            }

            if (is_singular_flag) {
                if (allow_fallback) {
                    log_info(paste("Model is singular for", singular_msg, label, "- falling back to weighted lm"), verbose = runtime_options$verbose)
                    if (model_family == "gaussian") {
                        model <- fit_weighted_lm(subset_data, extra_pcs = extra_pcs)
                        summary_stats <- extract_lm_consumption_summary(model, runtime_options$boot_nsim, subset_data)
                    } else if (model_family == "binomial") {
                        model <- fit_binomial_glm(subset_data)
                        summary_stats <- extract_glm_binomial_consumption_summary(model, runtime_options$boot_nsim, subset_data)
                    } else if (model_family == "ppml") {
                        model <- fit_ppml_glm(subset_data, extra_pcs = extra_pcs)
                        summary_stats <- extract_glm_ppml_consumption_summary(model, runtime_options$boot_nsim, subset_data, extra_pcs = extra_pcs)
                    }
                } else {
                    log_info(paste("Model is singular for", singular_msg, label, "- excluding"), verbose = runtime_options$verbose)
                    next
                }
            } else {
                if (model_family == "ppml") {
                    summary_stats <- extract_glmer_ppml_consumption_summary(model, runtime_options$boot_nsim, subset_data, extra_pcs = extra_pcs)
                } else {
                    summary_stats <- extract_lmer_consumption_summary(model, runtime_options$boot_nsim, subset_data)
                }
            }
        } else {
            if (allow_fallback) {
                log_info(paste("Only", n_levels, "levels for", singular_msg, label, "- falling back to weighted lm"), verbose = runtime_options$verbose)
                if (model_family == "gaussian") {
                    model <- fit_weighted_lm(subset_data, extra_pcs = extra_pcs)
                    summary_stats <- extract_lm_consumption_summary(model, runtime_options$boot_nsim, subset_data)
                } else if (model_family == "binomial") {
                    model <- fit_binomial_glm(subset_data)
                    summary_stats <- extract_glm_binomial_consumption_summary(model, runtime_options$boot_nsim, subset_data)
                } else if (model_family == "ppml") {
                    model <- fit_ppml_glm(subset_data, extra_pcs = extra_pcs)
                    summary_stats <- extract_glm_ppml_consumption_summary(model, runtime_options$boot_nsim, subset_data, extra_pcs = extra_pcs)
                }
            } else {
                log_info(paste("Only", n_levels, "levels for", singular_msg, label, "- excluding"), verbose = runtime_options$verbose)
                next
            }
        }
        RMSE_val <- sqrt(mean(residuals(model)^2))
        AIC_val <- AIC(model)
        BIC_val <- BIC(model)
        # R-Squared
        r_sq <- tryCatch({
            if (inherits(model, "merMod") || inherits(model, "glmmTMB")) {
                suppressWarnings(MuMIn::r.squaredGLMM(model)[1, "R2c"]) 
            } else if (model_family == "gaussian") {
                summary(model)$r.squared
            } else {
                1 - (model$deviance / model$null.deviance)
            }
        }, error = function(e) NA_real_)
        # White Test
        white_p <- tryCatch({
            if (inherits(model, "merMod")) {
                # Proxy White test for lmer: regress squared residuals on fitted values
                aux_mod <- lm(residuals(model)^2 ~ poly(fitted(model), 2))
                f_stat <- summary(aux_mod)$fstatistic
                pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
            } else {
                skedastic::white(model)$p.value
            }
        }, error = function(e) NA_real_)
        # VIF
        vifs <- tryCatch({
            car::vif(model)
        }, error = function(e) {
            NA # Return NA if aliasing/rank deficiency breaks VIF
        })
        # Centralized appendage handles both lm and lmer uniformly
        outdf <- build_bootstrap_df(label_var, label, summary_stats$bootstrap_values)
        accumulator <- append_random_effects_result(
            accumulator = accumulator,
            label = label,
            gradient = summary_stats$gradient,
            intercept = summary_stats$intercept,
            lower_ci = summary_stats$lower_ci,
            upper_ci = summary_stats$upper_ci,
            bootstrap_df = outdf,
            r_sq = r_sq,
            white_p = white_p,
            rmse = RMSE_val,
            aic = AIC_val,
            bic = BIC_val,
            vifs = vifs
        )
        # # plot and save model residuals
        # png_filename <- paste0("residuals_", mode, "_", label, ".png")
        # png(png_filename, width = 800, height = 600)
        # plot(residuals(model), main = paste("Residuals for", mode, label))
        # abline(h = 0, col = "red")
        # dev.off()

        print(summary(model))
    }

    model_gradients <- setNames(accumulator$gradients, accumulator$labels)
    model_intercepts <- setNames(accumulator$intercepts, accumulator$labels)
    model_lower_ci <- setNames(accumulator$lower_ci, accumulator$labels)
    model_upper_ci <- setNames(accumulator$upper_ci, accumulator$labels)

    # --- NEW: Format combo outputs if in combo mode ---
    if (mode == "country") {
        model_gradients <- format_combo_results(model_gradients, "Gradient")
        model_lower_ci  <- format_combo_results(model_lower_ci, "Lower_CI")
        model_upper_ci  <- format_combo_results(model_upper_ci, "Upper_CI")
        
        # The bootstraps object is already a dataframe, so we just separate the column directly
        accumulator$bootstraps <- accumulator$bootstraps %>%
            tidyr::separate(Combo_Label, into = c("Pathogen", "Antibiotic"), sep = "___")
    }
    # --------------------------------------------------

    write_random_effects_outputs(
        gradient_prefix = gradient_prefix,
        lower_prefix = lower_prefix,
        upper_prefix = upper_prefix,
        bootstrap_prefix = bootstrap_prefix,
        output_tag = output_tag,
        gradients = model_gradients,
        lower_ci = model_lower_ci,
        upper_ci = model_upper_ci,
        bootstraps = accumulator$bootstraps,
        r_squareds = accumulator$r_squareds,
        white_tests = accumulator$white_tests,
        rmses = accumulator$rmses,
        aics = accumulator$aics,
        bics = accumulator$bics,
        vif_list = accumulator$vif_list
    )
}

resolve_model_jobs <- function(scenario) {
    if (scenario == "main") {
        return(list(list(
            income = "all",
            summed_data_path = "summed_data_new.csv",
            merged_sums_path = "summed_data_sums_new.csv",
            class_output_tag = "all",
            pathogen_output_tag = "main",
            analysis_intent = "main_publication",
            apply_lagged_response = FALSE,
            apply_lagged_consumption = FALSE,
            data_source = "summed",
            output_prefix = "database"
        )))
    }
    
    if (scenario == "main_finer") {
        return(list(list(
            income = "all",
            summed_data_path = "finer_data_new.csv",
            merged_sums_path = "finer_data_sums_new.csv",
            antibiotic_col = "Antibiotic",
            consumption_class_col = "ATC.Class",
            class_output_tag = "all",
            pathogen_output_tag = "main_finer",
            analysis_intent = "main_publication",
            apply_lagged_response = FALSE,
            apply_lagged_consumption = FALSE,
            data_source = "summed",
            output_prefix = "database"
        )))
    }

    if (scenario == "main_2000") {
        return(list(list(
            income = "all",
            summed_data_path = "summed_data_new.csv",
            merged_sums_path = "summed_data_sums_new.csv",
            class_output_tag = "all_2000",
            pathogen_output_tag = "main_2000",
            analysis_intent = "main_publication",
            apply_lagged_response = FALSE,
            apply_lagged_consumption = FALSE,
            data_source = "summed",
            output_prefix = "database",
            reference_year = 2000
        )))
    }

    if (scenario == "main_binomial") {
        return(list(list(
            income = "all",
            summed_data_path = "summed_data_new.csv",
            merged_sums_path = "summed_data_sums_new.csv",
            class_output_tag = "all_binomial",
            pathogen_output_tag = "main_binomial",
            analysis_intent = "main_publication",
            apply_lagged_response = FALSE,
            apply_lagged_consumption = FALSE,
            data_source = "summed",
            output_prefix = "database",
            model_family = "binomial"
        )))
    }

    if (scenario == "main_ppml") {
        return(list(list(
            income = "all",
            summed_data_path = "summed_data_new.csv",
            merged_sums_path = "summed_data_sums_new.csv",
            class_output_tag = "all_ppml",
            pathogen_output_tag = "main_ppml",
            analysis_intent = "main_publication",
            apply_lagged_response = FALSE,
            apply_lagged_consumption = FALSE,
            data_source = "summed",
            output_prefix = "database",
            model_family = "ppml"
        )))
    }

    if (scenario == "main_ppml") {
        return(list(list(
            income = "all",
            summed_data_path = "summed_data_new.csv",
            merged_sums_path = "summed_data_sums_new.csv",
            class_output_tag = "all_ppml",
            pathogen_output_tag = "main_ppml",
            analysis_intent = "main_publication",
            apply_lagged_response = FALSE,
            apply_lagged_consumption = FALSE,
            data_source = "summed",
            output_prefix = "database",
            model_family = "ppml"
        )))
    }

    if (scenario == "extra_pcs") {
        return(list(list(
            income = "all",
            summed_data_path = "summed_data_new.csv",
            merged_sums_path = "summed_data_sums_new.csv",
            class_output_tag = "all_extra_pcs",
            pathogen_output_tag = "extra_pcs",
            analysis_intent = "main_publication",
            apply_lagged_response = FALSE,
            apply_lagged_consumption = FALSE,
            data_source = "summed",
            output_prefix = "database",
            extra_pcs = TRUE
        )))
    }

    if (scenario == "hic") {
        return(list(list(
            income = "HIC",
            summed_data_path = "summed_data_new.csv",
            merged_sums_path = "summed_data_sums_new.csv",
            class_output_tag = "HIC",
            pathogen_output_tag = "HIC",
            analysis_intent = "main_publication",
            apply_lagged_response = FALSE,
            apply_lagged_consumption = FALSE,
            data_source = "summed",
            output_prefix = "database"
        )))
    }

    if (scenario == "hic_ppml") {
        return(list(list(
            income = "HIC",
            summed_data_path = "summed_data_new.csv",
            merged_sums_path = "summed_data_sums_new.csv",
            class_output_tag = "HIC_ppml",
            pathogen_output_tag = "HIC_ppml",
            analysis_intent = "main_publication",
            apply_lagged_response = FALSE,
            apply_lagged_consumption = FALSE,
            data_source = "summed",
            output_prefix = "database",
            model_family = "ppml"
        )))
    }

    if (scenario == "lmic") {
        return(list(list(
            income = "LMIC",
            summed_data_path = "summed_data_new.csv",
            merged_sums_path = "summed_data_sums_new.csv",
            class_output_tag = "LMIC",
            pathogen_output_tag = "LMIC",
            analysis_intent = "main_publication",
            apply_lagged_response = FALSE,
            apply_lagged_consumption = FALSE,
            data_source = "summed",
            output_prefix = "database"
        )))
    }

    if (scenario == "lmic_ppml") {
        return(list(list(
            income = "LMIC",
            summed_data_path = "summed_data_new.csv",
            merged_sums_path = "summed_data_sums_new.csv",
            class_output_tag = "LMIC_ppml",
            pathogen_output_tag = "LMIC_ppml",
            analysis_intent = "main_publication",
            apply_lagged_response = FALSE,
            apply_lagged_consumption = FALSE,
            data_source = "summed",
            output_prefix = "database",
            model_family = "ppml"
        )))
    }

    if (scenario == "raw_iqvia") {
        return(list(
            list(
                income = "all",
                summed_data_path = "summed_data_new_IQVIA.csv",
                merged_sums_path = "summed_data_sums_new_IQVIA.csv",
                class_output_tag = "all_IQVIA",
                pathogen_output_tag = "IQVIA",
                analysis_intent = "main_publication",
                apply_lagged_response = FALSE,
                apply_lagged_consumption = FALSE,
                data_source = "summed",
                output_prefix = "database"
            ),
            list(
                income = "all",
                summed_data_path = "summed_data_new_IQVIAextrapolation.csv",
                merged_sums_path = "summed_data_sums_new_IQVIAextrapolation.csv",
                class_output_tag = "IQVIAextrapolation_all",
                pathogen_output_tag = "IQVIAextrapolation",
                analysis_intent = "main_publication",
                apply_lagged_response = FALSE,
                apply_lagged_consumption = FALSE,
                data_source = "summed",
                output_prefix = "database"
            )
        ))
    }

    if (scenario == "hospital_nagorsen") {
        return(list(list(
            income = "all",
            summed_data_path = "",
            merged_sums_path = "",
            class_output_tag = "hospital_to_all_filtered",
            pathogen_output_tag = "hospital_to_all_filtered",
            analysis_intent = "supplementary_publication",
            apply_lagged_response = FALSE,
            apply_lagged_consumption = FALSE,
            data_source = "nagorsen",
            output_prefix = "Nagorsen",
                min_entries_per_combo = 20,
            prepared_data_path = "summed_data_Nagorsen_hospital_to_all_filtered.csv"
        )))
    }

    if (scenario == "exploratory_lagged") {
        # Fetch the custom lag value first
        custom_lag <- get_integer_env("AMR_LAG_N", default = 1)
        
        # Define the tag suffix: e.g., "lagged" for 1, or "lagged_3y" for 3
        tag_suffix <- if (custom_lag == 1) "lagged" else paste0("lagged_", custom_lag, "y")

        return(list(list(
            income = "all",
            summed_data_path = "summed_data_new.csv",
            merged_sums_path = "summed_data_sums_new.csv",
            class_output_tag = paste0("all_", tag_suffix),
            pathogen_output_tag = tag_suffix,
            analysis_intent = "exploratory_only",
            apply_lagged_response = TRUE,
            apply_lagged_consumption = FALSE,
            lag_n = custom_lag,
            data_source = "summed",
            output_prefix = "database"
        )))
    }

    if (scenario == "consumption_lagged") {
        # Fetch the custom lag value first
        custom_lag <- get_integer_env("AMR_LAG_N", default = 1)
        
        # Define the tag suffix: e.g., "lagged" for 1, or "lagged_3y" for 3
        tag_suffix <- if (custom_lag == 1) "clagged" else paste0("clagged_", custom_lag, "y")

        return(list(list(
            income = "all",
            summed_data_path = "summed_data_new.csv",
            merged_sums_path = "summed_data_sums_new.csv",
            class_output_tag = paste0("all_", tag_suffix),
            pathogen_output_tag = tag_suffix,
            analysis_intent = "exploratory_only",
            apply_lagged_response = FALSE,
            apply_lagged_consumption = TRUE,
            lag_n = custom_lag,
            data_source = "summed",
            output_prefix = "database"
        )))
    }

    if (scenario == "consumption_lagged_ppml") {
        custom_lag <- get_integer_env("AMR_LAG_N", default = 1)
        tag_suffix <- if (custom_lag == 1) "clagged_ppml" else paste0("clagged_ppml_", custom_lag, "y")

        return(list(list(
            income = "all",
            summed_data_path = "summed_data_new.csv",
            merged_sums_path = "summed_data_sums_new.csv",
            class_output_tag = paste0("all_", tag_suffix),
            pathogen_output_tag = tag_suffix,
            analysis_intent = "exploratory_only",
            apply_lagged_response = FALSE,
            apply_lagged_consumption = TRUE,
            lag_n = custom_lag,
            data_source = "summed",
            output_prefix = "database",
            model_family = "ppml"
        )))
    }

    # Permutation scenario: shuffle Consumption for one antibiotic class.
    # Set AMR_PERMUTATION_CLASS to one of J01A, J01C, J01D, J01E, J01F, J01G, J01M.
    # Only fit_combined_pathogen_drug_lm is run (no class/random-effects models).
    # Output bootstrap CSV: database_gradients_bootstraps_pathogen_ATC3_PCA_canonical_weighted_permutation{AB}.csv
    if (scenario == "permutation") {
        perm_class <- Sys.getenv("AMR_PERMUTATION_CLASS", unset = "")
        if (identical(perm_class, "")) {
            stop("[permutation] AMR_PERMUTATION_CLASS environment variable must be set (e.g. J01A)", call. = FALSE)
        }
        valid_classes <- c("J01A", "J01C", "J01D", "J01E", "J01F", "J01G", "J01M")
        if (!perm_class %in% valid_classes) {
            stop(sprintf("[permutation] AMR_PERMUTATION_CLASS='%s' is not a recognised class. Valid: %s",
                         perm_class, paste(valid_classes, collapse = ", ")), call. = FALSE)
        }
        return(list(list(
            income = "all",
            summed_data_path = "summed_data_new.csv",
            merged_sums_path = "summed_data_sums_new.csv",
            class_output_tag = NA_character_,            # not used for permutation
            pathogen_output_tag = paste0("permutation", perm_class),
            analysis_intent = "permutation",
            apply_lagged_response = FALSE,
            apply_lagged_consumption = FALSE,
            data_source = "summed",
            output_prefix = "database",
            permutation_class = perm_class
        )))
    }

    # PPML permutation scenario: same permutation workflow as above, but fit
    # with PPML models and write distinct output tags.
    if (scenario == "permutations_ppml") {
        perm_class <- Sys.getenv("AMR_PERMUTATION_CLASS", unset = "")
        if (identical(perm_class, "")) {
            stop("[permutations_ppml] AMR_PERMUTATION_CLASS environment variable must be set (e.g. J01A)", call. = FALSE)
        }
        valid_classes <- c("J01A", "J01C", "J01D", "J01E", "J01F", "J01G", "J01M")
        if (!perm_class %in% valid_classes) {
            stop(sprintf("[permutations_ppml] AMR_PERMUTATION_CLASS='%s' is not a recognised class. Valid: %s",
                         perm_class, paste(valid_classes, collapse = ", ")), call. = FALSE)
        }
        return(list(list(
            income = "all",
            summed_data_path = "summed_data_new.csv",
            merged_sums_path = "summed_data_sums_new.csv",
            class_output_tag = NA_character_,            # not used for permutation
            pathogen_output_tag = paste0("permutations_ppml", perm_class),
            analysis_intent = "permutation",
            apply_lagged_response = FALSE,
            apply_lagged_consumption = FALSE,
            data_source = "summed",
            output_prefix = "database",
            model_family = "ppml",
            permutation_class = perm_class
        )))
    }
    
    if (scenario == "mi") {
        return(list(list(
            income = "all",
            summed_data_path = "summed_data_new.csv",
            merged_sums_path = "summed_data_sums_new.csv",
            class_output_tag = "all_mi",
            pathogen_output_tag = "mi",
            analysis_intent = "main_mi",
            apply_lagged_response = FALSE,
            apply_lagged_consumption = FALSE,
            data_source = "summed",
            output_prefix = "database"
        )))
    }

    list()
}

run_class_model_pipeline <- function(job) {
    runtime_options <- get_runtime_options()
    model_family <- job_or_default(job, "model_family", "gaussian")

    if (!model_family %in% c("gaussian", "binomial", "ppml")) {
        stop("Unsupported model_family: ", model_family, call. = FALSE)
    }

    # Permutation outputs are consumed as bootstrap distributions for Wilcoxon tests.
    # Running permutation in smoke mode (boot_nsim=0) silently degrades to point estimates.
    if (identical(job$analysis_intent, "permutation") && runtime_options$boot_nsim <= 0) {
        stop(
            "[permutation] boot_nsim=0 (likely smoke mode enabled via AMR_DEV_SMOKE/AMR_SMOKE). ",
            "Permutation Wilcoxon inputs require bootstrap distributions. Re-run with AMR_DEV_SMOKE=0 AMR_SMOKE=0.",
            call. = FALSE
        )
    }

    set.seed(runtime_options$random_seed)

    inputs <- if (identical(job$data_source, "nagorsen")) {
        load_nagorsen_model_inputs(
            prepared_data_path = job_or_default(job, "prepared_data_path", "summed_data_Nagorsen_hospital_to_all_filtered.csv"),
            min_entries_per_combo = job_or_default(job, "min_entries_per_combo", 20)
        )
    } else {
        load_model_inputs(
            summed_data_path = job$summed_data_path,
            merged_sums_path = job$merged_sums_path,
            antibiotic_col = job_or_default(job, "antibiotic_col", "ATC.Class"),
            consumption_class_col = job_or_default(job, "consumption_class_col", "ATC.Class")
        )
    }
    data <- select_income_slice(inputs, job$income)

    ref_year <- resolve_reference_year(scenario = get_amr_scenario(), job = job)
    global_consumption <- build_global_consumption_reference(year = ref_year)
    data <- scale_and_log_transform(
        data,
        global_consumption,
        model_family = model_family
    )

    # if ref_year != 2018, add it to output tags
    if (ref_year != "2018") {
        job$class_output_tag <- paste0(job$class_output_tag, "_", ref_year)
        job$pathogen_output_tag <- paste0(job$pathogen_output_tag, "_", ref_year)
    }

    if (isTRUE(job$apply_lagged_response)) {
        # Retrieve the custom lag value, defaulting to 1
        custom_lag <- job_or_default(job, "lag_n", 1) 
        
        data <- data %>%
            group_by(ISO3, Antibiotic, Pathogen) %>%
            arrange(Year) %>%
            mutate(Resistance = (Resistance - lag(Resistance, n = custom_lag)) / lag(Resistance, n = custom_lag)) %>%
            ungroup() %>%
            filter(is.finite(Resistance))
    }

    if (isTRUE(job$apply_lagged_consumption)) {
        # Retrieve the custom lag value, defaulting to 1
        custom_lag <- job_or_default(job, "lag_n", 1) 
        
        data <- data %>%
            group_by(ISO3, Antibiotic, Pathogen) %>%
            arrange(Year) %>%
            mutate(Consumption = lag(Consumption, n = custom_lag)) %>%
            ungroup() %>%
            filter(is.finite(Consumption))
    }

    # Permutation: shuffle Consumption within each country-year for one antibiotic class only.
    # This breaks the real signal for that class while preserving covariate structure.
    if (!is.null(job$permutation_class) && !is.na(job$permutation_class)) {
        perm_class <- job$permutation_class
        message("[permutation] Shuffling Consumption for class: ", perm_class)
        perm_rows <- data$Antibiotic == perm_class
        data$Consumption[perm_rows] <- sample(data$Consumption[perm_rows])
    }

    data_ <- limit_for_smoke_mode(data, runtime_options)

    log_info(
        "[ddd-linear-model] runtime_options: smoke_mode=", runtime_options$smoke_mode,
        ", boot_nsim=", runtime_options$boot_nsim,
        ", random_seed=", runtime_options$random_seed,
        ", max_classes=", runtime_options$smoke_max_classes,
        ", max_pathogens=", runtime_options$smoke_max_pathogens,
        ", max_pairs=", runtime_options$smoke_max_pairs,
        ", max_rows_per_pair=", runtime_options$smoke_max_rows_per_pair,
        verbose = runtime_options$verbose
    )

    is_extra_pcs <- isTRUE(job$extra_pcs) # Safely handles NULLs

    fit_combined_pathogen_drug_lm(
        data_,
        output_tag = job$pathogen_output_tag,
        runtime_options = runtime_options,
        output_prefix = job$output_prefix,
        model_family = model_family,
        extra_pcs = is_extra_pcs
    )

    allow_fallback <- TRUE

    # Detect how many cores your machine has, minus 1 for safety
    n_cores <- parallel::detectCores() - 1

    # Tell R to resolve tasks in parallel using independent R sessions
    plan(multisession, workers = n_cores)

    # Class and random-effects models are not needed for permutation runs.
    if (!identical(job$analysis_intent, "permutation")) {
        fit_random_effects_models(
            data_,
            output_tag = job$class_output_tag,
            runtime_options = runtime_options,
            mode = "class",
            output_prefix = job$output_prefix,
            allow_fallback = allow_fallback,
            model_family = model_family,
            extra_pcs = is_extra_pcs
        )
        fit_random_effects_models(
            data_,
            output_tag = job$class_output_tag,
            runtime_options = runtime_options,
            mode = "pathogen",
            output_prefix = job$output_prefix,
            allow_fallback = allow_fallback,
            model_family = model_family,
            extra_pcs = is_extra_pcs
        )
        # fit_random_effects_models(
        #     data_,
        #     output_tag = job$class_output_tag,
        #     runtime_options = runtime_options,
        #     mode = "country",
        #     output_prefix = job$output_prefix,
        #     allow_fallback = allow_fallback,
        #     model_family = model_family,
        #     extra_pcs = is_extra_pcs
        # )
    }
}

run_mi_class_model_pipeline <- function(job) {
    runtime_options <- get_runtime_options()

    runtime_options$boot_nsim <- 1000 
    model_family <- job_or_default(job, "model_family", "gaussian")
    
    set.seed(runtime_options$random_seed)
    global_consumption <- build_global_consumption_reference()
    
    # NEW: List to hold all 20 processed datasets
    mi_datasets <- list()
    
    for (i in 1:20) {
        imputed_file <- paste0("imputed_data_", i, ".csv")
        if (!file.exists(imputed_file)) next
        
        # Temporarily overwrite job path for the loader
        job$summed_data_path <- imputed_file 
        
        inputs <- load_model_inputs(
            summed_data_path = job$summed_data_path,
            merged_sums_path = job$merged_sums_path
        )
        data <- select_income_slice(inputs, job$income)
        
        data <- scale_and_log_transform(
            data,
            global_consumption,
            model_family = model_family
        )

        # if global option MI_EQUAL=1, then set Weight to 1 for all rows
        if (get_integer_env("MI_EQUAL", default = 0) == 1) {
            message("[ddd-linear-model] MI_EQUAL=1: Setting all weights to 1 for imputed dataset ", i)
            data$Weight <- 1
        }
        
        data_ <- limit_for_smoke_mode(data, runtime_options)
        mi_datasets[[i]] <- data_
    }
    
    if (length(mi_datasets) == 0) {
        stop("No imputed datasets found matching 'imputed_data1.csv', etc.")
    }

    is_extra_pcs <- isTRUE(job$extra_pcs)
    
    # Pass the entire list of datasets to the fitting function
    fit_combined_pathogen_drug_mi_lm(
        mi_datasets,
        output_tag = paste0(job$pathogen_output_tag, "_MI_pooled"),
        runtime_options = runtime_options,
        output_prefix = job$output_prefix,
        model_family = model_family,
        extra_pcs = is_extra_pcs
    )
}

require(lme4)
scenario <- get_amr_scenario()
jobs <- resolve_model_jobs(scenario)

if (length(jobs) == 0) {
    message("[ddd-linear-model] No model jobs configured for scenario: ", scenario)
} else {
    for (job in jobs) {
        message(
            "[ddd-linear-model] Running job with data=", job$summed_data_path,
            ", income=", job$income,
            ", analysis_intent=", job$analysis_intent,
            ", apply_lagged_response=", job$apply_lagged_response,
            ", apply_lagged_consumption=", job$apply_lagged_consumption,
            ", class_output_tag=", job$class_output_tag,
            ", pathogen_output_tag=", job$pathogen_output_tag
        )
        if (identical(job$analysis_intent, "main_mi")) {
            run_mi_class_model_pipeline(job)
        } else {
            run_class_model_pipeline(job)
        }
    }
}