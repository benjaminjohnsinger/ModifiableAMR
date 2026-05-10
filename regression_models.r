# ## Linear regression on database DDD data
# ## BJS March 2025
library(tidyverse)
library(dplyr)
library(glmnet)
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
    prepared_data_path = "merged_data_Nagorsen_hospital_to_all_filtered.csv",
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
    merged_data_path = "merged_data_new.csv",
    merged_sums_path = "merged_data_sums_new.csv",
    min_entries_per_combo = 10,
    min_isolates_per_combo = 100
) {
    data <- read.csv(merged_data_path)

    runtime_options <- get_runtime_options()
    log_info("[ddd-linear-model] Input file: ", merged_data_path, verbose = runtime_options$verbose)
    log_info("[ddd-linear-model] Unique ISO3 count: ", length(table(data$ISO3)), verbose = runtime_options$verbose)
    log_info("[ddd-linear-model] Max rows for one ISO3: ", max(table(data$ISO3)), verbose = runtime_options$verbose)
    log_info("[ddd-linear-model] Years included: ", paste(sort(unique(data$Year)), collapse = ", "), verbose = runtime_options$verbose)

    data <- data %>%
        rename(
            Antibiotic = ATC.Class,
            Consumption = Antibiotic.Consumption,
            Resistance = Percent.Resistant.Isolates,
            Pathogen = Pathogen,
            Weight = Total.Isolates
        )

    country_covariates <- c("PC1", "PC2", "PC3", "GDP", "Year")
    data$GDP <- data$GDP / mean(data$GDP, na.rm = TRUE)

    data <- data[
        !is.na(data$Consumption) & !is.na(data$Resistance) &
            !is.na(data$Pathogen) & !is.na(data$Antibiotic) & !is.na(data$Weight),
    ]

    pathogen_drug_counts <- table(paste(data$Pathogen, data$Antibiotic))
    pathogen_drug_to_remove <- names(pathogen_drug_counts[pathogen_drug_counts <= min_entries_per_combo])
    data <- data[!paste(data$Pathogen, data$Antibiotic) %in% pathogen_drug_to_remove, ]

    merged_data_sums <- read.csv(merged_sums_path) %>%
        rename(Antibiotic = ATC.Class)

    pathogen_drug_counts <- merged_data_sums %>%
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

    countries_per_pathogen <- data %>%
        group_by(Pathogen) %>%
        summarise(Countries = paste(unique(ISO3), collapse = ", "), .groups = "drop")

    for (i in seq_len(nrow(countries_per_pathogen))) {
        log_info(paste0(countries_per_pathogen$Pathogen[i], ": ", countries_per_pathogen$Countries[i]), verbose = runtime_options$verbose)
    }

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
        select(Consumption, Resistance, Pathogen, Antibiotic, Weight, ISO3, 
               all_of(country_covariates), ends_with(".Consumption"))

    list(data = data, data_HIC = data_HIC, data_LMIC = data_LMIC)
}
build_global_consumption_reference <- function(consumption_path = "antibiotic_consumption_by_ATC3.csv") {
    consumption <- read.csv(consumption_path)
    global_consumption <- consumption[consumption$Location == "Global", ]
    global_consumption <- global_consumption[global_consumption$Year == "2018", ]
    global_consumption <- select(global_consumption, ATC.level.3.class, Antibiotic.consumption..DDD.1.000.day.)
    global_consumption <- global_consumption %>%
        rename(
            Antibiotic = ATC.level.3.class,
            Global.Consumption = Antibiotic.consumption..DDD.1.000.day.
        )
    global_consumption$Antibiotic <- sub("-.*", "", global_consumption$Antibiotic)
    global_consumption
}

scale_and_log_transform <- function(df, global_consumption) {
    df <- df %>%
        left_join(global_consumption, by = c("Antibiotic")) %>%
        mutate(Consumption = Consumption / Global.Consumption) %>%
        select(-Global.Consumption)

    df <- df %>%
        filter(!is.na(Consumption) & !is.na(Resistance) & !is.na(Weight))

    df$Consumption <- log(df$Consumption + 1)
    df$Resistance <- log(df$Resistance + 1)
    df$Weight <- df$Weight / max(df$Weight, na.rm = TRUE)
    df
}

get_fixed_effects_formula <- function() {
    Resistance ~ Consumption + PC1 + PC2 + PC3 + GDP + Year
}

fit_weighted_lm <- function(data_subset) {
    lm(
        formula = get_fixed_effects_formula(),
        data = data_subset,
        weights = Weight
    )
}

fit_random_lmer <- function(data_subset, random_effect_var) {
    formula_str <- paste0(
        "Resistance ~ Consumption + (Consumption||", random_effect_var,
        ") + PC1 + PC2 + PC3 + GDP + Year"
    )
    lmer(
        formula = as.formula(formula_str),
        data = data_subset,
        weights = Weight
    )
}

build_output_path <- function(prefix, output_tag) {
    paste0("Outputs/", prefix, "_", output_tag, ".csv")
}

write_random_effects_outputs <- function(
    gradient_prefix,
    lower_prefix,
    upper_prefix,
    bootstrap_prefix,
    output_tag,
    gradients,
    lower_ci,
    upper_ci,
    bootstraps
) {
    write.csv(gradients, build_output_path(gradient_prefix, output_tag), row.names = TRUE)
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
        bootstraps = data.frame()
    )
}

append_random_effects_result <- function(
    accumulator,
    label,
    gradient,
    intercept,
    lower_ci,
    upper_ci,
    bootstrap_df
) {
    accumulator$labels <- c(accumulator$labels, label)
    accumulator$gradients <- c(accumulator$gradients, gradient)
    accumulator$intercepts <- c(accumulator$intercepts, intercept)
    accumulator$lower_ci <- c(accumulator$lower_ci, lower_ci)
    accumulator$upper_ci <- c(accumulator$upper_ci, upper_ci)
    accumulator$bootstraps <- rbind(accumulator$bootstraps, bootstrap_df)
    accumulator
}

extract_lm_consumption_summary <- function(model, boot_nsim) {
    gradient <- summary(model)$coefficients["Consumption", 1]
    intercept <- summary(model)$coefficients["(Intercept)", 1]
    intervals <- confint(model)
    bootstrap_values <- if (boot_nsim > 0) {
        tryCatch(
            {
                bs <- car::Boot(model, R = boot_nsim)
                bs$t[, "Consumption"]
            },
            error = function(e) {
                warning("car::Boot failed for LM Consumption term; using point estimate fallback. Error: ", conditionMessage(e))
                gradient
            }
        )
    } else {
        gradient
    }

    list(
        gradient = gradient,
        intercept = intercept,
        lower_ci = intervals["Consumption", 1],
        upper_ci = intervals["Consumption", 2],
        bootstrap_values = bootstrap_values
    )
}

extract_lmer_consumption_summary <- function(model, boot_nsim) {
    gradient <- summary(model)$coefficients["Consumption", 1]
    intercept <- summary(model)$coefficients["(Intercept)", 1]

    if (boot_nsim > 0) {
        boot_result <- tryCatch(
            {
                bs <- bootMer(model, FUN = function(x) fixef(x)["Consumption"], nsim = boot_nsim)
                intervals <- confint(bs)
                list(
                    lower_ci = intervals["Consumption", "2.5 %"],
                    upper_ci = intervals["Consumption", "97.5 %"],
                    bootstrap_values = as.vector(bs$t)
                )
            },
            error = function(e) {
                warning("bootMer failed for LMER Consumption term; using Wald CI fallback. Error: ", conditionMessage(e))
                intervals <- confint(model, parm = "Consumption", method = "Wald")
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
        intervals <- confint(model, parm = "Consumption", method = "Wald")
        lower_ci <- intervals[1, 1]
        upper_ci <- intervals[1, 2]
        bootstrap_values <- gradient
    }

    list(
        gradient = gradient,
        intercept = intercept,
        lower_ci = lower_ci,
        upper_ci = upper_ci,
        bootstrap_values = bootstrap_values
    )
}

build_bootstrap_df <- function(label_name, label_value, bootstrap_values) {
    outdf <- data.frame(Gradient = bootstrap_values)
    outdf[[label_name]] <- label_value
    outdf[, c(label_name, "Gradient"), drop = FALSE]
}

fit_combined_pathogen_drug_lm <- function(data_, output_tag = "lagged", runtime_options = get_runtime_options(), output_prefix = "database") {
    gradients <- c()
    conf_intervals <- c()
    pathogens <- c()
    abs <- c()
    bootstraps <- data.frame()
    r_squareds <- c()
    variation_explained_list <- list()

    antibiotics_to_fit <- sort(unique(data_$Antibiotic))
    pathogens_to_fit <- sort(unique(data_$Pathogen))

    if (isTRUE(runtime_options$smoke_mode)) {
        antibiotics_to_fit <- head(antibiotics_to_fit, runtime_options$smoke_max_classes)
        pathogens_to_fit <- head(pathogens_to_fit, runtime_options$smoke_max_pathogens)
    }

    for (antibiotic in antibiotics_to_fit) {
        for (pathogen in pathogens_to_fit) {
            data_subset <- data_[data_$Pathogen == pathogen & data_$Antibiotic == antibiotic, ]
            vars_needed <- c("Resistance", "Consumption", "PC1", "PC2", "PC3", "GDP", "Year", "Weight")
            data_subset <- data_subset[complete.cases(data_subset[, vars_needed]), ]
            data_subset <- data_subset[!is.infinite(data_subset$Consumption), ]
            data_subset <- data_subset[!is.infinite(data_subset$Resistance), ]
            if (nrow(data_subset) <= 1) {
                next
            }

            model <- fit_weighted_lm(data_subset)

            # Compute R-squared and per-variable variance explained (leave-one-out)
            r_sq <- summary(model)$r.squared
            ve_vars <- c("Consumption", "PC1", "PC2", "PC3", "GDP", "Year")
            ve <- setNames(rep(NA_real_, length(ve_vars)), ve_vars)
            for (ve_var in ve_vars) {
                other_vars <- ve_vars[ve_vars != ve_var]
                formula_null <- as.formula(paste("Resistance ~", paste(other_vars, collapse = " + ")))
                model_null <- tryCatch(
                    lm(formula_null, data = data_subset, weights = data_subset$Weight),
                    error = function(e) NULL
                )
                if (!is.null(model_null)) {
                    ve[ve_var] <- r_sq - summary(model_null)$r.squared
                }
            }

            if (is.na(coef(model)["Consumption"])) {
                next
            }

            gradients <- c(gradients, coef(model)["Consumption"])
            conf_intervals <- c(conf_intervals, confint(model)["Consumption", ])
            pathogens <- c(pathogens, pathogen)
            abs <- c(abs, antibiotic)
            r_squareds <- c(r_squareds, r_sq)
            variation_explained_list[[length(variation_explained_list) + 1]] <- ve

            if (runtime_options$boot_nsim > 0) {
                # Manual bootstrap: resample rows and refit, avoiding car::Boot scope issues
                boot_values <- tryCatch(
                    {
                        local({
                            boot_df <- data_subset
                            n <- nrow(boot_df)
                            replicate(runtime_options$boot_nsim, {
                                idx <- sample(n, n, replace = TRUE)
                                m_b <- lm(get_fixed_effects_formula(), data = boot_df[idx, ], weights = boot_df$Weight[idx])
                                coef(m_b)["Consumption"]
                            })
                        })
                    },
                    error = function(e) {
                        warning(
                            "Bootstrap failed for combined LM model ",
                            pathogen,
                            " x ",
                            antibiotic,
                            "; using point estimate fallback. Error: ",
                            conditionMessage(e)
                        )
                        coef(model)["Consumption"]
                    }
                )
                outdf <- data.frame(Pathogen = pathogen, Antibiotic = antibiotic, Gradient = boot_values)
            } else {
                outdf <- data.frame(Pathogen = pathogen, Antibiotic = antibiotic, Gradient = coef(model)["Consumption"])
            }
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
    results <- data.frame(
        Antibiotic = abs,
        Pathogen = pathogens,
        Response = gradients,
        Lower_CI = conf_intervals[, 1],
        Upper_CI = conf_intervals[, 2],
        R_squared = r_squareds
    )
    results <- cbind(results, variation_explained_df)

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
            
            # Find the index of our target "Consumption" column and set its penalty to 0
            target_idx <- which(colnames(x_matrix) == "Consumption")
            if (length(target_idx) > 0) {
                p_fac[target_idx] <- 0 
            } else {
                warning("Consumption column not found in x_matrix for ", pathogen, " x ", antibiotic, "; proceeding without unpenalized term.")
            }

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

fit_random_effects_models <- function(
    data_,
    output_tag,
    runtime_options = get_runtime_options(),
    mode = c("class", "pathogen"),
    output_prefix = "database",
    allow_fallback = FALSE
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
    } else {
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
        if (n_levels > 2) {
            model <- fit_random_lmer(subset_data, random_effect_var = random_effect_var)
            log_info(capture.output(print(model)), verbose = runtime_options$verbose)

            if (isSingular(model)) {
                if (allow_fallback) {
                    log_info(paste("Model is singular for", singular_msg, label, "- falling back to weighted lm"), verbose = runtime_options$verbose)
                    model <- fit_weighted_lm(subset_data)
                    summary_stats <- extract_lm_consumption_summary(model, runtime_options$boot_nsim)
                } else {
                    log_info(paste("Model is singular for", singular_msg, label, "- excluding"), verbose = runtime_options$verbose)
                    next
                }
            } else {
                summary_stats <- extract_lmer_consumption_summary(model, runtime_options$boot_nsim)
            }
        } else {
            if (allow_fallback) {
                log_info(paste("Only", n_levels, "levels for", singular_msg, label, "- falling back to weighted lm"), verbose = runtime_options$verbose)
                model <- fit_weighted_lm(subset_data)
                summary_stats <- extract_lm_consumption_summary(model, runtime_options$boot_nsim)
            } else {
                log_info(paste("Only", n_levels, "levels for", singular_msg, label, "- excluding"), verbose = runtime_options$verbose)
                next
            }
        }
        # Centralized appendage handles both lm and lmer uniformly
        outdf <- build_bootstrap_df(label_var, label, summary_stats$bootstrap_values)
        accumulator <- append_random_effects_result(
            accumulator = accumulator,
            label = label,
            gradient = summary_stats$gradient,
            intercept = summary_stats$intercept,
            lower_ci = summary_stats$lower_ci,
            upper_ci = summary_stats$upper_ci,
            bootstrap_df = outdf
        )
    }

    model_gradients <- setNames(accumulator$gradients, accumulator$labels)
    model_intercepts <- setNames(accumulator$intercepts, accumulator$labels)
    model_lower_ci <- setNames(accumulator$lower_ci, accumulator$labels)
    model_upper_ci <- setNames(accumulator$upper_ci, accumulator$labels)

    write_random_effects_outputs(
        gradient_prefix = gradient_prefix,
        lower_prefix = lower_prefix,
        upper_prefix = upper_prefix,
        bootstrap_prefix = bootstrap_prefix,
        output_tag = output_tag,
        gradients = model_gradients,
        lower_ci = model_lower_ci,
        upper_ci = model_upper_ci,
        bootstraps = accumulator$bootstraps
    )
}

resolve_model_jobs <- function(scenario) {
    if (scenario == "main") {
        return(list(list(
            income = "all",
            merged_data_path = "merged_data_new.csv",
            merged_sums_path = "merged_data_sums_new.csv",
            class_output_tag = "all",
            pathogen_output_tag = "main",
            random_pathogen_output_tag = "all",
            analysis_intent = "main_publication",
            apply_lagged_response = FALSE,
            data_source = "merged",
            output_prefix = "database"
        )))
    }

    if (scenario == "hic") {
        return(list(list(
            income = "HIC",
            merged_data_path = "merged_data_new.csv",
            merged_sums_path = "merged_data_sums_new.csv",
            class_output_tag = "HIC",
            pathogen_output_tag = "HIC",
            random_pathogen_output_tag = "HIC",
            analysis_intent = "main_publication",
            apply_lagged_response = FALSE,
            data_source = "merged",
            output_prefix = "database"
        )))
    }

    if (scenario == "lmic") {
        return(list(list(
            income = "LMIC",
            merged_data_path = "merged_data_new.csv",
            merged_sums_path = "merged_data_sums_new.csv",
            class_output_tag = "LMIC",
            pathogen_output_tag = "LMIC",
            random_pathogen_output_tag = "LMIC",
            analysis_intent = "main_publication",
            apply_lagged_response = FALSE,
            data_source = "merged",
            output_prefix = "database"
        )))
    }

    if (scenario == "raw_iqvia") {
        return(list(
            list(
                income = "all",
                merged_data_path = "merged_data_new_IQVIA.csv",
                merged_sums_path = "merged_data_sums_new_IQVIA.csv",
                class_output_tag = "all_IQVIA",
                pathogen_output_tag = "IQVIA",
                random_pathogen_output_tag = "all_IQVIA",
                analysis_intent = "main_publication",
                apply_lagged_response = FALSE,
                data_source = "merged",
                output_prefix = "database"
            ),
            list(
                income = "all",
                merged_data_path = "merged_data_new_IQVIAextrapolation.csv",
                merged_sums_path = "merged_data_sums_new_IQVIAextrapolation.csv",
                class_output_tag = "IQVIAextrapolation_all",
                pathogen_output_tag = "IQVIAextrapolation",
                random_pathogen_output_tag = "IQVIAextrapolation_all",
                analysis_intent = "main_publication",
                apply_lagged_response = FALSE,
                data_source = "merged",
                output_prefix = "database"
            )
        ))
    }

    if (scenario == "hospital_nagorsen") {
        return(list(list(
            income = "all",
            merged_data_path = "",
            merged_sums_path = "",
            class_output_tag = "hospital_to_all_filtered",
            pathogen_output_tag = "hospital_to_all_filtered",
            random_pathogen_output_tag = "hospital_to_all_filtered",
            analysis_intent = "supplementary_publication",
            apply_lagged_response = FALSE,
            data_source = "nagorsen",
            output_prefix = "Nagorsen",
                min_entries_per_combo = 20,
            prepared_data_path = "merged_data_Nagorsen_hospital_to_all_filtered.csv"
        )))
    }

    if (scenario == "exploratory_lagged") {
        return(list(list(
            income = "all",
            merged_data_path = "merged_data_new.csv",
            merged_sums_path = "merged_data_sums_new.csv",
            class_output_tag = "all_lagged",
            pathogen_output_tag = "lagged",
            random_pathogen_output_tag = "all_lagged",
            analysis_intent = "exploratory_only",
            apply_lagged_response = TRUE,
            data_source = "merged",
            output_prefix = "database"
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
            merged_data_path = "merged_data_new.csv",
            merged_sums_path = "merged_data_sums_new.csv",
            class_output_tag = NA_character_,            # not used for permutation
            pathogen_output_tag = paste0("permutation", perm_class),
            random_pathogen_output_tag = NA_character_,  # not used for permutation
            analysis_intent = "permutation",
            apply_lagged_response = FALSE,
            data_source = "merged",
            output_prefix = "database",
            permutation_class = perm_class
        )))
    }

    list()
}

run_class_model_pipeline <- function(job) {
    runtime_options <- get_runtime_options()
    set.seed(runtime_options$random_seed)

    inputs <- if (identical(job$data_source, "nagorsen")) {
        load_nagorsen_model_inputs(
            prepared_data_path = job_or_default(job, "prepared_data_path", "merged_data_Nagorsen_hospital_to_all_filtered.csv"),
            min_entries_per_combo = job_or_default(job, "min_entries_per_combo", 20)
        )
    } else {
        load_model_inputs(
            merged_data_path = job$merged_data_path,
            merged_sums_path = job$merged_sums_path
        )
    }
    data <- select_income_slice(inputs, job$income)

    global_consumption <- build_global_consumption_reference()
    data <- scale_and_log_transform(data, global_consumption)

    if (isTRUE(job$apply_lagged_response)) {
        # Exploratory-only pathway: lagged response is not part of the main or supplementary publication analyses.
        data <- data %>%
            group_by(ISO3, Antibiotic, Pathogen) %>%
            arrange(Year) %>%
            mutate(Resistance = (Resistance - lag(Resistance)) / lag(Resistance)) %>%
            ungroup() %>%
            filter(is.finite(Resistance))
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

    fit_combined_pathogen_drug_lm(
        data_,
        output_tag = job$pathogen_output_tag,
        runtime_options = runtime_options,
        output_prefix = job$output_prefix
    )

    # Determine fallback based on whether this is an income-stratified job
    allow_fallback <- job$income != "all"

    # Class and random-effects models are not needed for permutation runs.
    if (!identical(job$analysis_intent, "permutation")) {
        fit_random_effects_models(
            data_,
            output_tag = job$class_output_tag,
            runtime_options = runtime_options,
            mode = "class",
            output_prefix = job$output_prefix,
            allow_fallback = allow_fallback
        )
        fit_random_effects_models(
            data_,
            output_tag = job$random_pathogen_output_tag,
            runtime_options = runtime_options,
            mode = "pathogen",
            output_prefix = job$output_prefix,
            allow_fallback = allow_fallback
        )
    }
}
require(lme4)
scenario <- get_amr_scenario()
jobs <- resolve_model_jobs(scenario)

if (length(jobs) == 0) {
    message("[ddd-linear-model] No model jobs configured for scenario: ", scenario)
} else {
    for (job in jobs) {
        message(
            "[ddd-linear-model] Running job with data=", job$merged_data_path,
            ", income=", job$income,
            ", analysis_intent=", job$analysis_intent,
            ", apply_lagged_response=", job$apply_lagged_response,
            ", class_output_tag=", job$class_output_tag,
            ", pathogen_output_tag=", job$pathogen_output_tag,
            ", random_pathogen_output_tag=", job$random_pathogen_output_tag
        )
        run_class_model_pipeline(job)
    }
}

