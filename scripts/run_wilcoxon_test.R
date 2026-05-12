# =============================================================================
# run_wilcoxon_test.R
# Wilcoxon rank-sum test comparing bootstrapped gradients to permuted gradients
# per antibiotic class and pathogen, to assess significance of model associations.
#
# Requires:
#   - Bootstrapped gradient CSV (main scenario output from regression_models.r)
#   - Per-class permuted bootstrap CSVs (from a permutation model run)
#     Pattern: Outputs/database_gradients_bootstraps_pathogen_ATC3_PCA_canonical_weighted_permutation{AB}.csv
#
# Writes:
#   - Outputs/wilcoxon_bootstrapped_vs_permuted_gradients.csv
# =============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(data.table)
})

# Resolve paths from options (set by stage_tables.R) or use defaults
bootstrapped_path <- getOption(
    "amr_wilcoxon_bootstrapped_path",
    "Outputs/database_gradients_bootstraps_pathogen_ATC3_PCA_canonical_weighted_main.csv"
)
permuted_pattern <- getOption(
    "amr_wilcoxon_permuted_pattern",
    "Outputs/database_gradients_bootstraps_pathogen_ATC3_PCA_canonical_weighted_permutation%s.csv"
)
output_path <- getOption(
    "amr_wilcoxon_output_path",
    "Outputs/wilcoxon_bootstrapped_vs_permuted_gradients.csv"
)
permutation_pairing <- getOption(
    "amr_wilcoxon_permutation_pairing",
    "mismatch_only"
)
gradient_path <- getOption(
    "amr_wilcoxon_gradient_path",
    "Outputs/database_gradients_pathogen_ATC3_PCA_canonical_weighted_main.csv"
)
permuted_gradient_pattern <- getOption(
    "amr_wilcoxon_permuted_gradient_pattern",
    "Outputs/database_gradients_pathogen_ATC3_PCA_canonical_weighted_permutation%s.csv"
)

antibiotic_classes <- c("J01A", "J01C", "J01D", "J01E", "J01F", "J01G", "J01M")

.report_gradient_summary <- function(main_gradients, permuted_gradients, wilcoxon_results = NULL) {
    if (is.null(main_gradients) || is.null(permuted_gradients)) {
        return(invisible(NULL))
    }

    required_cols <- c("Response", "Lower_CI", "Upper_CI")
    if (!all(required_cols %in% names(main_gradients)) || !all(required_cols %in% names(permuted_gradients))) {
        message("[wilcoxon] Skipping gradient summary: missing Response/Lower_CI/Upper_CI columns.")
        return(invisible(NULL))
    }

    sig_pos_main <- main_gradients$Response > 0 & main_gradients$Lower_CI > 0
    sig_neg_main <- main_gradients$Response < 0 & main_gradients$Upper_CI < 0
    sig_pos_perm <- permuted_gradients$Response > 0 & permuted_gradients$Lower_CI > 0
    sig_neg_perm <- permuted_gradients$Response < 0 & permuted_gradients$Upper_CI < 0

    message("[wilcoxon] Gradient summary using 95% CIs as a two-sided p<0.05 proxy:")
    message(sprintf(
        "[wilcoxon]   Primary analysis: %d/%d positive (%.1f%%) and %d/%d negative (%.1f%%) coefficients were significant.",
        sum(sig_pos_main), nrow(main_gradients), 100 * mean(sig_pos_main),
        sum(sig_neg_main), nrow(main_gradients), 100 * mean(sig_neg_main)
    ))
    message(sprintf(
        "[wilcoxon]   Permuted analysis: %d/%d positive (%.1f%%) and %d/%d negative (%.1f%%) coefficients were significant.",
        sum(sig_pos_perm), nrow(permuted_gradients), 100 * mean(sig_pos_perm),
        sum(sig_neg_perm), nrow(permuted_gradients), 100 * mean(sig_neg_perm)
    ))

    if (!is.null(wilcoxon_results)) {
        positive_pairs <- main_gradients[sig_pos_main, c("Antibiotic", "Pathogen")]
        if (nrow(positive_pairs) == 0) {
            message("[wilcoxon]   No positive significant primary coefficients were found for the Wilcoxon comparison.")
            return(invisible(NULL))
        }

        matched_results <- merge(
            positive_pairs,
            wilcoxon_results,
            by = c("Antibiotic", "Pathogen"),
            all.x = TRUE,
            sort = FALSE
        )

        lower_perm_count <- sum(matched_results$p_value < 0.05, na.rm = TRUE)
        message(sprintf(
            "[wilcoxon]   Of the %d positive significant primary coefficients, %d (%.1f%%) had significantly lower permuted estimates in the Wilcoxon test (p<0.05).",
            nrow(positive_pairs), lower_perm_count, 100 * lower_perm_count / nrow(positive_pairs)
        ))
    }
    invisible(NULL)
}

# Guard helper that works both when run as a script and when source()'d
.wilcoxon_early_exit <- function(msg) {
    message(msg)
    if (sys.nframe() <= 1L) {
        quit(save = "no", status = 0)
    } else {
        stop(".wilcoxon_skip", call. = FALSE)
    }
}

# --- Guard: bootstrapped file must exist ---
if (!file.exists(bootstrapped_path)) {
    .wilcoxon_early_exit(paste0(
        "[wilcoxon] Skipping: bootstrapped gradient file not found:\n  ", bootstrapped_path,
        "\n[wilcoxon]   Run the main regression stage first."
    ))
}

# --- Guard: check at least one permuted file exists ---
permuted_files   <- sprintf(permuted_pattern, antibiotic_classes)
missing_permuted <- permuted_files[!file.exists(permuted_files)]
if (length(missing_permuted) == length(permuted_files)) {
    .wilcoxon_early_exit(paste0(
        "[wilcoxon] Skipping: no permuted bootstrap files found.\n",
        "[wilcoxon]   Run the permutation model to generate these files:\n  ",
        paste(permuted_files, collapse = "\n  ")
    ))
}
if (length(missing_permuted) > 0) {
    message("[wilcoxon] Warning: ", length(missing_permuted),
            " permuted file(s) missing — skipping those classes:")
    message("[wilcoxon]   ", paste(missing_permuted, collapse = "\n  "))
}

# --- Load bootstrapped gradients ---
message("[wilcoxon] Loading bootstrapped gradients: ", bootstrapped_path)
bootstraps <- fread(bootstrapped_path)

# --- Load point-estimate gradient summaries for manuscript-style summary messages ---
if (file.exists(gradient_path)) {
    message("[wilcoxon] Loading gradient summary file: ", gradient_path)
    main_gradients <- fread(gradient_path)
    permuted_gradients <- rbindlist(
        lapply(antibiotic_classes, function(ab) {
            perm_path <- sprintf(permuted_gradient_pattern, ab)
            if (!file.exists(perm_path)) {
                return(NULL)
            }
            perm_dt <- fread(perm_path)
            perm_dt$Antibiotic_permuted <- ab
            perm_dt
        }),
        fill = TRUE
    )
}

# --- Load permuted bootstraps ---
bootstraps_permuted <- data.frame()
for (ab in antibiotic_classes) {
    perm_path <- sprintf(permuted_pattern, ab)
    if (!file.exists(perm_path)) next
    message("[wilcoxon] Loading permuted gradients for ", ab)
    ab_df <- fread(perm_path)
    ab_df$Antibiotic_permuted <- ab
    bootstraps_permuted <- rbind(bootstraps_permuted, ab_df)
}

if (nrow(bootstraps_permuted) == 0) {
    message("[wilcoxon] No permuted data loaded. Aborting.")
    quit(save = "no", status = 0)
}

# --- Run Wilcoxon test per antibiotic class x pathogen ---
# Column holding the gradient value in the bootstrap CSV is "Gradient"
gradient_col <- if ("Gradient" %in% colnames(bootstraps)) "Gradient" else "Gradient.Consumption"
perm_gradient_col <- if ("Gradient" %in% colnames(bootstraps_permuted)) "Gradient" else "Gradient.Consumption"

valid_pairings <- c("mismatch_only", "legacy_pool", "matched_only")
if (!permutation_pairing %in% valid_pairings) {
    stop(
        "[wilcoxon] Invalid amr_wilcoxon_permutation_pairing='", permutation_pairing,
        "'. Valid values: ", paste(valid_pairings, collapse = ", "),
        call. = FALSE
    )
}

pathogens <- unique(bootstraps$Pathogen)
available_classes <- unique(bootstraps_permuted$Antibiotic_permuted)

results <- data.frame(
    Antibiotic         = character(),
    Pathogen           = character(),
    p_value            = numeric(),
    wilcoxon_statistic = numeric(),
    stringsAsFactors   = FALSE
)

for (ab in available_classes) {
    for (pathogen in pathogens) {
        bootstrapped_values <- bootstraps[bootstraps$Antibiotic == ab &
                                          bootstraps$Pathogen == pathogen, ][[gradient_col]]
        perm_idx <- bootstraps_permuted$Antibiotic_permuted == ab &
                    bootstraps_permuted$Pathogen == pathogen
        if (permutation_pairing == "mismatch_only") {
            perm_idx <- perm_idx & bootstraps_permuted$Antibiotic != ab
        } else if (permutation_pairing == "matched_only") {
            perm_idx <- perm_idx & bootstraps_permuted$Antibiotic == ab
        }
        permuted_values <- bootstraps_permuted[perm_idx, ][[perm_gradient_col]]

        if (length(bootstrapped_values) == 0 || length(permuted_values) == 0) next

        test_result <- wilcox.test(bootstrapped_values, permuted_values, alternative = "greater")
        results <- rbind(results, data.frame(
            Antibiotic         = ab,
            Pathogen           = pathogen,
            p_value            = test_result$p.value,
            wilcoxon_statistic = as.numeric(test_result$statistic),
            stringsAsFactors   = FALSE
        ))
    }
}

results$significance <- ifelse(results$p_value < 0.001, "X", "")

message("[wilcoxon] Writing results to: ", output_path)
fwrite(results, output_path)

if (exists("main_gradients", inherits = FALSE) && exists("permuted_gradients", inherits = FALSE) &&
    nrow(permuted_gradients) > 0) {
    .report_gradient_summary(main_gradients, permuted_gradients, results)
}

message("[wilcoxon] Done. ", nrow(results), " tests written.")
