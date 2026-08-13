#!/usr/bin/env Rscript

compute_bootstrap_zero_p_value <- function(bootstrap_values) {
  x <- suppressWarnings(as.numeric(bootstrap_values))
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  prob_left <- mean(x <= 0, na.rm = TRUE)
  prob_right <- mean(x >= 0, na.rm = TRUE)
  p_value <- 2 * min(prob_left, prob_right)
  min(1, max(0, p_value))
}

load_gradient_table_with_labels <- function(path, label_col = NULL) {
  raw <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  if ("Gradient" %in% names(raw)) {
    return(raw)
  }

  rowname_df <- tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1),
    error = function(e) NULL
  )

  if (!is.null(rowname_df) && "Gradient" %in% names(rowname_df)) {
    return(as.data.frame(rowname_df, stringsAsFactors = FALSE))
  }

  if (!is.null(rowname_df) && nrow(rowname_df) > 0L) {
    label_name <- if (!is.null(label_col) && nzchar(label_col)) label_col else "Label"
    out <- as.data.frame(rowname_df, stringsAsFactors = FALSE)
    out[[label_name]] <- rownames(rowname_df)
    return(out)
  }

  raw
}

add_zero_effect_p_value <- function(gradients_path, bootstrap_path, output_path = gradients_path, label_col = NULL) {
  gradients <- load_gradient_table_with_labels(gradients_path, label_col = label_col)
  bootstraps <- read.csv(bootstrap_path, stringsAsFactors = FALSE)

  if (!"Gradient" %in% names(gradients)) {
    stop(sprintf("Gradient file is missing a 'Gradient' column: %s", gradients_path), call. = FALSE)
  }
  if (!"Gradient" %in% names(bootstraps)) {
    stop(sprintf("Bootstrap file is missing a 'Gradient' column: %s", bootstrap_path), call. = FALSE)
  }

  if (is.null(label_col)) {
    for (candidate in c("Antibiotic", "Pathogen", "Combo_Label", "Label")) {
      if (candidate %in% names(bootstraps) && candidate %in% names(gradients)) {
        label_col <- candidate
        break
      }
    }
  }

  if (is.null(label_col)) {
    label_col <- "Label"
  }

  bootstrap_label_col <- NULL
  if (!is.null(label_col) && label_col %in% names(bootstraps)) {
    bootstrap_label_col <- label_col
  } else {
    for (candidate in c("Antibiotic", "Pathogen", "Combo_Label", "Label")) {
      if (candidate %in% names(bootstraps)) {
        bootstrap_label_col <- candidate
        break
      }
    }
  }

  if (is.null(bootstrap_label_col)) {
    gradients$P_value <- compute_bootstrap_zero_p_value(bootstraps$Gradient)
  } else {
    p_by_label <- vapply(
      split(bootstraps$Gradient, as.character(bootstraps[[bootstrap_label_col]])),
      compute_bootstrap_zero_p_value,
      numeric(1)
    )

    if (label_col %in% names(gradients)) {
      label_values <- as.character(gradients[[label_col]])
      if (length(label_values) == length(p_by_label) && all(!is.na(suppressWarnings(as.numeric(label_values))))) {
        gradients[[label_col]] <- names(p_by_label)
        label_values <- as.character(gradients[[label_col]])
      }
    } else {
      gradients[[label_col]] <- names(p_by_label)
      label_values <- as.character(gradients[[label_col]])
    }

    gradients$P_value <- unname(p_by_label[match(label_values, names(p_by_label))])
  }

  if ("X" %in% names(gradients)) {
    gradients$X <- NULL
  }

  write.csv(gradients, output_path, row.names = FALSE)
  invisible(gradients)
}

parse_cli <- function(args) {
  if (length(args) == 0L) {
    return(list(
      tasks = list(
        list(
          gradients = "Outputs/database_gradients_ATC3_PCA_canonical_weighted_all.csv",
          bootstrap = "Outputs/database_gradients_bootstraps_ATC3_PCA_canonical_weighted_all.csv",
          output = "Outputs/database_gradients_ATC3_PCA_canonical_weighted_all.csv",
          label_col = "Antibiotic"
        ),
        list(
          gradients = "Outputs/database_gradients_pathogen_PCA_canonical_weighted_main.csv",
          bootstrap = "Outputs/database_gradients_bootstraps_pathogen_PCA_canonical_weighted_main.csv",
          output = "Outputs/database_gradients_pathogen_PCA_canonical_weighted_main.csv",
          label_col = "Pathogen"
        )
      )
    ))
  }

  named <- list(gradients = NULL, bootstrap = NULL, output = NULL, label_col = NULL)
  positional <- character()
  i <- 1L

  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--help") {
      message("Usage: Rscript scripts/add_bootstrap_zero_p_values.R [--gradients path --bootstrap path [--output path] [--label-col Name]]")
      quit(status = 0)
    }
    if (arg %in% c("--gradients", "--bootstrap", "--output", "--label-col")) {
      key <- sub("^--", "", arg)
      key <- gsub("-", "_", key)
      if (i == length(args)) {
        stop(sprintf("Missing value for %s", arg), call. = FALSE)
      }
      named[[key]] <- args[[i + 1L]]
      i <- i + 2L
    } else {
      positional <- c(positional, arg)
      i <- i + 1L
    }
  }

  if (!is.null(named$gradients) || !is.null(named$bootstrap) || !is.null(named$output) || !is.null(named$label_col)) {
    if (is.null(named$gradients) || is.null(named$bootstrap)) {
      stop("Both --gradients and --bootstrap must be supplied together.", call. = FALSE)
    }
    return(list(tasks = list(list(
      gradients = named$gradients,
      bootstrap = named$bootstrap,
      output = if (!is.null(named$output)) named$output else named$gradients,
      label_col = named$label_col
    ))))
  }

  if (length(positional) != 2L) {
    stop("Expected either no arguments or exactly two positional arguments: gradients_path bootstrap_path.", call. = FALSE)
  }

  list(tasks = list(list(
    gradients = positional[[1L]],
    bootstrap = positional[[2L]],
    output = positional[[1L]],
    label_col = NULL
  )))
}

main <- function() {
  cli <- parse_cli(commandArgs(trailingOnly = TRUE))

  for (task in cli$tasks) {
    message(sprintf("Computing zero-effect p-values for %s using %s", task$gradients, task$bootstrap))
    add_zero_effect_p_value(
      gradients_path = task$gradients,
      bootstrap_path = task$bootstrap,
      output_path = task$output,
      label_col = task$label_col
    )
    message(sprintf("Saved updated CSV to %s", task$output))
  }
}

main()
