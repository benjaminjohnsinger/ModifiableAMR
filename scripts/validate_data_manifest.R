#!/usr/bin/env Rscript

source("config.R")

manifest_path <- file.path(AMR_CONFIG$paths$manifests, "data_manifest.csv")
if (!file.exists(manifest_path)) {
  stop("Data manifest not found: ", manifest_path)
}

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- manifest[manifest$status == "required", , drop = FALSE]

missing <- character(0)
for (i in seq_len(nrow(required))) {
  p <- required$path[i]
  if (!file.exists(p)) {
    missing <- c(missing, p)
  }
}

if (length(missing) > 0) {
  message("[validate-inputs] Missing required inputs:")
  for (p in missing) message("  - ", p)
  quit(status = 1)
}

checksum_issues <- character(0)
for (i in seq_len(nrow(required))) {
  expected <- required$checksum_sha256[i]
  p <- required$path[i]
  if (identical(expected, "TBD") || expected == "") {
    next
  }

  out <- tryCatch(
    system2("shasum", c("-a", "256", shQuote(p)), stdout = TRUE, stderr = TRUE),
    error = function(e) e$message
  )
  if (!length(out) || grepl("not found|Error|cannot", out[1], ignore.case = TRUE)) {
    checksum_issues <- c(checksum_issues, sprintf("Could not compute checksum for %s", p))
    next
  }

  actual <- strsplit(out[1], "\\s+")[[1]][1]
  if (!identical(tolower(actual), tolower(expected))) {
    checksum_issues <- c(
      checksum_issues,
      sprintf("Checksum mismatch for %s (expected %s, actual %s)", p, expected, actual)
    )
  }
}

if (length(checksum_issues) > 0) {
  message("[validate-inputs] Checksum issues:")
  for (m in checksum_issues) message("  - ", m)
  quit(status = 1)
}

message("[validate-inputs] Required inputs present and checksum checks passed/skipped.")
