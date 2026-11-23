# Set the working directory to assessment2 before running.
a2_dir <- "."
raw <- read.csv(file.path(a2_dir, "SHDS_Assessment2_data-1.csv"),
                na.strings = c("", "NA"))

continuous_labels <- c(
  HDL = "HDL (mmol/L)", Sleep = "Sleep (hours/day)",
  BMI = "BMI (kg/m^2)", Age = "Age at recruitment (years)"
)
category_names <- c(
  Irritability = "Irritability", Smoking = "Smoking status",
  Physical_activity = "Daily physical activity",
  Dementia = "True dementia status", Dementia_diag = "Doctor's diagnosis"
)
category_labels <- list(
  Irritability = c("No", "Yes"),
  Smoking = c("Non-smoker", "Ex-smoker", "Current smoker"),
  Physical_activity = c("15 min/day or less", "More than 15 min/day"),
  Dementia = c("No dementia", "Dementia"),
  Dementia_diag = c("Diagnosis negative", "Diagnosis positive")
)
required <- c("ID", names(continuous_labels), names(category_labels))
stopifnot(nrow(raw) > 0, all(required %in% names(raw)))
if (anyDuplicated(na.omit(raw$ID))) stop("Check duplicate patient IDs.")

cat("Missing values by column\n")
print(colSums(is.na(raw[required])))

for (v in names(continuous_labels)) {
  x <- raw[[v]]
  stopifnot(is.numeric(x), all(is.finite(x[!is.na(x)])))
}
range_flags <- with(raw, c(
  Sleep = sum(Sleep < 0 | Sleep > 24, na.rm = TRUE),
  BMI = sum(BMI <= 0, na.rm = TRUE),
  HDL = sum(HDL <= 0, na.rm = TRUE),
  Age = sum(Age <= 0, na.rm = TRUE)
))
if (any(range_flags > 0)) {
  print(range_flags)
  stop("Review values outside the basic plausibility bounds.")
}

for (v in names(category_labels)) {
  labels <- category_labels[[v]]
  codes <- seq_along(labels) - 1L
  if (any(!is.na(raw[[v]]) & !raw[[v]] %in% codes)) {
    stop(paste("Unexpected category code in", v))
  }
  raw[[v]] <- factor(raw[[v]], levels = codes, labels = labels)
}

# Fix one complete-case sample for all planned analyses.
patients <- raw[complete.cases(raw[required]), required]
stopifnot(nrow(patients) > 0)
cat("Original:", nrow(raw), "| Excluded:", nrow(raw) - nrow(patients),
    "| Analysis sample:", nrow(patients), "\n")

# Flag extremes for review without dropping them. The brief gives no age cutoff.
range_checks <- do.call(rbind, lapply(names(continuous_labels), function(v) {
  x <- patients[[v]]
  data.frame(variable = v, minimum = min(x), maximum = max(x))
}))
print(range_checks, row.names = FALSE)
# HDL review thresholds used by Q1 and Q2 (mmol/L).
hdl_review_flag <- with(patients, HDL < 0.5 | HDL > 3.5)
