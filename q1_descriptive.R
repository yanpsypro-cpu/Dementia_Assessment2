# Set the working directory to assessment2 before running.
source("00_prepare.R", local = TRUE)

output_dir <- file.path(a2_dir, "q1_output")
dir.create(output_dir, showWarnings = FALSE)

summary_fmt <- c(HDL = "median", Sleep = "median", BMI = "mean", Age = "mean")

summary_digits <- c(HDL = 2, Sleep = 0, BMI = 1, Age = 1)

numeric_stats <- function(x) {
  quartiles <- quantile(x, c(.25, .5, .75), names = FALSE)
  centered <- x - mean(x)
  variance <- mean(centered^2)
  c(n = length(x), mean = mean(x), sd = sd(x), median = quartiles[2],
    q1 = quartiles[1], q3 = quartiles[3],
    skewness = if (variance > 0) mean(centered^3) / variance^(3/2) else NA_real_)
}
diagnostics <- t(vapply(patients[names(summary_fmt)], numeric_stats, numeric(7)))
print(round(diagnostics, 3))

write.csv(data.frame(variable = rownames(diagnostics), diagnostics, row.names = NULL),
          file.path(output_dir, "numeric_diagnostics.csv"), row.names = FALSE)
write.csv(range_checks, file.path(output_dir, "range_checks.csv"), row.names = FALSE)

# Descriptive review thresholds. Keep these observations in the analysis sample.
hdl_review <- patients[hdl_review_flag, c("ID", "HDL", "Smoking")]
write.csv(hdl_review, file.path(output_dir, "hdl_extreme_values.csv"), row.names = FALSE)

format_numeric <- function(x, method, digits) {
  s <- numeric_stats(x)
  if (method == "median") {
    quantiles <- s[c("median", "q1", "q3")]
    if (digits == 0 && any(quantiles != round(quantiles))) digits <- 2
    sprintf("%.*f (%.*f to %.*f)", digits, quantiles[1],
            digits, quantiles[2], digits, quantiles[3])
  } else {
    sprintf("%.*f (%.*f)", digits, s["mean"], digits, s["sd"])
  }
}

subgroups <- list(
  Overall = patients,
  `No dementia` = patients[patients$Dementia == "No dementia", ],
  Dementia = patients[patients$Dementia == "Dementia", ]
)
stopifnot(all(vapply(subgroups, nrow, integer(1)) > 0))
headers <- sprintf("%s (n=%d)", names(subgroups),
                   vapply(subgroups, nrow, integer(1)))

numeric_rows <- do.call(rbind, lapply(names(continuous_labels), function(v) {
  method <- summary_fmt[[v]]
  values <- vapply(subgroups, function(group) format_numeric(group[[v]], method, summary_digits[[v]]), "")
  data.frame(variable = continuous_labels[[v]], category = "",
             statistic = if (method == "median") "Median (Q1 to Q3)" else "Mean (SD)",
             t(values), check.names = FALSE)
}))

category_rows <- do.call(rbind, lapply(names(category_labels), function(v) {
  values <- lapply(subgroups, function(group) {
    counts <- table(group[[v]])
    sprintf("%d (%.1f%%)", as.integer(counts), 100 * as.numeric(counts) / nrow(group))
  })
  # Report marginal status and diagnosis frequencies only.
  if (v %in% c("Dementia", "Dementia_diag")) {
    values[2:3] <- list(rep("", length(values[[1]])))
  }
  data.frame(variable = category_names[[v]], category = category_labels[[v]],
             statistic = "n (%)", values, check.names = FALSE)
}))

age_range <- vapply(subgroups, function(group) {
  paste(range(group$Age), collapse = " to ")
}, "")
age_range_row <- data.frame(variable = continuous_labels[["Age"]], category = "",
                            statistic = "Observed range", t(age_range),
                            check.names = FALSE)
table1 <- rbind(numeric_rows, age_range_row, category_rows)
names(table1)[4:6] <- headers
print(table1, row.names = FALSE)
write.csv(table1, file.path(output_dir, "sample_summary.csv"), row.names = FALSE)

png(file.path(output_dir, "numeric_distributions.png"),
    width = 1400, height = 1900, res = 150)
par(mfrow = c(4, 2), mar = c(4.2, 4.3, 2.3, 1))
for (v in names(continuous_labels)) {
  x <- patients[[v]]
  if (v == "Sleep" && all(x == floor(x))) {
    barplot(table(x), main = v, xlab = continuous_labels[[v]],
            ylab = "Frequency", col = "lightblue")
  } else {
    hist(x, breaks = "FD", main = v, xlab = continuous_labels[[v]],
         col = "lightblue")
  }
  qqnorm(x, main = paste(v, "Q-Q plot"), pch = 16, cex = 0.5)
  qqline(x)
}
invisible(dev.off())

png(file.path(output_dir, "hdl_by_smoking.png"), width = 1200, height = 800, res = 145)
smoking_labels <- sprintf("%s\n(n = %d)", levels(patients$Smoking),
                          as.integer(table(patients$Smoking)))
boxplot(HDL ~ Smoking, data = patients, names = smoking_labels,
        ylab = "HDL (mmol/L)", col = "lightsteelblue")
smoking_means <- tapply(patients$HDL, patients$Smoking, mean)
points(seq_along(smoking_means), smoking_means, pch = 19, col = "steelblue4")
legend("topright", legend = "Mean", pch = 19, col = "steelblue4", bty = "n")
invisible(dev.off())
