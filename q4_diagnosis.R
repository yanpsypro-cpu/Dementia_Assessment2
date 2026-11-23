# Q4: Assess the doctor's diagnosis against the recorded true status.
source("00_prepare.R", local = TRUE)
output_dir <- file.path(a2_dir, "q4_output")
dir.create(output_dir, showWarnings = FALSE)

tab <- with(patients, table(Diagnosis = Dementia_diag, Truth = Dementia))
tab <- tab[c("Diagnosis positive", "Diagnosis negative"), c("Dementia", "No dementia")]
stopifnot(sum(tab) == nrow(patients))
print(addmargins(tab))
write.csv(addmargins(tab), file.path(output_dir, "diagnosis_table.csv"))

tp <- tab["Diagnosis positive", "Dementia"]
fp <- tab["Diagnosis positive", "No dementia"]
fn <- tab["Diagnosis negative", "Dementia"]
tn <- tab["Diagnosis negative", "No dementia"]

# Exact binomial confidence intervals.
# Sensitivity has 52 observations, including four false negatives.
metric <- function(name, x, n) {
  interval <- binom.test(x, n)$conf.int
  data.frame(metric = name, numerator = x, denominator = n,
             estimate = x / n, lower_95 = interval[1], upper_95 = interval[2])
}
accuracy <- rbind(
  metric("Sensitivity", tp, tp + fn),
  metric("Specificity", tn, tn + fp),
  metric("PPV", tp, tp + fp),
  metric("NPV", tn, tn + fn),
  metric("Accuracy", tp + tn, sum(tab))
)
print(accuracy, row.names = FALSE)
write.csv(accuracy, file.path(output_dir, "diagnostic_accuracy.csv"), row.names = FALSE)

# Use unrounded estimates to verify PPV and NPV.
sensitivity <- tp / (tp + fn)
specificity <- tn / (tn + fp)
prevalence <- (tp + fn) / sum(tab)
ppv_from_prevalence <- function(p) {
  sensitivity * p / (sensitivity * p + (1 - specificity) * (1 - p))
}
npv_from_prevalence <- function(p) {
  specificity * (1 - p) /
    (specificity * (1 - p) + (1 - sensitivity) * p)
}
stopifnot(
  isTRUE(all.equal(ppv_from_prevalence(prevalence), tp / (tp + fp))),
  isTRUE(all.equal(npv_from_prevalence(prevalence), tn / (tn + fn)))
)
# Illustrative scenarios assume sensitivity and specificity remain constant.
scenarios <- data.frame(scenario = c("Observed", "Illustrative", "Illustrative"),
                        prevalence = c(prevalence, 0.10, 0.20))
scenarios$PPV <- ppv_from_prevalence(scenarios$prevalence)
scenarios$NPV <- npv_from_prevalence(scenarios$prevalence)
print(scenarios, row.names = FALSE)
write.csv(scenarios, file.path(output_dir, "ppv_by_prevalence.csv"), row.names = FALSE)
cat("Recorded true prevalence:", prevalence,
    "| Diagnosis-positive proportion:", (tp + fp) / sum(tab), "\n")
