# Q3: Explore associations with true dementia status.
source("00_prepare.R", local = TRUE)

output_dir <- file.path(a2_dir, "q3_output")
dir.create(output_dir, showWarnings = FALSE)
stopifnot(identical(levels(patients$Dementia), c("No dementia", "Dementia")))
cat("Patients:", nrow(patients), " Dementia cases:", sum(patients$Dementia == "Dementia"), "\n")

# Rescale continuous predictors so each coefficient has the reporting unit.
patients$Age10 <- patients$Age / 10
patients$HDL01 <- patients$HDL / 0.1
full <- glm(Dementia ~ Age10 + BMI + HDL01 + Sleep + Irritability +
              Smoking + Physical_activity,
            data = patients, family = binomial(), na.action = na.fail)
null <- update(full, . ~ 1)
aic_model <- step(full, scope = list(lower = ~ Age10), direction = "backward", trace = 0)
bic_model <- step(full, scope = list(lower = ~ Age10), direction = "backward",
                  k = log(nobs(full)), trace = 0)
models <- list(Null = null, Full = full, AIC = aic_model, BIC = bic_model)
model_summary <- do.call(rbind, lapply(names(models), function(label) {
  m <- models[[label]]
  data.frame(model = label, formula = paste(deparse(formula(m)), collapse = " "),
             n = nobs(m), parameters_including_intercept = m$rank, deviance = deviance(m),
             AIC = AIC(m), BIC = BIC(m))
}))
print(model_summary, row.names = FALSE)
write.csv(model_summary, file.path(output_dir, "models.csv"), row.names = FALSE)
write.csv(aic_model$anova, file.path(output_dir, "aic_path.csv"), row.names = FALSE)
write.csv(bic_model$anova, file.path(output_dir, "bic_path.csv"), row.names = FALSE)

# Wald intervals on the log-odds scale, then exponentiate.
or_table <- function(m) {
  estimates <- coef(m)
  intervals <- confint.default(m)
  data.frame(term = names(estimates), OR = exp(unname(estimates)),
             lower = exp(intervals[, 1]), upper = exp(intervals[, 2]),
             row.names = NULL)[names(estimates) != "(Intercept)", ]
}
terms_full <- attr(terms(full), "term.labels")
crude <- do.call(rbind, lapply(terms_full, function(v) {
  or_table(glm(reformulate(v, response = "Dementia"), data = patients,
               family = binomial(), na.action = na.fail))
}))
or_results <- data.frame(term = crude$term)
or_models <- list(crude = crude, full = or_table(full), aic = or_table(aic_model))
for (label in names(or_models)) {
  estimates <- or_models[[label]]
  or_results[paste0(label, c("_OR", "_lower", "_upper"))] <-
    estimates[match(or_results$term, estimates$term), c("OR", "lower", "upper")]
}
print(or_results, row.names = FALSE)
write.csv(or_results, file.path(output_dir, "odds_ratios.csv"), row.names = FALSE)

# Compare only nested models fitted to identical rows.
lrt_row <- function(smaller, larger, label) {
  if (!all(attr(terms(smaller), "term.labels") %in% attr(terms(larger), "term.labels")) ||
      smaller$rank >= larger$rank) return(NULL)
  stopifnot(identical(rownames(model.frame(smaller)), rownames(model.frame(larger))))
  result <- anova(smaller, larger, test = "LRT")
  data.frame(comparison = label, chi_squared = result$Deviance[2],
             df = result$Df[2], p_value = result[["Pr(>Chi)"]][2])
}
comparisons <- rbind(
  lrt_row(null, full, "Null vs full"),
  lrt_row(null, aic_model, "Null vs AIC"),
  lrt_row(aic_model, full, "AIC vs full"),
  lrt_row(bic_model, full, "BIC vs full"),
  lrt_row(bic_model, aic_model, "BIC vs AIC")
)
print(comparisons, row.names = FALSE)
write.csv(comparisons, file.path(output_dir, "model_lrt.csv"), row.names = FALSE)

# Drop each term as a whole. Smoking has two degrees of freedom.
term_tests <- do.call(rbind, lapply(c("Full", "AIC"), function(label) {
  tests <- drop1(models[[label]], test = "LRT")[-1, , drop = FALSE]
  data.frame(model = label, term = rownames(tests), df = tests$Df,
             chi_squared = tests$LRT, p_value = tests[["Pr(>Chi)"]])
}))
print(term_tests, row.names = FALSE)
write.csv(term_tests, file.path(output_dir, "term_lrt.csv"), row.names = FALSE)

# Exploratory sleep comparison: separate the curvature test from the joint test.
sleep_linear <- update(aic_model, . ~ . + Sleep)
sleep_quadratic <- update(sleep_linear, . ~ . + I(Sleep^2))
sleep_tests <- rbind(
  lrt_row(sleep_linear, sleep_quadratic, "Linear sleep vs quadratic sleep"),
  lrt_row(aic_model, sleep_quadratic, "AIC model vs both sleep terms")
)
print(sleep_tests, row.names = FALSE)
write.csv(sleep_tests, file.path(output_dir, "sleep_lrt.csv"), row.names = FALSE)

age_by_smoking <- aggregate(Age ~ Smoking, data = patients, FUN = mean)
print(age_by_smoking, row.names = FALSE)
write.csv(age_by_smoking, file.path(output_dir, "age_by_smoking.csv"), row.names = FALSE)

# Inspect the three highest Cook's distances and refit the fixed AIC formula.
cook <- cooks.distance(aic_model)
top_rows <- order(cook, decreasing = TRUE)[1:3]
influential <- data.frame(row = top_rows,
                          patients[top_rows, c("ID", "Dementia", "Smoking")],
                          cooks_distance = unname(cook[top_rows]), row.names = NULL)
print(influential, row.names = FALSE)
write.csv(influential, file.path(output_dir, "influence_cases.csv"), row.names = FALSE)

png(file.path(output_dir, "cooks_distance.png"), width = 1200, height = 750, res = 145)
plot(aic_model, which = 4, col = "darkcyan")
invisible(dev.off())

# Inspect the largest points individually and together; do not reselect terms.
smoking_results <- function(model, scenario) {
  estimates <- or_table(model)
  data.frame(scenario = scenario, n = nobs(model),
             subset(estimates, startsWith(term, "Smoking")),
             row.names = NULL)
}
smoking_sensitivity <- rbind(
  smoking_results(aic_model, "All patients"),
  do.call(rbind, lapply(top_rows, function(i) {
    refit <- update(aic_model, data = patients[-i, ])
    smoking_results(refit, paste("Omit observation", i))
  })),
  smoking_results(update(aic_model, data = patients[-top_rows, ]),
                  "Omit all three")
)
print(smoking_sensitivity, row.names = FALSE)
write.csv(smoking_sensitivity, file.path(output_dir, "smoking_influence.csv"),
          row.names = FALSE)

# Compare crude estimates with estimates adjusted for all full-model terms.
forest_labels <- c(
  Age10 = "Age (per 10 years)", BMI = "BMI (per 1 kg/m²)",
  HDL01 = "HDL (per 0.1 mmol/L)", Sleep = "Sleep (per 1 hour)",
  IrritabilityYes = "Irritability: yes vs no",
  `SmokingEx-smoker` = "Ex-smoker vs non-smoker",
  `SmokingCurrent smoker` = "Current vs non-smoker",
  `Physical_activityMore than 15 min/day` = "Activity: >15 vs ≤15 min/day"
)
forest <- or_results[match(names(forest_labels), or_results$term), ]
stopifnot(identical(forest$term, names(forest_labels)))
y <- rev(seq_len(nrow(forest)))
limits <- range(1, forest[c("crude_lower", "crude_upper", "full_lower", "full_upper")])

png(file.path(output_dir, "or_forest.png"), width = 1500, height = 850, res = 150)
par(mar = c(4.5, 13, 1, 2))
plot(NA, xlim = limits, ylim = c(0, max(y) + 0.5), log = "x",
     xaxt = "n", yaxt = "n", xlab = "Odds ratio (log scale)", ylab = "", bty = "n")
axis(1, at = c(0.2, 0.5, 1, 2, 5, 10))
axis(2, at = y, labels = unname(forest_labels), las = 1, tick = FALSE, cex.axis = 0.85)
abline(v = 1, lty = 2, col = "gray50")
segments(forest$crude_lower, y + 0.15, forest$crude_upper, y + 0.15,
         col = "gray45", lwd = 1.5)
points(forest$crude_OR, y + 0.15, pch = 1, col = "gray35", lwd = 1.5)
segments(forest$full_lower, y - 0.15, forest$full_upper, y - 0.15,
         col = "darkcyan", lwd = 2)
points(forest$full_OR, y - 0.15, pch = 19, col = "darkcyan")
legend("bottomright", c("Crude", "Adjusted (full model)"), pch = c(1, 19),
       col = c("gray35", "darkcyan"), bty = "n", horiz = TRUE, cex = 0.86)
invisible(dev.off())

# Exploratory checks added during report revision; keep the original model.
patients$Current_smoking <- factor(patients$Smoking == "Current smoker",
                                  levels = c(FALSE, TRUE), labels = c("No", "Yes"))
collapsed <- update(full, . ~ . - Smoking + Current_smoking)
extra_models <- list(
  Collapsed_AIC = step(collapsed, scope = list(lower = ~ Age10), direction = "backward", trace = 0),
  Collapsed_BIC = step(collapsed, scope = list(lower = ~ Age10), direction = "backward",
                       k = log(nrow(patients)), trace = 0),
  Unforced_AIC = step(full, direction = "backward", trace = 0),
  Doctor_diagnosis = update(aic_model, Dementia_diag ~ .))
without_1790 <- subset(patients, ID != 1790)
bmi_base <- update(aic_model, data = without_1790)
extra_lrt <- rbind(
  lrt_row(aic_model, update(aic_model, . ~ . + I(BMI^2)), "BMI curvature, all patients"),
  lrt_row(bmi_base, update(bmi_base, . ~ . + I(BMI^2)), "BMI curvature, omit ID 1790"))
young <- with(patients, Dementia == "Dementia" & Smoking == "Current smoker" & Age < 60)
extra_models$Omit_current_smoker_cases_under60 <- update(aic_model, data = patients[!young, ])
extra_or <- do.call(rbind, lapply(names(extra_models), function(label) {
  data.frame(model = label, or_table(extra_models[[label]]))
}))
counts <- data.frame(
  measure = c("Cases under60", "Patients under60", "Cases 60plus", "Patients 60plus",
              "BMI 25plus", "Current smoker cases under60"),
  value = c(sum(patients$Dementia == "Dementia" & patients$Age < 60), sum(patients$Age < 60),
            sum(patients$Dementia == "Dementia" & patients$Age >= 60), sum(patients$Age >= 60),
            sum(patients$BMI >= 25), sum(young)))
write.csv(extra_lrt, file.path(output_dir, "additional_lrt.csv"), row.names = FALSE)
write.csv(extra_or, file.path(output_dir, "additional_or.csv"), row.names = FALSE)
write.csv(counts, file.path(output_dir, "additional_counts.csv"), row.names = FALSE)
print(extra_lrt)
print(extra_or)
print(counts)

# Q5: observed prevalence within each age group.
age_group <- factor(patients$Age >= 60, levels = c(FALSE, TRUE),
                    labels = c("Under 60", "60 and over"))
age_table <- table(age_group, patients$Dementia)
age_prevalence <- data.frame(
  age_group = rownames(age_table),
  cases = as.integer(age_table[, "Dementia"]),
  n = as.integer(rowSums(age_table)),
  prevalence = as.numeric(age_table[, "Dementia"] / rowSums(age_table))
)
print(age_prevalence, row.names = FALSE)
write.csv(age_prevalence, file.path(output_dir, "age_prevalence.csv"), row.names = FALSE)
