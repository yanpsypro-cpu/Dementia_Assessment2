# H0: all three population means are equal. H1: at least two differ.
source("00_prepare.R", local = TRUE)
output_dir <- file.path(a2_dir, "q2_output")
dir.create(output_dir, showWarnings = FALSE)

group_stats <- do.call(rbind, lapply(split(patients$HDL, patients$Smoking), function(x) {
  q <- quantile(x, c(.25, .75), names = FALSE)
  data.frame(n = length(x), mean = mean(x), SD = sd(x), median = median(x),
             Q1 = q[1], Q3 = q[2], IQR = IQR(x))
}))
group_stats <- data.frame(Smoking = rownames(group_stats), group_stats, row.names = NULL)
ci_margin <- with(group_stats, qt(.975, n - 1) * SD / sqrt(n))
group_stats$mean_ci_lower <- group_stats$mean - ci_margin
group_stats$mean_ci_upper <- group_stats$mean + ci_margin
print(group_stats)
write.csv(group_stats, file.path(output_dir, "group_summary.csv"), row.names = FALSE)

png(file.path(output_dir, "group_mean_intervals.png"),
    width = 1200, height = 780, res = 145)
par(mar = c(6, 5, 1.5, 2))
group_position <- seq_len(nrow(group_stats))
plot(group_position, group_stats$mean, type = "n", xaxt = "n",
     xlab = "", ylab = "Mean HDL (mmol/L)",
     xlim = c(0.5, 3.5),
     ylim = range(group_stats$mean_ci_lower, group_stats$mean_ci_upper))
axis(1, at = group_position,
     labels = sprintf("%s\n(n = %d)", group_stats$Smoking, group_stats$n))
arrows(group_position, group_stats$mean_ci_lower, group_position,
       group_stats$mean_ci_upper, angle = 90, code = 3, length = 0.07,
       col = "steelblue", lwd = 2)
points(group_position, group_stats$mean, pch = 19, col = "steelblue4")
invisible(dev.off())

# Independence is assessed from the sampling design. Check the residuals
# graphically and use median-centred Levene's test to assess equal variances.
fit <- aov(HDL ~ Smoking, data = patients)
png(file.path(output_dir, "anova_diagnostics.png"), width = 1950, height = 720, res = 145)
par(mfrow = c(1, 3))
qqnorm(residuals(fit), pch = 16, cex = 0.55, col = "steelblue")
qqline(residuals(fit), col = "gray40")
plot(fitted(fit), residuals(fit), xlab = "Fitted HDL (mmol/L)",
     ylab = "Residual HDL (mmol/L)", main = "Residuals vs fitted",
     pch = 16, cex = 0.55, col = "steelblue")
abline(h = 0, lty = 2, col = "gray40")
plot(fit, which = 4, col = "chocolate3")
invisible(dev.off())

# Requires the car package (install once with install.packages("car")).
levene <- car::leveneTest(HDL ~ Smoking, data = patients)
cat("\nMedian-centred Levene test\n")
print(levene)

cat("\nPrimary one-way ANOVA\n")
primary <- anova(fit)
print(primary)
cat("\nWelch ANOVA on the same sample\n")
welch <- oneway.test(HDL ~ Smoking, data = patients, var.equal = FALSE)
print(welch)
cat("\nKruskal-Wallis rank test (supplementary)\n")
rank_test <- kruskal.test(HDL ~ Smoking, data = patients)
print(rank_test)
write.csv(data.frame(test = "Kruskal-Wallis", n = nrow(patients),
                     chi_squared = unname(rank_test$statistic),
                     df = unname(rank_test$parameter), p_value = rank_test$p.value),
          file.path(output_dir, "rank_test.csv"), row.names = FALSE)
sum_sq <- primary[["Sum Sq"]]
eta_squared <- sum_sq[1] / sum(sum_sq)
cat("Eta squared:", eta_squared, "\n")
write.csv(data.frame(eta_squared = eta_squared, percent = 100 * eta_squared),
          file.path(output_dir, "effect_size.csv"), row.names = FALSE)

# Exclude only the two observations flagged in Q1 for this sensitivity check.
# They remain in the primary analysis.
excluded <- hdl_review_flag
stopifnot(sum(excluded) == 2)
cat("\nObservations omitted in sensitivity analysis\n")
print(patients[excluded, c("ID", "HDL", "Smoking")], row.names = FALSE)
sensitivity_data <- patients[!excluded, ]
sensitivity_fit <- aov(HDL ~ Smoking, data = sensitivity_data)
sensitivity_anova <- anova(sensitivity_fit)
sensitivity_welch <- oneway.test(HDL ~ Smoking, data = sensitivity_data,
                               var.equal = FALSE)
print(sensitivity_anova)
print(sensitivity_welch)

anova_row <- function(result, label, n) {
  data.frame(test = label, n = n, F = result[["F value"]][1],
             df1 = result$Df[1], df2 = result$Df[2],
             p_value = result[["Pr(>F)"]][1])
}
welch_row <- function(result, label, n) {
  data.frame(test = label, n = n, F = unname(result$statistic),
             df1 = unname(result$parameter[1]), df2 = unname(result$parameter[2]),
             p_value = result$p.value)
}
test_results <- rbind(
  anova_row(levene, "Median-centred Levene (Brown-Forsythe)", nrow(patients)),
  anova_row(primary, "Classical ANOVA (primary)", nrow(patients)),
  welch_row(welch, "Welch ANOVA (full sample)", nrow(patients)),
  anova_row(sensitivity_anova, "Classical ANOVA (sensitivity)", nrow(sensitivity_data)),
  welch_row(sensitivity_welch, "Welch ANOVA (sensitivity)", nrow(sensitivity_data))
)
write.csv(test_results, file.path(output_dir, "test_results.csv"), row.names = FALSE)

if (primary[["Pr(>F)"]][1] < 0.05) {
  primary_pairs <- TukeyHSD(fit)$Smoking
  cat("\nPrimary Tukey comparisons with 95% family-wise intervals\n")
  print(primary_pairs)
  write.csv(primary_pairs, file.path(output_dir, "tukey_primary.csv"))
  
  # Compare the same contrasts and simultaneous intervals across samples.
  sensitivity_pairs <- TukeyHSD(sensitivity_fit)$Smoking
  comparison <- data.frame(
    contrast = rownames(primary_pairs),
    primary_diff = primary_pairs[, "diff"],
    primary_lower = primary_pairs[, "lwr"],
    primary_upper = primary_pairs[, "upr"],
    sensitivity_diff = sensitivity_pairs[, "diff"],
    sensitivity_lower = sensitivity_pairs[, "lwr"],
    sensitivity_upper = sensitivity_pairs[, "upr"],
    row.names = NULL
  )
  print(comparison)
  write.csv(comparison, file.path(output_dir, "contrast_sensitivity.csv"), row.names = FALSE)
}
