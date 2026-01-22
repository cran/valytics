## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4,
  fig.align = "center"
)

## ----load---------------------------------------------------------------------
library(valytics)
library(ggplot2)

## ----ba-example---------------------------------------------------------------
data("creatinine_serum")
ba <- ba_analysis(
  x = creatinine_serum$enzymatic,
  y = creatinine_serum$jaffe
)

## ----bias-output--------------------------------------------------------------
cat("Bias:", round(ba$results$bias, 3), "mg/dL\n")
cat("95% CI:", round(ba$results$bias_ci["lower"], 3), "to",
    round(ba$results$bias_ci["upper"], 3), "\n")

## ----loa-output---------------------------------------------------------------
cat("Lower LoA:", round(ba$results$loa_lower, 3), "\n")
cat("Upper LoA:", round(ba$results$loa_upper, 3), "\n")
cat("Width:", round(ba$results$loa_upper - ba$results$loa_lower, 3), "\n")

## ----ba-plot, fig.cap = "Bland-Altman plot showing differences vs. averages."----
plot(ba)

## ----normality----------------------------------------------------------------
summ <- summary(ba)
if (!is.null(summ$normality_test)) {
  cat("Shapiro-Wilk p-value:", round(summ$normality_test$p.value, 4), "\n")
}

## ----histogram, fig.cap = "Distribution of differences."----------------------
ggplot(data.frame(diff = ba$results$differences), aes(x = diff)) +
  geom_histogram(aes(y = after_stat(density)), bins = 15,
                 fill = "steelblue", alpha = 0.7) +
  geom_density(linewidth = 1) +
  labs(x = "Difference (Jaffe - Enzymatic)", y = "Density") +
  theme_minimal()

## ----pb-example---------------------------------------------------------------
pb <- pb_regression(
  x = creatinine_serum$enzymatic,
  y = creatinine_serum$jaffe
)

## ----pb-output----------------------------------------------------------------
cat("Slope:", round(pb$results$slope, 4), "\n")
cat("  95% CI:", round(pb$results$slope_ci["lower"], 4), "to",
    round(pb$results$slope_ci["upper"], 4), "\n")
cat("Intercept:", round(pb$results$intercept, 4), "\n")
cat("  95% CI:", round(pb$results$intercept_ci["lower"], 4), "to",
    round(pb$results$intercept_ci["upper"], 4), "\n")

## ----translation--------------------------------------------------------------
# At various concentrations, what's the expected difference?
concentrations <- c(0.8, 1.3, 3.0, 6.0)

for (conc in concentrations) {
  expected_y <- pb$results$intercept + pb$results$slope * conc
  difference <- expected_y - conc
  cat(sprintf("At X = %.1f: expected Y = %.3f, difference = %.3f\n",
              conc, expected_y, difference))
}

## ----cusum--------------------------------------------------------------------
cat("CUSUM statistic:", round(pb$cusum$statistic, 4), "\n")
cat("p-value:", round(pb$cusum$p_value, 4), "\n")

## ----cusum-plot, fig.cap = "CUSUM plot for linearity assessment."-------------
plot(pb, type = "cusum")

## ----correlation--------------------------------------------------------------
r <- cor(creatinine_serum$enzymatic, creatinine_serum$jaffe)
cat("Correlation coefficient:", round(r, 4), "\n")

## ----context------------------------------------------------------------------
# Example: Is a bias of X clinically meaningful?
# This depends entirely on YOUR application
bias_value <- ba$results$bias

cat("Observed bias:", round(bias_value, 3), "mg/dL\n")
cat("\nWhether this is 'acceptable' depends on:\n")
cat("- Your specific clinical decision thresholds\n")
cat("- Regulatory requirements for your application\n")
cat("- Intended use of the measurement\n")
cat("- Established performance goals (CLIA, biological variation, etc.)\n")

## ----report-------------------------------------------------------------------
# Bland-Altman summary
cat("=== Bland-Altman Analysis ===\n")
cat(sprintf("n = %d\n", ba$input$n))
cat(sprintf("Bias: %.3f (95%% CI: %.3f to %.3f)\n",
            ba$results$bias,
            ba$results$bias_ci["lower"],
            ba$results$bias_ci["upper"]))
cat(sprintf("SD of differences: %.3f\n", ba$results$sd_diff))
cat(sprintf("LoA: %.3f to %.3f\n\n",
            ba$results$loa_lower,
            ba$results$loa_upper))

# Passing-Bablok summary
cat("=== Passing-Bablok Regression ===\n")
cat(sprintf("Slope: %.4f (95%% CI: %.4f to %.4f)\n",
            pb$results$slope,
            pb$results$slope_ci["lower"],
            pb$results$slope_ci["upper"]))
cat(sprintf("Intercept: %.4f (95%% CI: %.4f to %.4f)\n",
            pb$results$intercept,
            pb$results$intercept_ci["lower"],
            pb$results$intercept_ci["upper"]))
cat(sprintf("CUSUM p-value: %.4f\n", pb$cusum$p_value))

## ----comparison-table, echo = FALSE-------------------------------------------
comparison_df <- data.frame(
  Aspect = c(
    "Primary question",
    "Statistical approach",
    "Error assumption",
    "Outlier handling",
    "Output focus",
    "Sample size",
    "Best when"
  ),
  `Bland-Altman` = c(
    "How well do methods agree?",
    "Descriptive statistics",
    "Differences ~ Normal",
    "Sensitive",
    "Bias, limits of agreement",
    "n >= 30 recommended",
    "Defining acceptable agreement"
  ),
  `Passing-Bablok` = c(
    "Is there systematic bias?",
    "Non-parametric regression",
    "Distribution-free",
    "Robust",
    "Slope, intercept CIs",
    "n >= 30 for stable CIs",
    "Outliers present, unknown error"
  ),
  Deming = c(
    "Is there systematic bias?",
    "Parametric regression",
    "Errors ~ Normal",
    "Sensitive",
    "Slope, intercept, SEs",
    "n >= 10 feasible",
    "Known error ratio, small n"
  ),
  check.names = FALSE
)

knitr::kable(comparison_df, caption = "Comparison of method comparison approaches")

## ----multiple-methods, eval = FALSE-------------------------------------------
#  # Complete method comparison workflow
#  ba <- ba_analysis(reference ~ test, data = mydata)
#  pb <- pb_regression(reference ~ test, data = mydata)
#  dm <- deming_regression(reference ~ test, data = mydata)
#  
#  # Bland-Altman for agreement assessment
#  summary(ba)
#  plot(ba)
#  
#  # Compare regression methods
#  cat("Passing-Bablok slope:", pb$results$slope, "\n")
#  cat("Deming slope:", dm$results$slope, "\n")

