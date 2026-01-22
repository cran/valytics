## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4,
  fig.align = "center"
)

## ----load-package-------------------------------------------------------------
library(valytics)
library(ggplot2)

## ----attenuation-demo, fig.cap = "Demonstration of attenuation bias: OLS slope is attenuated when X has measurement error."----
set.seed(42)

# True relationship: Y = 1.0 * X (perfect agreement)
true_values <- seq(50, 150, length.out = 100)

# Both methods have measurement error
x_observed <- true_values + rnorm(100, sd = 10)
y_observed <- true_values + rnorm(100, sd = 10)

# Compare OLS vs Deming
ols_fit <- lm(y_observed ~ x_observed)
dm_fit <- deming_regression(x_observed, y_observed)

# Visualize
df <- data.frame(x = x_observed, y = y_observed)

ggplot(df, aes(x = x, y = y)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", 
              color = "gray50", linewidth = 1) +
  geom_abline(intercept = coef(ols_fit)[1], slope = coef(ols_fit)[2],
              color = "red", linewidth = 0.8) +
  geom_abline(intercept = dm_fit$results$intercept, slope = dm_fit$results$slope,
              color = "blue", linewidth = 0.8) +
  annotate("text", x = 60, y = 145, label = "True (slope = 1)", color = "gray30") +
  annotate("text", x = 60, y = 138, 
           label = sprintf("OLS (slope = %.3f)", coef(ols_fit)[2]), color = "red") +
  annotate("text", x = 60, y = 131,
           label = sprintf("Deming (slope = %.3f)", dm_fit$results$slope), color = "blue") +
  labs(title = "Attenuation Bias in OLS Regression",
       x = "Method X", y = "Method Y") +
  theme_minimal()

## ----basic-usage--------------------------------------------------------------
# Load example data
data("glucose_methods")

# Deming regression with default settings (lambda = 1)
dm <- deming_regression(reference ~ poc_meter, data = glucose_methods)
dm

## ----summary------------------------------------------------------------------
summary(dm)

## ----scatter-plot, fig.cap = "Deming regression scatter plot with confidence band."----
plot(dm)

## ----residual-plot, fig.cap = "Residual plot for assessing model fit."--------
plot(dm, type = "residuals")

## ----error-ratio--------------------------------------------------------------
# If POC meter has twice the error variance of the reference
dm_lambda2 <- deming_regression(
  reference ~ poc_meter, 
  data = glucose_methods,
  error_ratio = 2
)

dm_lambda2

## ----cv-estimation, eval = FALSE----------------------------------------------
#  # Example: Reference CV = 2.5%, POC CV = 4.5%
#  cv_reference <- 0.025
#  cv_poc <- 0.045
#  
#  lambda_estimated <- (cv_poc / cv_reference)^2
#  # lambda_estimated = 3.24
#  
#  dm_cv <- deming_regression(
#    reference ~ poc_meter,
#    data = glucose_methods,
#    error_ratio = lambda_estimated
#  )

## ----jackknife----------------------------------------------------------------
dm_jack <- deming_regression(
  reference ~ poc_meter,
  data = glucose_methods,
  ci_method = "jackknife"
)

# Standard errors are available
cat("Slope SE:", round(dm_jack$results$slope_se, 4), "\n")
cat("Intercept SE:", round(dm_jack$results$intercept_se, 4), "\n")

## ----bootstrap, eval = FALSE--------------------------------------------------
#  dm_boot <- deming_regression(
#    reference ~ poc_meter,
#    data = glucose_methods,
#    ci_method = "bootstrap",
#    boot_n = 1999
#  )

## ----comparison, fig.cap = "Comparison of Deming and Passing-Bablok regression."----
# Fit both models
dm <- deming_regression(reference ~ poc_meter, data = glucose_methods)
pb <- pb_regression(reference ~ poc_meter, data = glucose_methods)

# Compare coefficients
comparison <- data.frame(
  Method = c("Deming", "Passing-Bablok"),
  Slope = c(dm$results$slope, pb$results$slope),
  Slope_Lower = c(dm$results$slope_ci["lower"], pb$results$slope_ci["lower"]),
  Slope_Upper = c(dm$results$slope_ci["upper"], pb$results$slope_ci["upper"]),
  Intercept = c(dm$results$intercept, pb$results$intercept),
  Int_Lower = c(dm$results$intercept_ci["lower"], pb$results$intercept_ci["lower"]),
  Int_Upper = c(dm$results$intercept_ci["upper"], pb$results$intercept_ci["upper"])
)

print(comparison, digits = 3)

## ----comparison-plot, fig.cap = "Visual comparison of Deming and Passing-Bablok regression lines."----
# Visual comparison
ggplot(glucose_methods, aes(x = reference, y = poc_meter)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_abline(intercept = dm$results$intercept, slope = dm$results$slope,
              color = "#2166AC", linewidth = 1) +
  geom_abline(intercept = pb$results$intercept, slope = pb$results$slope,
              color = "#B2182B", linewidth = 1) +
  annotate("text", x = 80, y = 340, label = "Identity", color = "gray50") +
  annotate("text", x = 80, y = 320, label = "Deming", color = "#2166AC") +
  annotate("text", x = 80, y = 300, label = "Passing-Bablok", color = "#B2182B") +
  labs(title = "Regression Method Comparison",
       x = "Reference (mg/dL)",
       y = "POC Meter (mg/dL)") +
  theme_minimal()

## ----workflow-----------------------------------------------------------------
# Load data
data("creatinine_serum")

# 1. Deming regression
dm <- deming_regression(enzymatic ~ jaffe, data = creatinine_serum)

# 2. View summary
summary(dm)

## ----workflow-plot, fig.cap = "Complete Deming regression analysis for creatinine methods."----
# 3. Create visualization
plot(dm, title = "Creatinine: Jaffe vs Enzymatic Method")

## ----workflow-residuals, fig.cap = "Residual diagnostics for creatinine comparison."----
# 4. Check residuals
plot(dm, type = "residuals")

## ----extraction---------------------------------------------------------------
# Coefficients
slope <- dm$results$slope
intercept <- dm$results$intercept

# Confidence intervals
slope_ci <- dm$results$slope_ci
intercept_ci <- dm$results$intercept_ci

# Standard errors
slope_se <- dm$results$slope_se
intercept_se <- dm$results$intercept_se

# Formatted output for reporting
cat(sprintf("Slope: %.4f (95%% CI: %.4f to %.4f)\n", 
            slope, slope_ci["lower"], slope_ci["upper"]))
cat(sprintf("Intercept: %.4f (95%% CI: %.4f to %.4f)\n",
            intercept, intercept_ci["lower"], intercept_ci["upper"]))

