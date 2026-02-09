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

## ----load-data----------------------------------------------------------------
data("troponin_precision")
head(troponin_precision)

## ----design-summary-----------------------------------------------------------
# Summarize the experimental design
table(troponin_precision$level, troponin_precision$day)

## ----precision-study----------------------------------------------------------
prec <- precision_study(
 data = troponin_precision,
 value = "value",
 sample = "level",
 day = "day",
 run = "run"
)

print(prec)

## ----by-sample----------------------------------------------------------------
# Access results for each concentration level
names(prec$by_sample)

# Example: Level 1 (5 ng/L)
prec$by_sample$L1$precision

## ----precision-profile--------------------------------------------------------
profile <- precision_profile(
 prec,
 model = "hyperbolic",
 cv_targets = c(10, 20)
)

print(profile)

## ----parameters---------------------------------------------------------------
# Extract model parameters
a <- profile$model$parameters["a"]
b <- profile$model$parameters["b"]

cat(sprintf("Asymptotic CV (a): %.2f%%\n", a))
cat(sprintf("Concentration component (b): %.2f\n", b))
cat(sprintf("\nModel equation: %s\n", profile$model$equation))

## ----profile-plot, fig.cap = "Precision profile showing CV versus concentration with fitted hyperbolic model."----
plot(profile)

## ----profile-log, fig.cap = "Precision profile with logarithmic concentration scale."----
plot(profile, log_x = TRUE)

## ----functional-sensitivity---------------------------------------------------
profile$functional_sensitivity

## ----unachievable-------------------------------------------------------------
# Try a target CV of 2% (below asymptotic ~3%)
profile_strict <- precision_profile(
 prec,
 cv_targets = c(2, 5, 10)
)

profile_strict$functional_sensitivity

## ----bootstrap, eval = FALSE--------------------------------------------------
#  # Note: This takes longer to run
#  profile_boot <- precision_profile(
#   prec,
#   cv_targets = c(10, 20),
#   bootstrap = TRUE,
#   boot_n = 999
#  )
#  
#  profile_boot$functional_sensitivity

## ----linear-model-------------------------------------------------------------
profile_linear <- precision_profile(
 prec,
 model = "linear",
 cv_targets = c(10, 20)
)

print(profile_linear)

## ----compare-models-----------------------------------------------------------
cat("Hyperbolic model R²:", round(profile$fit_quality$r_squared, 4), "\n")
cat("Linear model R²:", round(profile_linear$fit_quality$r_squared, 4), "\n")

## ----summary, eval = FALSE----------------------------------------------------
#  summary(profile)

## ----dataframe-input----------------------------------------------------------
# Create a data frame with concentration and CV values
cv_data <- data.frame(
 concentration = c(5, 10, 25, 50, 100, 500),
 cv_pct = c(5.8, 4.2, 3.5, 3.2, 3.1, 3.0)
)

profile_df <- precision_profile(
 cv_data,
 concentration = "concentration",
 cv = "cv_pct"
)

print(profile_df)

## ----clinical-check-----------------------------------------------------------
# Assume 99th percentile URL = 20 ng/L
url_99th <- 20

# Check functional sensitivity at 10% CV
fs_10pct <- profile$functional_sensitivity$concentration[
 profile$functional_sensitivity$cv_target == 10
]

cat(sprintf("99th percentile URL: %.0f ng/L\n", url_99th))
cat(sprintf("Functional sensitivity (10%% CV): %.1f ng/L\n", fs_10pct))
cat(sprintf("Ratio (FS / URL): %.2f\n", fs_10pct / url_99th))

if (fs_10pct <= 0.5 * url_99th) {
 cat("\n✓ Meets criterion: FS ≤ 50% of 99th percentile URL\n")
} else {
 cat("\n✗ Does not meet criterion: FS > 50% of 99th percentile URL\n")
}

## ----workflow-summary, eval = FALSE-------------------------------------------
#  # Complete workflow
#  data("troponin_precision")
#  
#  # Step 1: Precision study
#  prec <- precision_study(
#   data = troponin_precision,
#   value = "value",
#   sample = "level",
#   day = "day",
#   run = "run"
#  )
#  
#  # Step 2: Precision profile
#  profile <- precision_profile(
#   prec,
#   model = "hyperbolic",
#   cv_targets = c(10, 20)
#  )
#  
#  # Step 3: Interpret and visualize
#  summary(profile)
#  plot(profile, log_x = TRUE)

