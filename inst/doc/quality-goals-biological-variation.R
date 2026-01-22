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

## ----ate-basic----------------------------------------------------------------
# Example: Glucose
# CV_I = 5.6%, CV_G = 7.5% (illustrative values)
ate_glucose <- ate_from_bv(cvi = 5.6, cvg = 7.5)
ate_glucose

## ----ate-summary--------------------------------------------------------------
summary(ate_glucose)

## ----ate-levels---------------------------------------------------------------
# Optimal (most stringent)
ate_optimal <- ate_from_bv(cvi = 5.6, cvg = 7.5, level = "optimal")
ate_optimal$specifications$tea

# Minimum (least stringent)
ate_minimum <- ate_from_bv(cvi = 5.6, cvg = 7.5, level = "minimum")
ate_minimum$specifications$tea

## ----ate-cvi-only-------------------------------------------------------------
ate_cv_only <- ate_from_bv(cvi = 5.6)
ate_cv_only

## ----sigma-basic--------------------------------------------------------------
# Assume observed: bias = 1.5%, CV = 2.5%
# Using TEa from biological variation
sm <- sigma_metric(
 bias = 1.5,
  cv = 2.5,
  tea = ate_glucose$specifications$tea
)
sm

## ----sigma-summary------------------------------------------------------------
summary(sm)

## ----assess-basic-------------------------------------------------------------
assess <- ate_assessment(
  bias = 1.5,
  cv = 2.5,
  tea = ate_glucose$specifications$tea
)
assess

## ----assess-full--------------------------------------------------------------
assess_full <- ate_assessment(
  bias = 1.5,
  cv = 2.5,
  tea = ate_glucose$specifications$tea,
  allowable_bias = ate_glucose$specifications$allowable_bias,
  allowable_cv = ate_glucose$specifications$allowable_cv
)
summary(assess_full)

## ----assess-fail--------------------------------------------------------------
# A method with poor performance
assess_poor <- ate_assessment(
  bias = 4.0,
  cv = 5.0,
  tea = ate_glucose$specifications$tea,
  allowable_bias = ate_glucose$specifications$allowable_bias,
  allowable_cv = ate_glucose$specifications$allowable_cv
)
summary(assess_poor)

## ----workflow-----------------------------------------------------------------
# Step 1: Define quality goals from biological variation
specs <- ate_from_bv(cvi = 5.6, cvg = 7.5, level = "desirable")
cat("Quality Specifications:\n")
cat(sprintf("  Allowable CV: %.2f%%\n", specs$specifications$allowable_cv))
cat(sprintf("  Allowable Bias: %.2f%%\n", specs$specifications$allowable_bias))
cat(sprintf("  TEa: %.2f%%\n\n", specs$specifications$tea))

# Step 2: Assume we measured method performance
# (In practice, from validation studies)
observed_bias <- 1.8
observed_cv <- 2.2

# Step 3: Calculate sigma metric
sm <- sigma_metric(observed_bias, observed_cv, specs$specifications$tea)
cat(sprintf("Sigma Metric: %.2f (%s)\n\n", sm$sigma, sm$interpretation$category))

# Step 4: Full assessment
assessment <- ate_assessment(
  bias = observed_bias,
  cv = observed_cv,
  tea = specs$specifications$tea,
  allowable_bias = specs$specifications$allowable_bias,
  allowable_cv = specs$specifications$allowable_cv
)

# Step 5: Decision
if (assessment$assessment$overall) {
  cat("DECISION: Method acceptable for clinical use\n")
} else {
  cat("DECISION: Method requires improvement\n")
}

## ----other-sources------------------------------------------------------------
# Using a CLIA-based TEa for glucose (example: ±6 mg/dL or ±10%)
# For a sample at 100 mg/dL, 10% = 10 mg/dL
sm_clia <- sigma_metric(bias = 2, cv = 3, tea = 10)
sm_clia

