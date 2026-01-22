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
data("glucose_methods")
head(glucose_methods)

## ----scatter-raw, fig.cap = "Scatter plot of POC vs laboratory glucose measurements with identity line."----
ggplot(glucose_methods, aes(x = reference, y = poc_meter)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    x = "Reference Method (mg/dL)",
    y = "POC Meter (mg/dL)",
    title = "Glucose Method Comparison"
  ) +
  coord_fixed() +
  theme_minimal()

## ----ba-analysis--------------------------------------------------------------
# Vector interface
ba <- ba_analysis(
  x = glucose_methods$reference,
  y = glucose_methods$poc_meter
)

# Alternative: formula interface
# ba <- ba_analysis(reference ~ poc_meter, data = glucose_methods)

ba

## ----ba-summary---------------------------------------------------------------
summary(ba)

## ----ba-plot, fig.cap = "Bland-Altman plot showing bias and 95% limits of agreement."----
plot(ba)

## ----ba-autoplot, fig.cap = "Customized Bland-Altman plot."-------------------
autoplot(ba) +
  labs(title = "POC Meter vs Reference Analyzer Agreement") +
  theme_bw()

## ----ba-percent---------------------------------------------------------------
ba_pct <- ba_analysis(
  x = glucose_methods$reference,
  y = glucose_methods$poc_meter,
  type = "percent"
)
ba_pct

## ----ba-percent-plot, fig.cap = "Bland-Altman plot with percentage differences."----
plot(ba_pct)

## ----pb-regression------------------------------------------------------------
pb <- pb_regression(
  x = glucose_methods$reference,
  y = glucose_methods$poc_meter
)
pb

## ----pb-summary---------------------------------------------------------------
summary(pb)

## ----pb-scatter, fig.cap = "Passing-Bablok regression with 95% confidence band."----
plot(pb, type = "scatter")

## ----pb-residuals, fig.cap = "Perpendicular residuals from Passing-Bablok regression."----
plot(pb, type = "residuals")

## ----pb-cusum, fig.cap = "CUSUM plot for linearity assessment."---------------
plot(pb, type = "cusum")

## ----pb-bootstrap-------------------------------------------------------------
pb_boot <- pb_regression(
  x = glucose_methods$reference,
  y = glucose_methods$poc_meter,
  ci_method = "bootstrap",
  boot_n = 1999
)
summary(pb_boot)

## ----workflow-summary, eval = FALSE-------------------------------------------
#  # 1. Load and inspect data
#  data("glucose_methods")
#  
#  # 2. Bland-Altman analysis for agreement assessment
#  ba <- ba_analysis(reference ~ poc_meter, data = glucose_methods)
#  summary(ba)
#  plot(ba)
#  
#  # 3. Passing-Bablok regression for systematic differences
#  pb <- pb_regression(reference ~ poc_meter, data = glucose_methods)
#  summary(pb)
#  plot(pb, type = "scatter")
#  plot(pb, type = "cusum")
#  
#  # 4. Document conclusions
#  # - Bias and LoA from Bland-Altman
#  # - Slope and intercept CIs from Passing-Bablok
#  # - Clinical interpretation based on acceptable performance criteria

## ----missing-data-------------------------------------------------------------
# Create data with missing values for demonstration
glucose_missing <- glucose_methods
glucose_missing$poc_meter[c(5, 15, 25)] <- NA

# Default behavior: remove pairs with missing values
ba_complete <- ba_analysis(
  reference ~ poc_meter,
  data = glucose_missing,
  na_action = "omit"
)

# Require complete cases (will error if any NA present)
# ba_strict <- ba_analysis(
#   reference ~ poc_meter,
#   data = glucose_missing,
#   na_action = "fail"
# )

## ----other-datasets-----------------------------------------------------------
# Creatinine: enzymatic vs Jaffe methods
data("creatinine_serum")
head(creatinine_serum)

# High-sensitivity troponin: two immunoassay platforms
data("troponin_cardiac")
head(troponin_cardiac)

