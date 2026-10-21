## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 4)
library(flexMatch)

## ----data---------------------------------------------------------------------
data("flexmatch_example")
panel <- flexmatch_example[flexmatch_example$unit_id <= 100, ]
head(panel)

## ----pipeline-----------------------------------------------------------------
result <- flexmatch(panel, "unit_id", "year", "event",
  k_pre = 1, k_post = 1, k_freeze = 1,
  match_vars = c("size", "leverage"), method = "nearest",
  assess_quality = TRUE, plot_events = TRUE, verbose = FALSE)
result$diagnostics
head(result$final_expanded_data)

## ----timing-------------------------------------------------------------------
print(result$event_plot)

## ----balance------------------------------------------------------------------
result$quality_stats$balance_table
print(result$quality_stats$numeric_balance_panel)

## ----estimation---------------------------------------------------------------
if (requireNamespace("fixest", quietly = TRUE)) {
  fit <- fit_event_study(result$final_expanded_data, "outcome", "year")
  coef(fit)
}

