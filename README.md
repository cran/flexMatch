# flexMatch

Flexible matching for stacked difference-in-differences panels with staggered or repeated treatment events.

```r
library(flexMatch)
data("flexmatch_example")
d <- flexmatch_example[flexmatch_example$unit_id <= 100, ]
r <- flexmatch(d, "unit_id", "year", "event",
  k_pre = 1, k_post = 1, k_freeze = 1, match_vars = c("size", "leverage"),
  method = "nearest", assess_quality = TRUE, plot_events = TRUE, verbose = FALSE)
r$quality_stats$balance_table
print(r$event_plot)
vignette("introduction", package = "flexMatch")
```

The three stages are `label_candidates()`, `match_candidates()` and `expand_matched_stacks()`. Matching is exact on the event cohort period, with optional exact covariates and additional MatchIt options. Each matched set has a `cohort_id`; each unit within the set has a `stack_id`, including reused controls. `assess_match_quality()` produces numeric and categorical balance plots.

The example data are entirely synthetic and contain 65,000 unit-year observations. The HTML vignette shows a small runnable example. Use `help(package = "flexMatch")` or `?flexmatch` for the R help index and reference pages.

The illustrative `fit_event_study()` uses matching weights, unit-within-set fixed effects, matched-set-specific calendar-period effects, and original-unit clustering. Sample screening and balance do not prove parallel trends or identify a causal effect. Repeated events require a justified cooldown assumption. Retaining full histories with `include_all_years = TRUE` also retains observations outside the screened window.

Time distances refer to the ordered periods observed in the input. Entirely missing calendar periods cannot be inferred. `control_type = "all"` and `"nt"` currently select the same clean-window pool.
