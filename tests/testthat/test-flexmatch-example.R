example_slice <- function() {
  data("flexmatch_example", package = "flexMatch")
  flexmatch_example[flexmatch_example$unit_id <= 300L, ]
}

test_that("flexmatch_example has the expected structure", {
  data("flexmatch_example", package = "flexMatch")

  expect_s3_class(flexmatch_example, "data.frame")
  expect_equal(nrow(flexmatch_example), 65000)
  expect_equal(
    names(flexmatch_example),
    c("unit_id", "year", "industry", "region", "size", "leverage", "event", "outcome")
  )
  expect_equal(length(unique(flexmatch_example$unit_id)), 10000)
  expect_equal(sort(unique(flexmatch_example$year)), 2013:2022)
  expect_true(all(flexmatch_example$event %in% c(0L, 1L)))
  obs_per_unit <- table(flexmatch_example$unit_id)
  expect_equal(min(obs_per_unit), 3)
  expect_equal(max(obs_per_unit), 10)
  expect_false(length(unique(obs_per_unit)) == 1)
})

test_that("flexmatch returns windowed stacks by default", {
  flexmatch_example <- example_slice()

  res <- suppressWarnings(flexmatch(
    dt = flexmatch_example,
    id_col = "unit_id",
    time_col = "year",
    event_col = "event",
    k_pre = 2,
    k_post = 2,
    k_freeze = 1,
    match_vars = c("size", "leverage"),
    method = "nearest",
    ratio = 1,
    assess_quality = FALSE,
    plot_events = FALSE,
    verbose = FALSE
  ))

  expect_true(is.list(res))
  expect_true(all(c(
    "final_expanded_data",
    "matched_cross_section",
    "stack_candidates",
    "diagnostics",
    "match_object"
  ) %in% names(res)))

  expect_equal(sort(unique(res$final_expanded_data$rel_time)), -2:2)
  expect_true(all(res$final_expanded_data$rel_time >= -2))
  expect_true(all(res$final_expanded_data$rel_time <= 2))
  expect_gt(nrow(res$matched_cross_section), 0)
  expect_gt(res$diagnostics$matching$treated_matched, 0)
})

test_that("include_all_years keeps full unit histories and clips tails", {
  flexmatch_example <- example_slice()

  res_window <- suppressWarnings(flexmatch(
    dt = flexmatch_example,
    id_col = "unit_id",
    time_col = "year",
    event_col = "event",
    k_pre = 2,
    k_post = 2,
    k_freeze = 1,
    match_vars = c("size", "leverage"),
    method = "nearest",
    ratio = 1,
    include_all_years = FALSE,
    assess_quality = FALSE,
    plot_events = FALSE,
    verbose = FALSE
  ))

  res_all_years <- suppressWarnings(flexmatch(
    dt = flexmatch_example,
    id_col = "unit_id",
    time_col = "year",
    event_col = "event",
    k_pre = 2,
    k_post = 2,
    k_freeze = 1,
    match_vars = c("size", "leverage"),
    method = "nearest",
    ratio = 1,
    include_all_years = TRUE,
    assess_quality = FALSE,
    plot_events = FALSE,
    verbose = FALSE
  ))

  expect_gt(nrow(res_all_years$final_expanded_data), nrow(res_window$final_expanded_data))
  expect_equal(sort(unique(res_window$final_expanded_data$rel_time)), -2:2)
  expect_equal(sort(unique(res_all_years$final_expanded_data$rel_time)), c(-3:-2, -1:2, 3))
  expect_true(any(res_all_years$final_expanded_data$rel_time == -3))
  expect_true(any(res_all_years$final_expanded_data$rel_time == 3))
})

test_that("treated_type all yields at least as many treated candidates as first", {
  flexmatch_example <- example_slice()

  cands_first <- label_candidates(
    dt = flexmatch_example,
    id_col = "unit_id",
    time_col = "year",
    event_col = "event",
    k_pre = 1,
    k_post = 1,
    k_freeze = 1,
    treated_type = "first",
    control_type = "all"
  )

  cands_all <- label_candidates(
    dt = flexmatch_example,
    id_col = "unit_id",
    time_col = "year",
    event_col = "event",
    k_pre = 1,
    k_post = 1,
    k_freeze = 1,
    treated_type = "all",
    control_type = "all"
  )

  expect_gte(nrow(cands_all[treat == 1]), nrow(cands_first[treat == 1]))
})

test_that("control_type nyt uses only eventually treated units before first treatment", {
  flexmatch_example <- example_slice()

  cands_nyt <- label_candidates(
    dt = flexmatch_example,
    id_col = "unit_id",
    time_col = "year",
    event_col = "event",
    k_pre = 1,
    k_post = 1,
    k_freeze = 1,
    treated_type = "first",
    control_type = "nyt"
  )

  control_ids <- unique(cands_nyt[treat == 0, unit_id])
  ever_treated_ids <- unique(flexmatch_example[event == 1, unit_id])
  first_treat_year <- flexmatch_example[event == 1, .(first_treat_year = min(year)), by = unit_id]
  nyt_controls <- merge(
    cands_nyt[treat == 0, .(unit_id, cohort_year)],
    first_treat_year,
    by = "unit_id",
    all.x = TRUE
  )

  expect_true(length(control_ids) > 0)
  expect_true(all(control_ids %in% ever_treated_ids))
  expect_true(all(nyt_controls$cohort_year < nyt_controls$first_treat_year))
})

test_that("example dataset contains never-treated, once-treated, and repeat-treated units", {
  data("flexmatch_example", package = "flexMatch")

  event_counts <- flexmatch_example[, .(n_events = sum(event)), by = unit_id]

  expect_true(any(event_counts$n_events == 0))
  expect_true(any(event_counts$n_events == 1))
  expect_true(any(event_counts$n_events >= 2))
})

test_that("fit_event_study runs on the example stacked output", {
  skip_if_not_installed("fixest")
  flexmatch_example <- example_slice()

  res <- suppressWarnings(flexmatch(
    dt = flexmatch_example,
    id_col = "unit_id",
    time_col = "year",
    event_col = "event",
    k_pre = 2,
    k_post = 2,
    k_freeze = 1,
    match_vars = c("size", "leverage"),
    method = "nearest",
    ratio = 1,
    assess_quality = FALSE,
    plot_events = FALSE,
    verbose = FALSE
  ))

  est <- fit_event_study(
    final_data = res$final_expanded_data,
    outcome = "outcome",
    time_col = "year"
  )

  expect_s3_class(est, "fixest")
})

test_that("plot_event_dynamics can shade observed data availability", {
  flexmatch_example <- data.table::as.data.table(example_slice())
  flexmatch_example[unit_id %% 7L == 0L & year >= 2019L, outcome := NA_real_]
  flexmatch_example[unit_id %% 11L == 0L & year <= 2016L, outcome := NA_real_]

  plot_out <- plot_event_dynamics(
    dt = flexmatch_example,
    id_col = "unit_id",
    time_col = "year",
    event_col = "event",
    show_data = "outcome",
    sort_by_intensity = TRUE
  )

  expect_s3_class(plot_out$plot, "ggplot")
  expect_true(all(c("first_event_idx", "data_start_idx", "data_end_idx") %in% names(plot_out$order_table)))
  expect_true(any(is.finite(plot_out$order_table$data_start_idx)))
  expect_match(plot_out$plot$labels$subtitle, "Gray: non-missing outcome", fixed = TRUE)
})

test_that("assess_match_quality returns upgraded balance plots", {
  flexmatch_example <- example_slice()

  res <- suppressWarnings(flexmatch(
    dt = flexmatch_example,
    id_col = "unit_id",
    time_col = "year",
    event_col = "event",
    k_pre = 2,
    k_post = 2,
    k_freeze = 1,
    match_vars = c("size", "leverage", "industry"),
    method = "nearest",
    ratio = 1,
    assess_quality = TRUE,
    plot_events = FALSE,
    verbose = FALSE
  ))

  expect_s3_class(res$quality_stats$numeric_balance_panel, "ggplot")
  expect_s3_class(res$quality_stats$categorical_balance_panel, "ggplot")
  expect_match(
    res$quality_stats$numeric_balance_panel$labels$subtitle,
    "Vertical lines mark group means",
    fixed = TRUE
  )
  expect_match(
    res$quality_stats$categorical_balance_panel$labels$subtitle,
    "Levels are ordered by overall frequency",
    fixed = TRUE
  )
})
