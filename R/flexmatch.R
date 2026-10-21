#' Flexible Matching for Stacked Difference-in-Differences
#'
#' @details
#' Start with [flexmatch()] and [flexmatch_example]. The pipeline screens
#' event windows, matches covariates using 'MatchIt', expands matched sets,
#' and produces balance and timing diagnostics. See
#' \code{vignette("introduction", package = "flexMatch")} for a worked example.
#' Matching and timing screens do not establish parallel trends, rule out
#' anticipation, or identify a causal effect by themselves. Repeated treatments
#' require a substantive justification for the chosen cooldown period.
#' @keywords internal
"_PACKAGE"


#' @import data.table
#' @importFrom stats na.omit
NULL

as_flex_dt <- function(x) {
  data.table::as.data.table(data.table::copy(x))
}

validate_panel_columns <- function(dt, cols, object_name = "dt") {
  missing_cols <- setdiff(cols, names(dt))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "Missing required column(s) in %s: %s",
        object_name,
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

validate_flexmatch_inputs <- function(dt,
                                      id_col,
                                      time_col,
                                      event_col = NULL,
                                      k_pre = NULL,
                                      k_post = NULL,
                                      k_freeze = NULL,
                                      match_vars = NULL,
                                      match_exact = NULL) {
  validate_panel_columns(dt, c(id_col, time_col, event_col, match_vars, match_exact))

  if (anyNA(dt[[id_col]]) || anyNA(dt[[time_col]])) {
    stop("Unit and time identifiers must not contain missing values.", call. = FALSE)
  }

  if (anyDuplicated(dt[, c(id_col, time_col), with = FALSE]) > 0) {
    stop("Each id/time combination must be unique.", call. = FALSE)
  }

  if (!is.null(event_col)) {
    if (anyNA(dt[[event_col]])) stop("Treatment indicators must not contain missing values.", call. = FALSE)
    event_vals <- unique(stats::na.omit(dt[[event_col]]))
    if (length(setdiff(event_vals, c(0, 1))) > 0) {
      stop(sprintf("`%s` must be a binary 0/1 event indicator.", event_col), call. = FALSE)
    }
  }

  int_args <- list(k_pre = k_pre, k_post = k_post, k_freeze = k_freeze)
  for (arg_name in names(int_args)) {
    arg_val <- int_args[[arg_name]]
    if (is.null(arg_val)) {
      next
    }
    if (!is.numeric(arg_val) || length(arg_val) != 1L || is.na(arg_val) || arg_val < 0 || arg_val %% 1 != 0) {
      stop(sprintf("`%s` must be a single non-negative integer.", arg_name), call. = FALSE)
    }
  }

  if (!is.null(match_vars) && length(match_vars) == 0L) {
    stop("`match_vars` must contain at least one covariate.", call. = FALSE)
  }
}

build_time_index_map <- function(dt, id_col, time_col) {
  time_map <- dt[, .SD, .SDcols = c(id_col, time_col)]
  setnames(time_map, c(id_col, time_col), c("id", "period"))
  setorderv(time_map, c("id", "period"))

  global_periods <- sort(unique(time_map$period))
  time_map[, time_idx := match(period, global_periods)]

  list(
    time_map = time_map,
    global_periods = global_periods
  )
}

get_label_safe <- function(v, df) {
  if (!v %in% names(df)) {
    return(v)
  }

  lbl <- attr(df[[v]], "label")
  if (is.null(lbl) || length(lbl) == 0) {
    return(as.character(v))
  }

  as.character(lbl)
}

build_balance_plot_data <- function(matched_df, unmatched_pool, covariates, var_labs) {
  dt_matched <- melt(
    matched_df[, c("treat", "weights", covariates), with = FALSE],
    measure.vars = covariates,
    variable.name = "term"
  )
  dt_matched[, Stage := "Matched"]

  dt_unmatched <- melt(
    unmatched_pool[, c("treat", covariates), with = FALSE],
    measure.vars = covariates,
    variable.name = "term"
  )
  dt_unmatched[, weights := 1]
  dt_unmatched[, Stage := "Unmatched"]

  dt_plot <- rbind(dt_matched, dt_unmatched, fill = TRUE)
  dt_plot[, Variable_Label := var_labs[as.character(term)]]
  dt_plot[, Group := factor(treat, levels = c(0, 1), labels = c("Control", "Treated"))]

  dt_plot
}

build_numeric_balance_plot <- function(dt_plot) {
  mean_lines <- dt_plot[
    !is.na(value),
    list(mean_value = stats::weighted.mean(value, w = weights, na.rm = TRUE)),
    by = list(Variable_Label, Stage, Group)
  ]

  ggplot2::ggplot() +
    ggplot2::geom_density(
      data = dt_plot[Stage == "Unmatched"],
      ggplot2::aes(x = value, color = Group, linetype = "Unmatched"),
      linewidth = 0.6,
      alpha = 0.6
    ) +
    ggplot2::geom_density(
      data = dt_plot[Stage == "Matched"],
      ggplot2::aes(x = value, color = Group, weight = weights, linetype = "Matched"),
      linewidth = 1
    ) +
    ggplot2::geom_vline(
      data = mean_lines,
      ggplot2::aes(xintercept = mean_value, color = Group, linetype = Stage),
      linewidth = 0.5,
      alpha = 0.7
    ) +
    ggplot2::facet_wrap(~ Variable_Label, scales = "free", ncol = 3) +
    ggplot2::scale_linetype_manual(name = "Stage", values = c("Unmatched" = "dashed", "Matched" = "solid")) +
    ggplot2::scale_color_manual(values = c("Control" = "#377EB8", "Treated" = "#E41A1C")) +
    ggplot2::labs(
      title = "Covariate Balance: Full Pool (Dashed) vs Matched Set (Solid)",
      subtitle = "Solid lines should overlap. Vertical lines mark group means before and after matching.",
      x = "Covariate Value",
      y = "Density"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top", strip.text = ggplot2::element_text(face = "bold"))
}

build_categorical_balance_plot <- function(dt_plot) {
  dt_plot[, value_chr := as.character(value)]
  dt_counts <- dt_plot[
    ,
    list(weight = sum(weights, na.rm = TRUE)),
    by = list(Variable_Label, term, Stage, Group, value_chr)
  ]
  dt_counts[, share := weight / sum(weight), by = list(Variable_Label, Stage, Group)]
  level_order <- dt_counts[
    ,
    list(total_weight = sum(weight, na.rm = TRUE)),
    by = list(Variable_Label, term, value_chr)
  ]
  setorderv(level_order, c("Variable_Label", "total_weight", "value_chr"), c(1, -1, 1))
  level_order[, level_rank := seq_len(.N), by = list(Variable_Label, term)]
  dt_counts <- merge(
    dt_counts,
    level_order[, list(Variable_Label, term, value_chr, level_rank)],
    by = c("Variable_Label", "term", "value_chr"),
    all.x = TRUE,
    sort = FALSE
  )
  setorderv(dt_counts, c("Variable_Label", "level_rank", "value_chr"))
  dt_counts[, value_chr := factor(value_chr, levels = unique(value_chr)), by = list(Variable_Label, term)]

  ggplot2::ggplot(
    dt_counts,
    ggplot2::aes(x = value_chr, y = share, fill = Group)
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_grid(Stage ~ Variable_Label, scales = "free_x", space = "free_x") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      title = "Categorical Covariate Balance",
      subtitle = "Share of each level by treatment group before and after matching. Levels are ordered by overall frequency.",
      x = "Level",
      y = "Share"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "top",
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

format_match_diagnostics <- function(stack_candidates, matched_cross_section, id_col, dropped_rows_missing_covariates) {
  treated_candidate_ids <- unique(stack_candidates[treat == 1, get(id_col)])
  matched_treated_ids <- unique(matched_cross_section[treat == 1, get(id_col)])
  unmatched_treated_ids <- setdiff(treated_candidate_ids, matched_treated_ids)

  list(
    candidate_counts = list(
      treated = nrow(stack_candidates[treat == 1]),
      control = nrow(stack_candidates[treat == 0])
    ),
    matching = list(
      treated_candidates = length(treated_candidate_ids),
      treated_matched = length(matched_treated_ids),
      treated_unmatched = length(unmatched_treated_ids),
      unmatched_treated_ids = unmatched_treated_ids,
      dropped_rows_missing_covariates = dropped_rows_missing_covariates
    )
  )
}

# ==============================================================================
# 1. SCAN: Label Candidates
# ==============================================================================
#' Label Treatment and Control Stack Candidates
#'
#' @description
#' Identifies valid stack candidates while respecting treatment timing,
#' clean windows, and exact matching period lookup.
#'
#' @inheritParams flexmatch
#' @return A data.table with eligible unit/cohort rows, a treatment flag, and
#'   the period at which matching covariates are observed.
#' @details Period distances use the ordered periods present anywhere in the
#'   panel. Missing calendar periods across the entire input cannot be inferred.
#'   Controls must have a clean full window; treated units must have clean
#'   pre-event observations and no additional event inside the post-event window.
#'   All unit/time pairs must be unique.
#' @examples
#' data("flexmatch_example")
#' d <- flexmatch_example[flexmatch_example$unit_id <= 100, ]
#' candidates <- label_candidates(d, "unit_id", "year", "event", 1, 1, 1)
#' table(candidates$treat)
#' @export
label_candidates <- function(dt,
                             id_col,
                             time_col,
                             event_col,
                             k_pre,
                             k_post,
                             k_freeze,
                             match_by = -1,
                             treated_type = c("first", "all"),
                             control_type = c("all", "nyt", "nt", "nevert")) {

  if (!is.numeric(match_by) || length(match_by) != 1L || is.na(match_by) ||
      !is.finite(match_by) || match_by %% 1 != 0 || match_by > 0) {
    stop("match_by must be a single non-positive integer.", call. = FALSE)
  }
  treated_type <- match.arg(treated_type)
  control_type <- match.arg(control_type)

  dt <- as_flex_dt(dt)
  validate_flexmatch_inputs(
    dt = dt,
    id_col = id_col,
    time_col = time_col,
    event_col = event_col,
    k_pre = k_pre,
    k_post = k_post,
    k_freeze = k_freeze
  )

  calc_dt <- dt[, .SD, .SDcols = c(id_col, time_col, event_col)]
  setnames(calc_dt, c(id_col, time_col, event_col), c("id", "period", "event"))
  setorderv(calc_dt, c("id", "period"))

  time_info <- build_time_index_map(calc_dt, "id", "period")
  calc_dt <- merge(calc_dt, time_info$time_map, by = c("id", "period"), all.x = TRUE, sort = FALSE)
  setorderv(calc_dt, c("id", "time_idx"))

  global_stats <- calc_dt[, list(
    is_ever_treated = max(event, na.rm = TRUE) == 1,
    first_treat_idx = if (max(event, na.rm = TRUE) == 1) min(time_idx[event == 1], na.rm = TRUE) else NA_integer_
  ), by = id]
  calc_dt <- merge(calc_dt, global_stats, by = "id", all.x = TRUE, sort = FALSE)

  calc_dt[, evt_loc := fifelse(event == 1, time_idx, NA_integer_)]
  calc_dt[, last_evt := nafill(evt_loc, type = "locf"), by = id]
  calc_dt[, next_evt := nafill(evt_loc, type = "nocb"), by = id]
  calc_dt[, dist_last := time_idx - last_evt]
  calc_dt[, dist_next := next_evt - time_idx]

  calc_dt[, is_dirty_past := !is.na(dist_last) & dist_last >= 0 & dist_last <= (k_post + k_freeze)]
  calc_dt[, is_dirty_future := !is.na(dist_next) & dist_next <= k_post]

  calc_dt[, lag_idx := shift(time_idx, n = k_pre, type = "lag"), by = id]
  calc_dt[, lead_idx := shift(time_idx, n = k_post, type = "lead"), by = id]
  calc_dt[, has_history := !is.na(lag_idx) & (time_idx - lag_idx == k_pre)]
  calc_dt[, has_future := !is.na(lead_idx) & (lead_idx - time_idx == k_post)]

  calc_dt[, is_dirty_past_int := as.integer(is_dirty_past)]
  if (k_pre == 0L) {
    calc_dt[, sum_dirty_pre := 0L]
  } else {
    calc_dt[, sum_dirty_pre := shift(
      frollsum(is_dirty_past_int, n = k_pre, align = "right", fill = NA), 1L
    ), by = id]
  }
  calc_dt[, next_treatment := shift(next_evt, type = "lead"), by = id]
  calc_dt[, is_valid_treated := (event == 1) & has_history & has_future &
    (sum_dirty_pre == 0) & (is.na(next_treatment) | next_treatment - time_idx > k_post)]

  calc_dt[, is_dirty_any := is_dirty_past | is_dirty_future]
  calc_dt[, is_dirty_any_int := as.integer(is_dirty_any)]
  win_len <- k_pre + 1L + k_post
  calc_dt[, sum_dirty_full := shift(
    frollsum(is_dirty_any_int, n = win_len, align = "right", fill = NA),
    n = k_post, type = "lead"
  ), by = id]

  calc_dt[, is_valid_control := has_history & has_future & (sum_dirty_full == 0)]

  treats <- calc_dt[
    is_valid_treated == TRUE,
    list(id, cohort_year = period, time_idx, treat = 1, first_treat_idx)
  ]

  if (treated_type == "first") {
    treats <- treats[time_idx == first_treat_idx]
    treats <- unique(treats, by = "id")
  }
  treats[, first_treat_idx := NULL]

  ctrls_raw <- calc_dt[
    is_valid_control == TRUE,
    list(id, cohort_year = period, time_idx, treat = 0, is_ever_treated, first_treat_idx)
  ]

  if (control_type == "nevert") {
    ctrls <- ctrls_raw[is_ever_treated == FALSE]
  } else if (control_type == "nt") {
    ctrls <- ctrls_raw
  } else if (control_type == "nyt") {
    ctrls <- ctrls_raw[(is_ever_treated == TRUE) & (time_idx < first_treat_idx)]
  } else {
    ctrls <- ctrls_raw
  }
  ctrls[, c("is_ever_treated", "first_treat_idx") := NULL]

  candidates <- rbind(treats, ctrls, fill = TRUE)
  candidates[, match_time_idx := time_idx + match_by]

  time_map <- unique(calc_dt[, list(id, time_idx, real_match_period = period)])
  candidates <- merge(
    candidates,
    time_map,
    by.x = c("id", "match_time_idx"),
    by.y = c("id", "time_idx"),
    all.x = FALSE,
    sort = FALSE
  )

  candidates[, match_period := real_match_period]
  candidates[, c("match_time_idx", "real_match_period", "time_idx") := NULL]
  setnames(candidates, "id", id_col)

  candidates[]
}

# ==============================================================================
# 2. MATCH: Cross-Sectional Matching
# ==============================================================================
#' Cross-Sectional Matching for Stacked DiD
#'
#' @description
#' Performs exact cohort matching plus flexible 1-to-K matching via \code{MatchIt}.
#'
#' @inheritParams flexmatch
#' @param stack_candidates Candidate table returned by [label_candidates()].
#' @return A list with matched_cross_section, the original MatchIt match_object,
#'   and diagnostics for missing matching covariates.
#' @details Matching is exact on the event cohort period, in addition to any
#'   exact covariates supplied. Matching with replacement can reuse a unit in
#'   different matched sets. Rows with missing matching covariates are removed.
#'   The methods 'cem' and 'subclass' are rejected because their pooled calls
#'   ignore the required exact-cohort argument. For method = 'exact', the
#'   cohort and exact covariates are included directly in the matching formula.
#' @examples
#' data("flexmatch_example")
#' d <- flexmatch_example[flexmatch_example$unit_id <= 100, ]
#' c <- label_candidates(d, "unit_id", "year", "event", 1, 1, 1)
#' m <- match_candidates(d, c, "unit_id", "year", "size", method = "nearest")
#' head(m$matched_cross_section)
#' @export
match_candidates <- function(dt, stack_candidates, id_col, time_col, match_vars, match_exact = NULL, ...) {
  dt <- as_flex_dt(dt)
  stack_candidates <- as_flex_dt(stack_candidates)

  validate_flexmatch_inputs(
    dt = dt,
    id_col = id_col,
    time_col = time_col,
    match_vars = match_vars,
    match_exact = match_exact
  )
  validate_panel_columns(stack_candidates, c(id_col, "match_period", "cohort_year", "treat"), "stack_candidates")

  cols_to_fetch <- unique(c(id_col, time_col, match_vars, match_exact))
  cov_data <- dt[, .SD, .SDcols = cols_to_fetch]
  setnames(cov_data, c(id_col, time_col), c("id", "match_period"))

  match_input_base <- copy(stack_candidates)
  setnames(match_input_base, id_col, "id")

  match_input <- merge(
    match_input_base,
    cov_data,
    by = c("id", "match_period"),
    all.x = FALSE,
    sort = FALSE
  )

  before_complete_cases <- nrow(match_input)
  match_input <- na.omit(match_input, cols = c(match_vars, match_exact))
  dropped_rows_missing_covariates <- before_complete_cases - nrow(match_input)

  if (nrow(match_input[treat == 1]) == 0) {
    stop("No treated candidates remain after dropping rows with missing matching covariates.", call. = FALSE)
  }
  if (nrow(match_input[treat == 0]) == 0) {
    stop("No control candidates remain after dropping rows with missing matching covariates.", call. = FALSE)
  }

  fm <- stats::as.formula(paste("treat ~", paste(match_vars, collapse = " + ")))
  args <- list(...)
  args$formula <- fm
  args$data <- match_input

  if (!is.null(match_exact)) {
    args$exact <- stats::as.formula(paste("~", paste(c(match_exact, "cohort_year"), collapse = " + ")))
  } else if ("exact" %in% names(args)) {
    cur <- args$exact
    if (inherits(cur, "formula")) {
      args$exact <- stats::as.formula(paste("~", as.character(cur)[2], "+ cohort_year"))
    } else {
      args$exact <- stats::as.formula(paste("~", paste(c(cur, "cohort_year"), collapse = " + ")))
    }
  } else {
    args$exact <- ~ cohort_year
  }

  if (identical(args$method, "exact")) {
    # Exact matching uses formula covariates rather than the ignored exact
    # argument. Include cohort membership directly in that formula.
    args$formula <- stats::as.formula(paste("treat ~",
      paste(unique(c(match_vars, match_exact, "cohort_year")), collapse = " + ")))
    args$exact <- NULL
  }
  if (!is.null(args$method) && args$method %in% c("cem", "subclass")) {
    stop("This matching method does not enforce the required exact cohort constraint; use a method with exact matching support or method = 'exact'.", call. = FALSE)
  }
  m_out <- do.call("matchit", args, envir = asNamespace("MatchIt"))

  matched_data <- if (!is.null(m_out$match.matrix)) {
    as.data.table(MatchIt::get_matches(m_out, id = "matchit_row_id"))
  } else {
    as.data.table(MatchIt::match_data(m_out))
  }
  if (!"subclass" %in% names(matched_data)) matched_data[, subclass := as.character(cohort_year)]
  if (any(matched_data[, uniqueN(cohort_year), by = subclass]$V1 > 1L)) {
    stop("A matched set contains multiple cohort periods; exact cohort matching is required.", call. = FALSE)
  }
  # A MatchIt subclass is one matched set, including methods with multiple
  # treated units. One identifier per set avoids a many-to-many expansion.
  matched_data[, cohort_id := paste(cohort_year, subclass, sep = "_")]
  matched_data[, stack_id := paste(id, cohort_id, sep = "_")]
  setnames(matched_data, "id", id_col)

  cols_to_keep <- c("stack_id", "cohort_id", id_col, "cohort_year", "treat", "match_period", "weights", match_vars, match_exact)
  cols_to_keep <- intersect(cols_to_keep, names(matched_data))
  final_output <- matched_data[, .SD, .SDcols = cols_to_keep]

  list(
    matched_cross_section = final_output,
    match_object = m_out,
    diagnostics = list(
      dropped_rows_missing_covariates = dropped_rows_missing_covariates,
      input_rows = before_complete_cases,
      complete_case_rows = nrow(match_input)
    )
  )
}

# ==============================================================================
# 3. BUILD: Expand Stacks
# ==============================================================================
#' Expand Matched Cross-Section to Full Stacks
#'
#' @param include_all_years Logical. If \code{FALSE}, keep only periods in the
#' standard event window from \code{T-k_pre} through \code{T+k_post}. If
#' \code{TRUE}, keep every observed period for each matched unit and collapse
#' periods before the window to \code{-k_pre-1} and after the window to
#' \code{k_post+1}.
#'
#' @inheritParams flexmatch
#' @param matched_data Matched cross section returned by [match_candidates()].
#' @return A data.table with one row per matched unit, matched set, and retained
#'   period. It includes stack_id, cohort_id, rel_time, treat, post, and weights.
#'   The original unit column name is stored as an attribute for clustering.
#' @examples
#' data("flexmatch_example")
#' d <- flexmatch_example[flexmatch_example$unit_id <= 100, ]
#' c <- label_candidates(d, "unit_id", "year", "event", 1, 1, 1)
#' m <- match_candidates(d, c, "unit_id", "year", "size", method = "nearest")
#' s <- expand_matched_stacks(d, m$matched_cross_section, "unit_id", "year", 1, 1)
#' table(s$rel_time)
#' @export
expand_matched_stacks <- function(dt,
                                  matched_data,
                                  id_col,
                                  time_col,
                                  k_pre,
                                  k_post,
                                  match_by = -1,
                                  include_all_years = FALSE) {
  dt <- as_flex_dt(dt)
  matched_data <- as_flex_dt(matched_data)

  validate_flexmatch_inputs(dt = dt, id_col = id_col, time_col = time_col, k_pre = k_pre, k_post = k_post)
  validate_panel_columns(matched_data, c(id_col, "cohort_year", "treat"), "matched_data")

  time_info <- build_time_index_map(dt, id_col, time_col)
  time_map <- time_info$time_map

  matched_base <- copy(matched_data)
  setnames(matched_base, id_col, "id")

  matched_base <- merge(
    matched_base,
    time_map[, list(id, period, cohort_idx = time_idx)],
    by.x = c("id", "cohort_year"),
    by.y = c("id", "period"),
    all.x = TRUE,
    sort = FALSE
  )

  meta_cols <- c("stack_id", "cohort_id", "id", "cohort_idx", "treat", "weights", "cohort_year")
  meta_cols <- intersect(meta_cols, names(matched_base))
  stack_defs <- unique(matched_base[, ..meta_cols])

  if (isTRUE(include_all_years)) {
    expanded_stacks <- merge(
      stack_defs,
      time_map[, list(id, period, time_idx)],
      by = "id",
      all.x = TRUE,
      sort = FALSE,
      allow.cartesian = TRUE
    )
    expanded_stacks[, rel_time_raw := time_idx - cohort_idx]
    expanded_stacks[
      ,
      rel_time := fifelse(
        rel_time_raw < -k_pre,
        -k_pre - 1L,
        fifelse(rel_time_raw > k_post, k_post + 1L, rel_time_raw)
      )
    ]
    expanded_stacks[, rel_time_raw := NULL]
  } else {
    rel_times <- seq(-k_pre, k_post)
    n_stacks <- nrow(stack_defs)
    expanded_stacks <- stack_defs[rep(seq_len(n_stacks), each = length(rel_times))]
    expanded_stacks[, rel_time := rep(rel_times, times = n_stacks)]
    expanded_stacks[, target_idx := cohort_idx + rel_time]

    expanded_stacks <- merge(
      expanded_stacks,
      time_map[, list(id, target_idx = time_idx, period)],
      by = c("id", "target_idx"),
      all.x = TRUE,
      sort = FALSE
    )
  }

  final_df <- merge(
    expanded_stacks,
    dt,
    by.x = c("id", "period"),
    by.y = c(id_col, time_col),
    all.x = TRUE,
    sort = FALSE
  )

  final_df[, post := as.integer(rel_time >= 0)]
  cols_to_drop <- intersect(c("cohort_idx", "target_idx", "time_idx"), names(final_df))
  final_df[, (cols_to_drop) := NULL]
  setnames(final_df, "period", time_col)
  setnames(final_df, "id", id_col)

  first_cols <- c("stack_id", "cohort_id", id_col, time_col, "rel_time", "treat", "post", "weights")
  rest_cols <- setdiff(names(final_df), first_cols)
  setcolorder(final_df, c(first_cols, rest_cols))

  attr(final_df, "flexmatch_id_col") <- id_col
  final_df[]
}

# ==============================================================================
# 4. AUDIT: Quality Assessment
# ==============================================================================
#' Assess Matching Quality
#'
#' @inheritParams flexmatch
#' @param match_res List returned by [match_candidates()].
#' @param stack_candidates Candidate table returned by [label_candidates()].
#' @return A list with balance_table and optional numeric_balance_panel and
#'   categorical_balance_panel ggplot objects. Threshold labels in the table
#'   are descriptive heuristics, not statistical tests of exchangeability.
#' @examples
#' data("flexmatch_example")
#' d <- flexmatch_example[flexmatch_example$unit_id <= 100, ]
#' c <- label_candidates(d, "unit_id", "year", "event", 1, 1, 1)
#' m <- match_candidates(d, c, "unit_id", "year", "size", method = "nearest")
#' q <- assess_match_quality(m, c, d, "size", "unit_id", "year")
#' q$balance_table
#' print(q$numeric_balance_panel)
#' @export
assess_match_quality <- function(match_res, stack_candidates, dt, match_vars, id_col, time_col) {
  dt <- as_flex_dt(dt)
  stack_candidates <- as_flex_dt(stack_candidates)
  matched_df <- as_flex_dt(match_res$matched_cross_section)

  covariates <- intersect(match_vars, names(matched_df))
  if (length(covariates) == 0) {
    stop("None of the `match_vars` were found in matched results.", call. = FALSE)
  }

  cols_fetch <- c(id_col, time_col, covariates)
  dt_subset <- dt[, ..cols_fetch]

  sc_copy <- copy(stack_candidates)
  setnames(sc_copy, id_col, "id")

  dt_sub_copy <- copy(dt_subset)
  setnames(dt_sub_copy, c(id_col, time_col), c("id", "match_period"))

  unmatched_pool <- merge(sc_copy, dt_sub_copy, by = c("id", "match_period"), all.x = TRUE, sort = FALSE)
  unmatched_pool <- na.omit(unmatched_pool, cols = covariates)

  var_labs <- vapply(covariates, get_label_safe, df = matched_df, FUN.VALUE = character(1))
  names(var_labs) <- covariates

  numeric_covariates <- covariates[vapply(matched_df[, ..covariates], is.numeric, FUN.VALUE = logical(1))]
  categorical_covariates <- setdiff(covariates, numeric_covariates)

  numeric_panel <- NULL
  if (length(numeric_covariates) > 0) {
    numeric_plot_data <- build_balance_plot_data(matched_df, unmatched_pool, numeric_covariates, var_labs)
    numeric_panel <- build_numeric_balance_plot(numeric_plot_data)
  }

  categorical_panel <- NULL
  if (length(categorical_covariates) > 0) {
    categorical_plot_data <- build_balance_plot_data(matched_df, unmatched_pool, categorical_covariates, var_labs)
    categorical_panel <- build_categorical_balance_plot(categorical_plot_data)
  }

  s_out <- summary(match_res$match_object, standardize = TRUE)
  stats <- as.data.table(s_out$sum.matched, keep.rownames = "term")
  stats <- stats[term %in% covariates]

  clean_stats <- stats[, list(
    Variable = var_labs[term],
    Mean_T = `Means Treated`,
    Mean_C = `Means Control`,
    SMD = `Std. Mean Diff.`,
    VR = `Var. Ratio`
  )]
  clean_stats[, Quality_SMD := fifelse(abs(SMD) < 0.1, "Excellent", fifelse(abs(SMD) < 0.25, "Acceptable", "Poor"))]
  clean_stats[, Quality_VR := fifelse(is.na(VR), "N/A", fifelse(VR >= 0.8 & VR <= 1.25, "Balanced", "Distorted"))]

  list(
    balance_table = clean_stats,
    numeric_balance_panel = numeric_panel,
    categorical_balance_panel = categorical_panel
  )
}

# ==============================================================================
# 5. VISUALIZE: Event Dynamics
# ==============================================================================
#' Plot Event Dynamics (Lasagna Plot with Overlay)
#'
#' @param show_data Optional character scalar naming a variable in \code{dt}.
#'   When supplied, periods with non-missing data for that variable are shown as
#'   a light gray background layer, and unit ordering also uses the first and
#'   last periods with observed data before treatment-intensity tie-breaks.
#'
#' @inheritParams flexmatch
#' @param dt_matched Optional matched cross section used to mark selected units.
#' @param n_bins Positive integer number of treatment-intensity categories.
#' @param sort_by_intensity Logical; order units by treatment intensity.
#' @return A list containing a ggplot object in plot and the unit order_table.
#' @examples
#' data("flexmatch_example")
#' d <- flexmatch_example[flexmatch_example$unit_id <= 50, ]
#' p <- plot_event_dynamics(d, id_col = "unit_id", time_col = "year",
#'   event_col = "event", show_data = "outcome")
#' print(p$plot)
#' @export
plot_event_dynamics <- function(dt,
                                dt_matched = NULL,
                                id_col,
                                time_col,
                                event_col,
                                n_bins = 4L,
                                sort_by_intensity = FALSE,
                                show_data = NULL) {
  dt <- as_flex_dt(dt)
  validate_flexmatch_inputs(dt = dt, id_col = id_col, time_col = time_col, event_col = event_col)
  if (!is.null(show_data)) {
    if (!is.character(show_data) || length(show_data) != 1L || is.na(show_data)) {
      stop("show_data must be NULL or a single column name.", call. = FALSE)
    }
    validate_panel_columns(dt, show_data, "dt")
  }

  plot_cols <- c(id_col, time_col, event_col)
  if (!is.null(show_data)) {
    plot_cols <- c(plot_cols, show_data)
  }
  d <- dt[, .SD, .SDcols = plot_cols]
  setnames(d, old = c(id_col, time_col, event_col), new = c("id", "time", "event"))
  if (!is.null(show_data)) {
    setnames(d, old = show_data, new = "show_data_value")
  }

  t_levels <- sort(unique(d$time))

  id_stats <- d[, {
    tr_times <- time[event == 1L]
    tr_idx <- match(tr_times, t_levels)
    list(
      n_event = length(tr_idx),
      first_event_idx = as.numeric(if (length(tr_idx) > 0) min(tr_idx) else Inf),
      tr_idx = list(tr_idx)
    )
  }, by = id]

  id_stats[, is_treated := n_event > 0]

  if (!is.null(show_data)) {
    data_coverage <- d[!is.na(show_data_value), list(
      data_start_idx = as.numeric(min(match(time, t_levels))),
      data_end_idx = as.numeric(max(match(time, t_levels)))
    ), by = id]
    id_stats <- merge(id_stats, data_coverage, by = "id", all.x = TRUE, sort = FALSE)
    id_stats[is.na(data_start_idx), data_start_idx := Inf]
    id_stats[is.na(data_end_idx), data_end_idx := -Inf]
  }

  if (isTRUE(sort_by_intensity)) {
    treated_mask <- id_stats$is_treated == TRUE
    if (sum(treated_mask) > 0) {
      qs <- unique(stats::quantile(id_stats$n_event[treated_mask], probs = seq(0, 1, length.out = n_bins + 1L), type = 1))
      if (length(qs) < 2) {
        id_stats[treated_mask, bin := 1L]
      } else {
        id_stats[treated_mask, bin := cut(n_event, breaks = qs, include.lowest = TRUE, labels = FALSE)]
      }
    }
    id_stats[is.na(bin), bin := 0L]
  }

  K <- max(id_stats$n_event, na.rm = TRUE)
  if (!is.finite(K) || K < 1L) {
    K <- 1L
  }

  pad_event_index <- function(v) {
    out <- rep(-1L, K)
    if (length(v) > 0) {
      out[seq_along(v)] <- v
    }
    as.integer(out)
  }

  padded_list <- lapply(id_stats$tr_idx, pad_event_index)
  keys <- as.data.table(do.call(rbind, padded_list))
  time_key_cols <- paste0("t", seq_len(K))
  setnames(keys, time_key_cols)
  id_stats <- cbind(id_stats, keys)

  sort_order <- c("first_event_idx")
  order_directions <- c(1)
  if (!is.null(show_data)) {
    sort_order <- c(sort_order, "data_start_idx", "data_end_idx")
    order_directions <- c(order_directions, 1, -1)
  }
  if (isTRUE(sort_by_intensity)) {
    sort_order <- c(sort_order, "bin")
    order_directions <- c(order_directions, 1)
  }
  sort_order <- c(sort_order, "n_event", time_key_cols, "id")
  order_directions <- c(order_directions, -1, rep(-1, length(time_key_cols)), 1)
  setorderv(id_stats, cols = sort_order, order = order_directions)
  id_levels_sorted <- id_stats$id

  subtitle_bits <- "Green: Events"
  if (!is.null(show_data)) {
    subtitle_bits <- c(sprintf("Gray: non-missing %s", show_data), subtitle_bits)
  }
  if (!is.null(dt_matched)) {
    subtitle_bits <- c(subtitle_bits, "Red: Matched Treated", "Blue: Matched Control Frequency")
  }

  p <- ggplot2::ggplot() +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::labs(
      title = "Event Dynamics & Matched Stacks",
      subtitle = paste(subtitle_bits, collapse = " | "),
      x = time_col,
      y = id_col
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  if (!is.null(show_data)) {
    plot_data_avail <- d[!is.na(show_data_value)]
    if (nrow(plot_data_avail) > 0) {
      plot_data_avail[, id_factor := factor(id, levels = id_levels_sorted)]
      plot_data_avail[, time_factor := factor(time, levels = t_levels)]
      p <- p + ggplot2::geom_raster(
        data = plot_data_avail,
        ggplot2::aes(x = time_factor, y = id_factor),
        fill = "grey40",
        alpha = 0.12
      )
    }
  }

  plot_data_base <- d[event == 1L]
  if (nrow(plot_data_base) > 0) {
    plot_data_base[, id_factor := factor(id, levels = id_levels_sorted)]
    plot_data_base[, time_factor := factor(time, levels = t_levels)]
    p <- p + ggplot2::geom_raster(data = plot_data_base, ggplot2::aes(x = time_factor, y = id_factor), fill = "green4")
  }

  if (!is.null(dt_matched)) {
    dt_matched <- as_flex_dt(dt_matched)
    use_id <- if (id_col %in% names(dt_matched)) id_col else "id"
    use_time <- if (time_col %in% names(dt_matched)) time_col else "period"

    mat_trt <- dt_matched[treat == 1, .N, by = c(use_id, use_time)]
    setnames(mat_trt, c("id", "time", "N"))
    mat_trt <- mat_trt[id %in% id_levels_sorted & time %in% t_levels]
    if (nrow(mat_trt) > 0) {
      mat_trt[, id_factor := factor(id, levels = id_levels_sorted)]
      mat_trt[, time_factor := factor(time, levels = t_levels)]
      p <- p + ggplot2::geom_raster(data = mat_trt, ggplot2::aes(x = time_factor, y = id_factor), fill = "firebrick")
    }

    mat_ctrl <- dt_matched[treat == 0, .N, by = c(use_id, use_time)]
    setnames(mat_ctrl, c("id", "time", "N"))
    mat_ctrl <- mat_ctrl[id %in% id_levels_sorted & time %in% t_levels]
    if (nrow(mat_ctrl) > 0) {
      mat_ctrl[, id_factor := factor(id, levels = id_levels_sorted)]
      mat_ctrl[, time_factor := factor(time, levels = t_levels)]
      p <- p + ggplot2::geom_raster(data = mat_ctrl, ggplot2::aes(x = time_factor, y = id_factor, alpha = N), fill = "steelblue")
      p <- p + ggplot2::scale_alpha_continuous(range = c(0.3, 1), name = "Ctrl Freq.")
    }
  }

  list(plot = p, order_table = id_stats)
}

# ==============================================================================
# 6. ESTIMATE: Event-Study Wrapper
# ==============================================================================
#' Fit an Event-Study Model with fixest
#'
#' @description
#' Convenience wrapper around \code{fixest::feols()} for the stacked dataset
#' returned by \code{flexmatch()}.
#'
#' @param final_data Analysis-ready stacked data from \code{flexmatch()}.
#' @param outcome Character scalar naming the outcome column.
#' @param time_col Character scalar naming the calendar-time fixed effect column.
#' @param ref Reference period for \code{i(rel_time, treat)}. Default is \code{-1}.
#' @param cluster Character vector of clustering variables. By default, cluster
#'   by the original unit identifier stored by [expand_matched_stacks()].
#'   Supply this explicitly for other data.
#' @param ... Additional arguments passed to \code{fixest::feols()}.
#'
#' @return A fixest model with matching weights, unit-within-set fixed effects,
#'   matched-set-by-period fixed effects, and clustered standard errors.
#' @details This is an illustrative regression wrapper. Identification and
#'   the interpretation of pooled event coefficients require appropriate
#'   parallel-trends, anticipation, treatment-effect, and weighting assumptions.
#'   Repeated units are clustered together by default. Pass estimator options
#'   through ...; do not supply a second weights argument.
#' @examples
#' if (requireNamespace("fixest", quietly = TRUE)) {
#'   data("flexmatch_example")
#'   d <- flexmatch_example[flexmatch_example$unit_id <= 100, ]
#'   r <- flexmatch(d, "unit_id", "year", "event", 1, 1, 1, "size",
#'     method = "nearest", assess_quality = FALSE, plot_events = FALSE,
#'     verbose = FALSE)
#'   fit <- fit_event_study(r$final_expanded_data, "outcome", "year")
#'   coef(fit)
#' }
#' @export
fit_event_study <- function(final_data, outcome, time_col, ref = -1, cluster = NULL, ...) {
  if (!requireNamespace("fixest", quietly = TRUE)) {
    stop("Package `fixest` must be installed to use `fit_event_study()`.", call. = FALSE)
  }

  original_id <- attr(final_data, "flexmatch_id_col")
  if (is.null(cluster)) {
    if (is.null(original_id)) stop("Supply cluster for data not created by expand_matched_stacks().", call. = FALSE)
    cluster <- original_id
  }
  final_data <- as_flex_dt(final_data)
  validate_panel_columns(final_data, c("stack_id", "cohort_id", "rel_time", "treat", outcome, time_col, cluster), "final_data")

  model_formula <- stats::as.formula(
    paste0(outcome, " ~ fixest::i(rel_time, treat, ref = ", ref, ") | stack_id + cohort_id^", time_col)
  )
  cluster_formula <- if (!is.null(cluster)) stats::as.formula(paste("~", paste(cluster, collapse = " + "))) else NULL

  fixest::feols(
    fml = model_formula,
    data = final_data,
    cluster = cluster_formula,
    weights = if ("weights" %in% names(final_data)) ~ weights else NULL,
    ...
  )
}

# ==============================================================================
# 7. MASTER WRAPPER: flexmatch()
# ==============================================================================
#' Flexmatch: Automated Stacked Difference-in-Differences with Flexible Matching Logic
#'
#' @description
#' \code{flexmatch} orchestrates a complete Stacked Difference-in-Differences (DiD) pipeline.
#' The package is designed for staggered-adoption panel settings where users want to:
#' \enumerate{
#'   \item identify valid treated and control observations around treatment events,
#'   \item match treated and control units on covariates measured at a chosen relative time,
#'   \item expand the matched sample into event-time stacks for estimation, and
#'   \item inspect diagnostics before fitting an event-study model.
#' }
#'
#' \code{flexmatch()} is the main high-level wrapper. It builds clean treatment stacks,
#' matches treated units to valid controls, expands the matched sample into event-time
#' panels, and returns diagnostics for review. The package also ships with
#' \code{flexmatch_example}, a synthetic unbalanced panel dataset that lets users try different
#' matching and expansion options immediately.
#'
#' @param dt A \code{data.table} (or data.frame) containing the panel data.
#' Each row should be one unit-period observation.
#' @param id_col Character. The name of the unique identifier column.
#' This identifies the panel unit, such as a firm, county, school, or person.
#' @param time_col Character. The name of the time identifier column.
#' This can be a year, quarter, date, or any sortable period variable.
#' @param event_col Character. The name of the 0/1 dummy indicating the start of treatment.
#' The variable should equal 1 only in the period when a treatment event begins.
#' @param k_pre Integer. Required number of pre-treatment periods.
#' These periods are retained in the standard event window and are also used when checking
#' that treated and control candidates have enough clean history.
#' @param k_post Integer. Required number of post-treatment periods.
#' These periods are retained in the standard event window and are also used when checking
#' whether recent or upcoming events contaminate a control window.
#' @param k_freeze Integer. Additional cool-down periods used to keep treated units out of the control pool.
#' Increasing this value makes the control pool stricter because recently treated units stay
#' "dirty" for longer after the event window ends.
#' @param match_vars Character vector of covariates used for matching.
#' These variables are measured at the period defined by \code{match_by} and are passed into
#' \code{MatchIt::matchit()}.
#' @param match_exact Optional character vector of exact-match covariates.
#' Use this when matched units must agree exactly on selected dimensions such as industry or region.
#' @param match_by Integer. Relative period used for matching covariates. Default \code{-1}.
#' For example, \code{-1} matches on the period immediately before the event and \code{-2}
#' matches on two periods before the event.
#' @param include_all_years Logical. If \code{FALSE} (default), stacks keep only the standard
#' event window from \code{T-k_pre} through \code{T+k_post}. If \code{TRUE}, each matched
#' unit keeps all periods available in the input data; periods before the standard window are
#' collapsed to \code{-k_pre-1} and periods after the window are collapsed to \code{k_post+1}.
#' @param treated_type Character controlling which treatment events are allowed to define stacks.
#' \itemize{
#'   \item \code{"first"} keeps only the first treatment event observed for each unit.
#'   \item \code{"all"} allows later treatment events too, provided they satisfy the clean-window rules.
#' }
#' @param control_type Character controlling which observations can act as controls.
#' \itemize{
#'   \item \code{"all"} uses any observation with a clean event window.
#'   \item \code{"nyt"} uses only eventually-treated units observed before their first treatment.
#'   \item \code{"nt"} uses observations that are not currently contaminated, including eventually-treated units.
#'   \item \code{"nevert"} restricts controls to units that are never treated in the full sample.
#' }
#' @param assess_quality Logical. If \code{TRUE}, generate covariate balance diagnostics.
#' This returns balance tables and numeric/categorical balance plots when possible.
#' @param plot_events Logical. If \code{TRUE}, generate the event dynamics plot.
#' This is useful for visually checking event timing and matched-sample coverage.
#' @param verbose Logical. If \code{TRUE}, print progress updates.
#' @param ... Additional arguments passed to \code{MatchIt::matchit()}.
#' Typical options include \code{method}, \code{ratio}, \code{replace}, and \code{caliper}.
#'
#' @return A list containing the stacked analysis data, matched cross section,
#' candidate pool, diagnostics, optional plots, and the original \code{matchit} object.
#' Important components are:
#' \itemize{
#'   \item \code{final_expanded_data}: the analysis-ready stacked panel.
#'   \item \code{matched_cross_section}: matched treated/control observations at the matching period.
#'   \item \code{stack_candidates}: all observations that passed the temporal eligibility rules.
#'   \item \code{diagnostics}: counts of candidates, matches, unmatched treated units, and rows dropped for missing covariates.
#'   \item \code{quality_stats}: balance outputs if \code{assess_quality = TRUE}.
#'   \item \code{event_plot}: the event-dynamics plot if \code{plot_events = TRUE}.
#' }
#'
#' @examples
#' data("flexmatch_example")
#' demo_panel <- flexmatch_example[flexmatch_example$unit_id <= 100, ]
#' results <- flexmatch(demo_panel, "unit_id", "year", "event",
#'   k_pre = 1, k_post = 1, k_freeze = 1, match_vars = "size",
#'   method = "nearest", assess_quality = TRUE, plot_events = TRUE,
#'   verbose = FALSE)
#' head(results$final_expanded_data)
#' results$quality_stats$balance_table
#' print(results$event_plot)
#'
#' @export
flexmatch <- function(dt,
                      id_col,
                      time_col,
                      event_col,
                      k_pre,
                      k_post,
                      k_freeze,
                      match_vars,
                      match_exact = NULL,
                      match_by = -1,
                      include_all_years = FALSE,
                      treated_type = "first",
                      control_type = "all",
                      assess_quality = TRUE,
                      plot_events = TRUE,
                      verbose = TRUE,
                      ...) {
  dt <- as_flex_dt(dt)
  validate_flexmatch_inputs(
    dt = dt,
    id_col = id_col,
    time_col = time_col,
    event_col = event_col,
    k_pre = k_pre,
    k_post = k_post,
    k_freeze = k_freeze,
    match_vars = match_vars,
    match_exact = match_exact
  )

  if (verbose) {
    message("\n=== Starting Flexmatch Pipeline ===")
    message(sprintf("1. Scanning Candidates (Treated = '%s', Control = '%s')...", treated_type, control_type))
  }

  cands <- label_candidates(
    dt = dt,
    id_col = id_col,
    time_col = time_col,
    event_col = event_col,
    k_pre = k_pre,
    k_post = k_post,
    k_freeze = k_freeze,
    match_by = match_by,
    treated_type = treated_type,
    control_type = control_type
  )

  n_treat <- nrow(cands[treat == 1])
  if (n_treat == 0) {
    stop("No valid treated candidates found.", call. = FALSE)
  }
  if (verbose) {
    message(sprintf("   Candidates found: %d treated, %d control.", n_treat, nrow(cands[treat == 0])))
    message(sprintf("2. Matching on [%s]...", paste(match_vars, collapse = ", ")))
  }

  match_res <- match_candidates(
    dt = dt,
    stack_candidates = cands,
    id_col = id_col,
    time_col = time_col,
    match_vars = match_vars,
    match_exact = match_exact,
    ...
  )

  diagnostics <- format_match_diagnostics(
    stack_candidates = cands,
    matched_cross_section = match_res$matched_cross_section,
    id_col = id_col,
    dropped_rows_missing_covariates = match_res$diagnostics$dropped_rows_missing_covariates
  )

  n_matched <- diagnostics$matching$treated_matched
  if (verbose) {
    message(sprintf("   Matched %d / %d treated units (%.1f%%).", n_matched, diagnostics$matching$treated_candidates, 100 * n_matched / diagnostics$matching$treated_candidates))
    message("3. Expanding Stacks...")
  }

  final_data <- expand_matched_stacks(
    dt = dt,
    matched_data = match_res$matched_cross_section,
    id_col = id_col,
    time_col = time_col,
    k_pre = k_pre,
    k_post = k_post,
    match_by = match_by,
    include_all_years = include_all_years
  )

  quality_stats <- NULL
  if (assess_quality) {
    if (verbose) {
      message("4. Assessing Quality...")
    }
    tryCatch({
      quality_stats <- assess_match_quality(match_res, cands, dt, match_vars, id_col, time_col)
    }, error = function(e) {
      warning(paste("Quality Check Failed:", e$message), call. = FALSE)
    })
  }

  event_plot <- NULL
  if (plot_events) {
    if (verbose) {
      message("5. Generating Event Dynamics Plot...")
    }
    tryCatch({
      ep <- plot_event_dynamics(
        dt = dt,
        dt_matched = final_data,
        id_col = id_col,
        time_col = time_col,
        event_col = event_col,
        sort_by_intensity = FALSE
      )
      event_plot <- ep$plot
    }, error = function(e) {
      warning(paste("Event Plot Failed:", e$message), call. = FALSE)
    })
  }

  if (verbose) {
    message("=== Pipeline Complete ===\n")
  }

  list(
    final_expanded_data = final_data,
    matched_cross_section = match_res$matched_cross_section,
    stack_candidates = cands,
    diagnostics = diagnostics,
    quality_stats = quality_stats,
    event_plot = event_plot,
    match_object = match_res$match_object
  )
}
