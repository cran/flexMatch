test_that("asymmetric and zero windows agree with an independent window scan", {
  d <- data.table::CJ(unit = 1:4, period = 1:12)
  d[, event := 0L]
  d[unit == 1 & period == 6, event := 1L]
  d[unit == 2 & period %in% c(4, 9), event := 1L]
  d[unit == 3 & period == 10, event := 1L]
  for (window in list(c(0,0),c(1,3),c(3,1))) {
    pre <- window[1]; post <- window[2]; freeze <- 1
    got <- label_candidates(d,"unit","period","event",pre,post,freeze,
      match_by=0,treated_type="all")
    expected <- character()
    for (u in unique(d$unit)) {
      events <- d[unit==u & event==1,period]
      dirty_past <- vapply(1:12,function(t) any(events<=t & t-events<=post+freeze),logical(1))
      dirty_future <- vapply(1:12,function(t) any(events>=t & events-t<=post),logical(1))
      for (t in 1:12) {
        if (t-pre<1 || t+post>12) next
        pre_rows <- if (pre>0) seq.int(t-pre,t-1) else integer()
        clean_pre <- !any(dirty_past[pre_rows])
        clean_future_events <- !any(events>t & events<=t+post)
        if (t %in% events && clean_pre && clean_future_events) expected <- c(expected,paste(u,t,1))
        if (!any((dirty_past|dirty_future)[seq.int(t-pre,t+post)])) expected <- c(expected,paste(u,t,0))
      }
    }
    expect_setequal(paste(got$unit,got$cohort_year,got$treat),expected)
  }
})

test_that("reused controls have distinct unit-within-set identifiers", {
  d <- data.table::CJ(unit=1:6,year=1:7)
  d[, event:=as.integer(unit<=3 & year==4)]
  d[, size:=as.numeric(unit)]
  c <- label_candidates(d,"unit","year","event",1,1,0)
  m <- match_candidates(d,c,"unit","year","size",method="nearest",replace=TRUE,distance="mahalanobis")
  z <- m$matched_cross_section
  expect_equal(nrow(z),6L)
  expect_equal(length(unique(z$cohort_id)),3L)
  expect_equal(anyDuplicated(z$stack_id),0L)
  expanded <- expand_matched_stacks(d,z,"unit","year",1,1)
  expect_equal(anyDuplicated(expanded[,.(stack_id,year)]),0L)
  expect_identical(attr(expanded,"flexmatch_id_col"),"unit")
})

test_that("invalid identifiers and missing treatment are rejected", {
  d <- data.table::CJ(unit=1:3,year=1:5)
  d[,event:=0L]
  d[1,event:=NA_integer_]
  expect_error(label_candidates(d,"unit","year","event",1,1,0),"missing")
  d[1,event:=0L]
  expect_error(label_candidates(d,"unit","year","event",1,1,0,match_by=1),"non-positive")
})

test_that("exact matching retains one set without a cartesian treated-map merge", {
  d <- data.table::CJ(unit=1:6,year=1:7)
  d[, event:=as.integer(unit<=3 & year==4)]
  d[, size:=1]
  c <- label_candidates(d,"unit","year","event",1,1,0)
  m <- match_candidates(d,c,"unit","year","size",method="exact")
  z <- m$matched_cross_section
  expect_equal(nrow(z),6L)
  expect_equal(length(unique(z$cohort_id)),1L)
  expect_equal(anyDuplicated(z$stack_id),0L)
  expect_equal(nrow(expand_matched_stacks(d,z,"unit","year",1,1)),18L)
})


test_that("the main wrapper returns drawable timing and balance plots", {
  data("flexmatch_example")
  d <- flexmatch_example[flexmatch_example$unit_id <= 100, ]
  r <- flexmatch(d,"unit_id","year","event",1,1,1,"size",
    method="nearest",verbose=FALSE)
  expect_s3_class(r$event_plot,"ggplot")
  expect_s3_class(r$quality_stats$numeric_balance_panel,"ggplot")
})
