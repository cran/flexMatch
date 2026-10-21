#' Example Panel Dataset for `flexmatch()`
#'
#' A large unbalanced unit-year panel dataset shipped with the package to demonstrate
#' treatment timing, matching, and stack expansion options.
#'
#' The dataset was designed to let users experiment with:
#' \itemize{
#'   \item \code{treated_type = "first"} versus \code{"all"},
#'   \item \code{control_type = "all"}, \code{"nyt"}, \code{"nt"}, and \code{"nevert"},
#'   \item exact matching with \code{match_exact},
#'   \item standard event windows versus \code{include_all_years = TRUE}, and
#'   \item downstream estimation with \code{fit_event_study()}.
#' }
#'
#' @format A data frame with 65000 rows and 8 variables:
#' \describe{
#'   \item{unit_id}{Integer unit identifier.}
#'   \item{year}{Calendar year.}
#'   \item{industry}{Industry category used for exact matching examples.}
#'   \item{region}{Geographic region label.}
#'   \item{size}{Continuous covariate for matching.}
#'   \item{leverage}{Continuous covariate for matching.}
#'   \item{event}{Binary indicator equal to 1 in periods when treatment starts.}
#'   \item{outcome}{Example outcome variable with treatment-driven dynamics.}
#' }
#'
#' @details
#' The dataset has 10,000 units observed between 2013 and 2022, but units enter
#' and exit the panel at different dates, so histories range from 3 to 10
#' years. Units are assigned deterministically to treatment patterns including
#' never treated, once treated, late treated, sparse repeated treatment, and
#' frequent treatment. This makes the dataset useful for illustrating the practical
#' differences between \code{treated_type}, \code{control_type},
#' \code{match_exact}, and \code{include_all_years}.
#'
#' @examples
#' data("flexmatch_example")
#' head(flexmatch_example)
#' with(flexmatch_example, table(event))
#'
"flexmatch_example"
