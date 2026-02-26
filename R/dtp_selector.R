# DTP (Dose Transition Pathways) decorator for TITE dose-finding models
#
# Wraps any TITE selector factory with a look-ahead projection that
# determines whether waiting for pending patients' follow-up could
# change the dose recommendation.

# -- Generics ----------------------------------------------------------

#' Should the trialist wait before dosing the next cohort?
#'
#' Returns `TRUE` if the DTP projection found that waiting could lead to
#' a higher dose recommendation within `t_max`.
#'
#' @param x A `dtp_selector` object.
#' @param ... Additional arguments (unused).
#' @return Logical scalar.
#' @export
dtp_should_wait <- function(x, ...) UseMethod("dtp_should_wait")

#' Minimum wait time before the dose recommendation changes.
#'
#' @param x A `dtp_selector` object.
#' @param type Character: `"time"` (default) returns wait in the user's
#'   time units (same as `t_max` and `obswin`); `"weight"` returns the
#'   corresponding weight increment.
#' @param ... Additional arguments (unused).
#' @return Numeric scalar (0 if no benefit from waiting).
#' @export
dtp_wait_time <- function(x, type = c("time", "weight"), ...) {
  UseMethod("dtp_wait_time")
}

#' Projected dose after waiting.
#'
#' @param x A `dtp_selector` object.
#' @param ... Additional arguments (unused).
#' @return Integer dose level, or `NA` if waiting has no benefit.
#' @export
dtp_projected_dose <- function(x, ...) UseMethod("dtp_projected_dose")

# -- Factory -----------------------------------------------------------

#' Apply the DTP look-ahead to a TITE selector factory.
#'
#' Wraps a parent TITE selector factory with the Dose Transition Pathways
#' projection engine. At fit time, the decorator checks whether waiting
#' for pending patients' follow-up (up to `t_max`) could change the dose
#' recommendation upward.
#'
#' @param parent_selector_factory A `selector_factory` (typically from
#'   [escalation::get_dfcrm_tite()] or similar TITE model).
#' @param t_max Maximum allowable wait time (numeric >= 0, or `Inf` for
#'   unlimited). Must use the same time unit as `obswin`.
#' @param obswin DLT observation window length (e.g. 56 days or 8 weeks).
#'   Defines the time-weight mapping: `weight = time_elapsed / obswin`.
#' @return A `dtp_selector_factory` object.
#' @export
apply_dtp <- function(parent_selector_factory, t_max, obswin) {
  stopifnot(is.numeric(t_max), length(t_max) == 1, t_max >= 0)
  stopifnot(is.numeric(obswin), length(obswin) == 1, obswin > 0)

  x <- list(
    parent = parent_selector_factory,
    t_max = t_max,
    obswin = obswin
  )
  class(x) <- c("dtp_selector_factory",
                 "derived_dose_selector_factory",
                 "selector_factory")
  x
}

# -- Fit ---------------------------------------------------------------

#' @export
fit.dtp_selector_factory <- function(selector_factory, outcomes, ...) {
  parent_selector <- selector_factory$parent |>
    fit(outcomes, ...)

  dose_recommended <- recommended_dose(parent_selector)

  all_doses <- doses_given(parent_selector)
  dose_current <- if (length(all_doses) > 0) all_doses[length(all_doses)] else NA

  projection <- list(
    wait_time = 0,
    weight_gain = 0,
    projected_dose = NA_integer_
  )

  if (!is.na(dose_recommended) && !is.na(dose_current) &&
      dose_recommended <= dose_current && selector_factory$t_max > 0) {
    projection <- .dtp_project(
      parent_factory = selector_factory$parent,
      outcomes = outcomes,
      dose_recommended = dose_recommended,
      t_max = selector_factory$t_max,
      obswin = selector_factory$obswin,
      ...
    )
  }

  dtp_selector(
    parent_selector = parent_selector,
    t_max = selector_factory$t_max,
    obswin = selector_factory$obswin,
    wait_time = projection$wait_time,
    weight_gain = projection$weight_gain,
    projected_dose = projection$projected_dose
  )
}

# -- Projection engine -------------------------------------------------

.dtp_project <- function(parent_factory, outcomes, dose_recommended,
                         t_max, obswin, ...) {
  is_pending <- outcomes$tox == 0 & outcomes$weight < 1

  if (!any(is_pending)) {
    return(list(wait_time = 0, weight_gain = 0, projected_dose = NA_integer_))
  }

  max_remaining_weight <- max(1 - outcomes$weight[is_pending])
  effective_max_time <- max_remaining_weight * obswin
  t_max_effective <- min(t_max, ceiling(effective_max_time))

  if (is.finite(t_max) && t_max / obswin >= max_remaining_weight) {
    message(
      "t_max (", t_max, ") covers all remaining follow-up (effective max: ",
      round(effective_max_time, 1), "); equivalent to unlimited waiting."
    )
  }

  weight_per_unit <- 1 / obswin

  for (step in seq_len(t_max_effective)) {
    cumulative_gain <- step * weight_per_unit

    projected_outcomes <- outcomes
    projected_outcomes$weight[is_pending] <- pmin(
      1, outcomes$weight[is_pending] + cumulative_gain
    )

    projected_fit <- parent_factory |> fit(projected_outcomes, ...)
    new_recommended_dose <- recommended_dose(projected_fit)

    if (!is.na(new_recommended_dose) &&
        new_recommended_dose > dose_recommended) {
      return(list(
        wait_time = step,
        weight_gain = cumulative_gain,
        projected_dose = as.integer(new_recommended_dose)
      ))
    }
  }

  list(wait_time = 0, weight_gain = 0, projected_dose = NA_integer_)
}

# -- Selector constructor ----------------------------------------------

dtp_selector <- function(parent_selector, t_max, obswin,
                         wait_time, weight_gain, projected_dose) {
  l <- list(
    parent = parent_selector,
    t_max = t_max,
    obswin = obswin,
    wait_time = wait_time,
    weight_gain = weight_gain,
    projected_dose = projected_dose
  )
  class(l) <- c("dtp_selector", "derived_dose_selector", "selector")
  l
}

# -- S3 methods: DTP decision-support ----------------------------------

#' @export
dtp_should_wait.dtp_selector <- function(x, ...) {
  x$wait_time > 0
}

#' @export
dtp_wait_time.dtp_selector <- function(x, type = c("time", "weight"), ...) {
  type <- match.arg(type)
  if (type == "time") x$wait_time else x$weight_gain
}

#' @export
dtp_projected_dose.dtp_selector <- function(x, ...) {
  x$projected_dose
}

# -- Forwarding through outer decorators ----------------------------------
# escalation's derived_dose_selector only forwards known generics.
# These methods ensure DTP generics propagate through any outer decorator.

#' @export
dtp_should_wait.derived_dose_selector <- function(x, ...) {
  dtp_should_wait(x$parent, ...)
}

#' @export
dtp_wait_time.derived_dose_selector <- function(x, type = c("time", "weight"), ...) {
  dtp_wait_time(x$parent, type = type, ...)
}

#' @export
dtp_projected_dose.derived_dose_selector <- function(x, ...) {
  dtp_projected_dose(x$parent, ...)
}
