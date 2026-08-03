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
#' @param projection_search Strategy used to locate the minimum wait time.
#'   `"binary"` (default) screens the top of the candidate range and then
#'   bisects, costing `O(log t_max)` model fits. `"linear"` scans every
#'   candidate day, costing `O(t_max)` fits. Both return the same answer
#'   whenever the projected dose is monotone non-decreasing in wait time,
#'   which holds for TITE-CRM and TITE-BOIN; `"linear"` is retained as a
#'   verification oracle and an escape hatch for non-monotone parents.
#' @param verbose If `TRUE` (default), emit an informational message when
#'   `t_max` covers all remaining follow-up. Set to `FALSE` in simulation,
#'   where this fires on nearly every projection. Affects messaging only,
#'   never the returned value.
#' @return A `dtp_selector_factory` object.
#' @export
apply_dtp <- function(parent_selector_factory, t_max, obswin,
                      projection_search = c("binary", "linear"),
                      verbose = TRUE) {
  stopifnot(is.numeric(t_max), length(t_max) == 1, t_max >= 0)
  stopifnot(is.numeric(obswin), length(obswin) == 1, obswin > 0)
  projection_search <- match.arg(projection_search)
  stopifnot(is.logical(verbose), length(verbose) == 1, !is.na(verbose))

  # Mutable settings (e.g. queue_size) stored in an environment so that

  # set_dtp_queue_size() can modify them through R's reference semantics.
  sim_settings <- new.env(parent = emptyenv())
  sim_settings$queue_size <- NULL

  x <- list(
    parent = parent_selector_factory,
    t_max = t_max,
    obswin = obswin,
    projection_search = projection_search,
    verbose = verbose,
    sim_settings = sim_settings
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

  # Nothing to project toward once the recommendation is already the top dose:
  # `.dtp_project()` only reports a wait when the projected dose rises strictly
  # above `dose_recommended`, which no amount of follow-up can do here. Given
  # the `dose_recommended <= dose_current` guard below, this is exactly the
  # "at the top dose and staying" case; a de-escalation off the top dose still
  # projects, since recovering to the top dose is a genuine hit.
  if (!is.na(dose_recommended) && !is.na(dose_current) &&
      dose_recommended <= dose_current && selector_factory$t_max > 0 &&
      dose_recommended < num_doses(parent_selector)) {
    projection <- .dtp_project(
      parent_factory = selector_factory$parent,
      outcomes = outcomes,
      dose_recommended = dose_recommended,
      t_max = selector_factory$t_max,
      obswin = selector_factory$obswin,
      # Defaults keep factories built before these fields existed working.
      search = selector_factory$projection_search %||% "binary",
      # sim_settings wins when set: the simulation driver silences the
      # message there (reference semantics) without rebuilding the chain.
      verbose = selector_factory$sim_settings$verbose %||%
        selector_factory$verbose %||% TRUE,
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
                         t_max, obswin, search = "binary", verbose = TRUE,
                         ...) {
  no_wait <- list(wait_time = 0, weight_gain = 0, projected_dose = NA_integer_)

  is_pending <- outcomes$tox == 0 & outcomes$weight < 1

  if (!any(is_pending)) {
    return(no_wait)
  }

  max_remaining_weight <- max(1 - outcomes$weight[is_pending])
  effective_max_time <- max_remaining_weight * obswin
  # Floored to an integer number of days: a fractional t_max would otherwise
  # let the search probe a step that `seq_len()` (which truncates) never
  # reached, so the two strategies would disagree.
  t_max_effective <- floor(min(t_max, ceiling(effective_max_time)))

  if (verbose && is.finite(t_max) && t_max / obswin >= max_remaining_weight) {
    message(
      "t_max (", t_max, ") covers all remaining follow-up (effective max: ",
      round(effective_max_time, 1), "); equivalent to unlimited waiting."
    )
  }

  if (t_max_effective < 1) {
    return(no_wait)
  }

  weight_per_unit <- 1 / obswin
  pending_weight <- outcomes$weight[is_pending]

  # Recommended dose after projecting `step` further days of DLT-free
  # follow-up onto every pending patient. NA propagates (a stopping rule in
  # the parent chain may decline to recommend anything).
  rec_at <- function(step) {
    projected_outcomes <- outcomes
    projected_outcomes$weight[is_pending] <- pmin(
      1, pending_weight + step * weight_per_unit
    )
    rec <- recommended_dose(parent_factory |> fit(projected_outcomes, ...))
    if (is.na(rec)) NA_integer_ else as.integer(rec)
  }

  hits <- function(rec) !is.na(rec) && rec > dose_recommended

  # Exhaustive scan: the reference behaviour, and the fallback below.
  linear_scan <- function() {
    for (step in seq_len(t_max_effective)) {
      rec <- rec_at(step)
      if (hits(rec)) {
        return(list(
          wait_time = step,
          weight_gain = step * weight_per_unit,
          projected_dose = rec
        ))
      }
    }
    no_wait
  }

  if (!identical(search, "binary")) {
    return(linear_scan())
  }

  # The projected dose is monotone non-decreasing in wait time for TITE-CRM
  # and TITE-BOIN (extra DLT-free follow-up can only push the toxicity
  # estimate down)
  top <- rec_at(t_max_effective)

  if (is.na(top)) {
    return(linear_scan())
  }

  if (!hits(top)) {
    return(no_wait)
  }

  # Invariant: !hits(lo), hits(hi). Bisect down to the first hitting step,
  # carrying hi's dose so no extra fit is needed to report it.
  lo <- 0L
  hi <- as.integer(t_max_effective)
  rec_hi <- top
  while (hi - lo > 1L) {
    mid <- (lo + hi) %/% 2L
    rec_mid <- rec_at(mid)
    if (hits(rec_mid)) {
      hi <- mid
      rec_hi <- rec_mid
    } else {
      lo <- mid
    }
  }

  list(
    wait_time = hi,
    weight_gain = hi * weight_per_unit,
    projected_dose = rec_hi
  )
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
