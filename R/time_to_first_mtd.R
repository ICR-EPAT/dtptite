# Time-to-first-MTD simulation metric
#
# Three user-facing pieces:
#
#   scenario_mtd(true_prob_tox, target)  — the scenario's true MTD
#   time_to_first_mtd(x)                 — per-replicate event times
#   summary() / plot() methods           — across-replicate roll-up and curve
#
# Unlike the dtp_wait_* metrics these apply to any design — DTP or not, TITE
# or not — because they read the patient-level dose/time frame that every
# escalation selector carries, rather than a DTP-specific event log. The frame
# is cumulative, so only the last fit of each replicate is inspected and
# `return_all_fits = TRUE` is not required.

# -- Internal: locate the patient-level dose/time frame ---------------------
#
# escalation's `model_frame()` drops the `time` column, so the frame has to be
# found by walking the decorator chain. Follows the recursion in
# `.find_dtp_factory()`.

.find_dose_time_df <- function(fit) {
  if (is.null(fit)) return(NULL)
  df <- fit$df
  if (!is.null(df) && all(c("dose", "time") %in% names(df))) {
    return(as.data.frame(df[, c("dose", "time")]))
  }
  .find_dose_time_df(fit$parent)
}

.patient_dose_times <- function(fit) {
  df <- .find_dose_time_df(fit)
  if (is.null(df)) {
    stop("Could not find a patient-level data frame carrying both `dose` and ",
         "`time` in this selector chain. `time_to_first_mtd()` needs the ",
         "simulation to have recorded patient times.", call. = FALSE)
  }
  df
}

# Last fit entry of one replicate (same idiom as dtp_wait_events).
.last_fit <- function(replicate_fits) {
  replicate_fits[[length(replicate_fits)]]$fit
}

# -- scenario_mtd ------------------------------------------------------

#' True MTD of a simulated scenario.
#'
#' Identifies the dose whose true toxicity probability lies within `tol` of the
#' target. When several doses qualify the lowest is returned. When none
#' qualifies the scenario has no MTD and `NA` is returned with a warning.
#'
#' @param true_prob_tox Numeric vector of true toxicity probabilities, one per
#'   dose.
#' @param target Target toxicity probability.
#' @param tol Half-width of the acceptance band. A dose qualifies when
#'   `abs(true_prob_tox - target) <= tol`. Default `0.05`.
#' @return Integer dose index, or `NA_integer_` when no dose qualifies.
#' @export
scenario_mtd <- function(true_prob_tox, target, tol = 0.05) {
  if (!is.numeric(true_prob_tox) || !length(true_prob_tox)) {
    stop("`true_prob_tox` must be a non-empty numeric vector.", call. = FALSE)
  }
  if (!is.numeric(target) || length(target) != 1L || is.na(target)) {
    stop("`target` must be a single non-missing number.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol < 0) {
    stop("`tol` must be a single non-negative number.", call. = FALSE)
  }

  band <- which(abs(true_prob_tox - target) <= tol)
  if (!length(band)) {
    if (all(true_prob_tox > target + tol)) {
      warning("No MTD: every dose exceeds target + tol (", target + tol, ").",
              call. = FALSE)
    } else if (all(true_prob_tox < target - tol)) {
      warning("No MTD: every dose falls below target - tol (", target - tol,
              ").", call. = FALSE)
    } else {
      warning("No MTD: no dose lies within tol (", tol, ") of the target (",
              target, "). Consider a wider `tol` or supplying `mtd` directly.",
              call. = FALSE)
    }
    return(NA_integer_)
  }
  as.integer(min(band))
}

# -- Internal: resolve the MTD for a simulations object --------------------

.resolve_mtd <- function(x, mtd, tol) {
  num_doses <- length(x$true_prob_tox)

  if (!is.null(mtd)) {
    if (length(mtd) != 1L) {
      stop("`mtd` must be a single value (or NULL to derive it).", call. = FALSE)
    }
    if (is.na(mtd)) return(NA_integer_)
    mtd <- as.integer(mtd)
    if (mtd < 1L || mtd > num_doses) {
      stop("`mtd` must lie between 1 and ", num_doses, ".", call. = FALSE)
    }
    return(mtd)
  }

  if (is.null(x$true_prob_tox)) {
    stop("This object carries no `true_prob_tox`; supply `mtd` explicitly.",
         call. = FALSE)
  }
  target <- tryCatch(tox_target(.last_fit(x$fits[[1]])),
                     error = function(e) NULL)
  if (is.null(target) || !is.numeric(target) || is.na(target)) {
    stop("Could not read the toxicity target from the fitted design; supply ",
         "`mtd` explicitly.", call. = FALSE)
  }
  scenario_mtd(x$true_prob_tox, target, tol = tol)
}

# -- Internal: first-MTD event for one replicate ---------------------------
#
# `n_before_mtd` counts patients dosed strictly earlier than the first MTD
# patient, so it does not depend on how ties are ordered — DTP doses a queued
# cohort at a single shared timestamp.

.first_mtd_event <- function(patients, mtd) {
  if (is.na(mtd)) {
    return(list(time = NA_real_, reached = NA, n_before = NA_integer_))
  }
  hits <- patients$time[patients$dose == mtd]
  if (!length(hits)) {
    return(list(time = NA_real_, reached = FALSE, n_before = NA_integer_))
  }
  first <- min(hits)
  list(time = first,
       reached = TRUE,
       n_before = as.integer(sum(patients$time < first)))
}

# -- time_to_first_mtd -------------------------------------------------

#' Time until the first patient is dosed at the true MTD.
#'
#' Returns one row per simulated replicate recording when, if ever, a patient
#' was first treated at the scenario's true MTD.
#'
#' @details
#' Times are dosing times measured from the start of the trial (`t = 0`). For
#' designs without DTP a patient is dosed on arrival, so dosing time equals
#' recruitment time; for patients held in a DTP queue it is the instant the
#' wait ends and they are dosed.
#'
#' Columns:
#'
#' - `replicate` — index of the simulated trial (and `design` for a
#'   `simulations_collection`)
#' - `mtd` — the dose treated as the true MTD
#' - `time_to_mtd` — dosing time of the first patient given `mtd`, or `NA` if
#'   no patient was given it before the trial ended
#' - `reached` — whether any patient was dosed at `mtd`; `NA` when the scenario
#'   has no MTD
#' - `n_before_mtd` — patients dosed strictly earlier than that patient, or
#'   `NA` if the MTD was never given
#'
#' When the scenario has no MTD (no dose within `tol` of the target) every row
#' carries `mtd = NA` and `reached = NA`, distinguishing "no MTD to reach" from
#' "an MTD existed and was not reached" (`reached = FALSE`).
#'
#' The number of replicates is recorded in a `num_sims` attribute and is used
#' as the denominator by [summary.time_to_first_mtd()].
#'
#' @param x A `simulations` or `simulations_collection` object.
#' @param mtd Dose index to treat as the true MTD. `NULL` (default) derives it
#'   with [scenario_mtd()]; `NA` declares the scenario to have no MTD.
#' @param tol Passed to [scenario_mtd()] when deriving the MTD. Default `0.05`.
#' @param ... Unused.
#' @return A tibble of class `time_to_first_mtd`. See Details for columns.
#' @export
time_to_first_mtd <- function(x, mtd = NULL, tol = 0.05, ...) {
  UseMethod("time_to_first_mtd")
}

#' @export
time_to_first_mtd.simulations <- function(x, mtd = NULL, tol = 0.05, ...) {
  num_sims <- length(x$fits)
  mtd <- .resolve_mtd(x, mtd, tol)

  events <- lapply(seq_len(num_sims), function(i) {
    .first_mtd_event(.patient_dose_times(.last_fit(x$fits[[i]])), mtd)
  })

  out <- tibble::tibble(
    replicate    = seq_len(num_sims),
    mtd          = rep(mtd, num_sims),
    time_to_mtd  = vapply(events, function(e) e$time, numeric(1)),
    reached      = vapply(events, function(e) e$reached, logical(1)),
    n_before_mtd = vapply(events, function(e) e$n_before, integer(1))
  )
  .as_time_to_first_mtd(out, num_sims)
}

#' @export
time_to_first_mtd.simulations_collection <- function(x, mtd = NULL,
                                                     tol = 0.05, ...) {
  parts <- lapply(x, time_to_first_mtd, mtd = mtd, tol = tol)
  num_sims <- vapply(parts, function(p) attr(p, "num_sims"), integer(1))
  out <- dplyr::bind_rows(lapply(parts, tibble::as_tibble), .id = "design")
  .as_time_to_first_mtd(out, num_sims)
}

.as_time_to_first_mtd <- function(x, num_sims) {
  x <- tibble::as_tibble(x)
  attr(x, "num_sims") <- num_sims
  class(x) <- c("time_to_first_mtd", class(x))
  x
}

# -- summary -----------------------------------------------------------

# Quantile of the cumulative incidence of reaching the MTD. Replicates that
# never reached are coded +Inf so that stats::quantile() returns Inf — and
# hence NA — for any probability beyond the reached fraction.
.cif_quantile <- function(times, num_sims, prob) {
  v <- ifelse(is.na(times), Inf, times)
  if (length(v) < num_sims) v <- c(v, rep(Inf, num_sims - length(v)))
  q <- stats::quantile(v, probs = prob, names = FALSE)
  if (is.infinite(q)) NA_real_ else q
}

.summarise_one <- function(d, num_sims, probs) {
  if (nrow(d) != num_sims) {
    warning("Object has ", nrow(d), " rows but the simulation ran ", num_sims,
            " replicates; it appears to have been subset. Absent replicates ",
            "are counted as not reaching the MTD.", call. = FALSE)
  }
  mtd <- d$mtd[1]
  num_reached <- sum(!is.na(d$time_to_mtd))

  out <- tibble::tibble(
    num_sims            = as.integer(num_sims),
    mtd                 = as.integer(mtd),
    num_reached         = as.integer(num_reached),
    reached_fraction    = if (is.na(mtd)) NA_real_ else num_reached / num_sims,
    median_time_to_mtd  = .cif_quantile(d$time_to_mtd, num_sims, 0.5),
    median_n_before_mtd = if (num_reached > 0L) {
      stats::median(d$n_before_mtd, na.rm = TRUE)
    } else {
      NA_real_
    }
  )

  for (p in probs) {
    out[[paste0("q", round(p * 100), "_time_to_mtd")]] <-
      .cif_quantile(d$time_to_mtd, num_sims, p)
  }
  out
}

#' Summarise time-to-first-MTD across replicates.
#'
#' @details
#' Columns:
#'
#' - `num_sims` — replicates simulated
#' - `mtd` — dose treated as the true MTD
#' - `num_reached` — replicates in which a patient was dosed at `mtd`
#' - `reached_fraction` — `num_reached / num_sims`; `NA` when the scenario has
#'   no MTD
#' - `median_time_to_mtd` — median of the cumulative incidence of reaching the
#'   MTD. `NA` unless more than half of replicates reached it, because the
#'   quantile is otherwise undefined — including when exactly half did
#' - `median_n_before_mtd` — median patients dosed before the MTD, over
#'   replicates that reached it
#'
#' Replicates that never reached the MTD remain in the denominator; they are
#' not dropped. Any probability supplied in `probs` adds a
#' `q<pct>_time_to_mtd` column, computed the same way and likewise `NA` beyond
#' `reached_fraction`.
#'
#' @param object A `time_to_first_mtd` object.
#' @param probs Optional numeric vector of probabilities for additional
#'   quantile columns. Default `NULL` adds none.
#' @param ... Unused.
#' @return A tibble with one row per design.
#' @export
summary.time_to_first_mtd <- function(object, probs = NULL, ...) {
  if (!is.null(probs)) {
    if (!is.numeric(probs) || any(is.na(probs)) ||
        any(probs < 0) || any(probs > 1)) {
      stop("`probs` must be numbers between 0 and 1.", call. = FALSE)
    }
  }
  num_sims <- attr(object, "num_sims")
  if (is.null(num_sims)) num_sims <- nrow(object)

  if (!"design" %in% names(object)) {
    return(.summarise_one(object, unname(num_sims)[1], probs))
  }

  designs <- unique(object$design)
  rows <- lapply(designs, function(d) {
    sub <- object[object$design == d, , drop = FALSE]
    n <- if (!is.null(names(num_sims)) && d %in% names(num_sims)) {
      num_sims[[d]]
    } else {
      nrow(sub)
    }
    cbind(design = d, .summarise_one(sub, n, probs))
  })
  dplyr::bind_rows(rows)
}

# -- plot --------------------------------------------------------------

# Wilson score interval for a binomial proportion — behaves near 0 and 1,
# where the cumulative incidence curve spends most of its length.
.wilson_interval <- function(k, n, conf = 0.95) {
  z <- stats::qnorm(1 - (1 - conf) / 2)
  phat <- k / n
  denom <- 1 + z^2 / n
  centre <- (phat + z^2 / (2 * n)) / denom
  half <- z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2)) / denom
  data.frame(lower = pmax(0, centre - half), upper = pmin(1, centre + half))
}

.cif_steps <- function(values, num_sims, conf) {
  finite <- sort(values[!is.na(values)])
  steps <- data.frame(
    x = c(0, finite),
    y = c(0, seq_along(finite)) / num_sims
  )
  ci <- .wilson_interval(c(0, seq_along(finite)), num_sims, conf)
  cbind(steps, ci)
}

#' Plot the cumulative incidence of reaching the MTD.
#'
#' Draws, for each design, the proportion of replicates in which a patient had
#' been dosed at the true MTD by a given point. The curve plateaus at
#' `reached_fraction`; replicates that never reached the MTD keep it below 1.
#'
#' @param x A `time_to_first_mtd` object.
#' @param axis `"time"` (default) plots calendar time since trial start;
#'   `"patients"` plots the number of patients dosed beforehand.
#' @param ci If `TRUE`, add pointwise Wilson intervals for Monte Carlo error.
#'   Default `FALSE`.
#' @param conf Confidence level for the intervals. Default `0.95`.
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
plot.time_to_first_mtd <- function(x, axis = c("time", "patients"),
                                   ci = FALSE, conf = 0.95, ...) {
  axis <- match.arg(axis)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed to plot; use summary() instead.",
         call. = FALSE)
  }

  num_sims <- attr(x, "num_sims")
  if (is.null(num_sims)) num_sims <- nrow(x)
  values <- if (axis == "time") x$time_to_mtd else as.numeric(x$n_before_mtd)
  # n_before_mtd is NA exactly when the MTD was not reached, so both axes
  # carry the same set of events.

  if (!"design" %in% names(x)) {
    curves <- cbind(design = "all",
                    .cif_steps(values, unname(num_sims)[1], conf))
  } else {
    curves <- dplyr::bind_rows(lapply(unique(x$design), function(d) {
      keep <- x$design == d
      n <- if (!is.null(names(num_sims)) && d %in% names(num_sims)) {
        num_sims[[d]]
      } else {
        sum(keep)
      }
      cbind(design = d, .cif_steps(values[keep], n, conf))
    }))
  }

  xlab <- if (axis == "time") "Time since trial start" else "Patients dosed"
  p <- ggplot2::ggplot(curves, ggplot2::aes(x = x, y = y, colour = design)) +
    ggplot2::geom_step() +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = xlab, y = "Proportion reaching MTD", colour = "Design")

  if (ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper, fill = design),
      alpha = 0.15, colour = NA, outline.type = "full"
    ) + ggplot2::labs(fill = "Design")
  }
  p
}

# Silence R CMD check NOTE: undefined globals used via NSE in ggplot2 aes
utils::globalVariables(c("x", "y", "design", "lower", "upper"))
