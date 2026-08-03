# DTP-specific simulation metrics
#
# Two accessors that dispatch on both `simulations` (single design) and
# `simulations_collection` (multiple arms from `simulate_compare`):
#
#   dtp_wait_events(x)           — raw per-wait event log
#   dtp_wait_summary(x, by_dose) — curated per-trial (or per-trial × dose) roll-up
#
# Both are driven by the event log attached to the last fit entry of each
# replicate by `phase1_dtp_tite_sim()`. Non-DTP designs (e.g. a plain
# TITE-CRM baseline through `simulate_compare`) carry no such field, which
# is how we distinguish "did not run under DTP" (all-NA metrics) from
# "ran under DTP and never waited" (num_waits = 0).

# -- dtp_wait_events ---------------------------------------------------

#' Raw DTP wait event log from simulated trials.
#'
#' Returns a tidy tibble with one row per DTP wait event fired across all
#' simulated replicates. This is the single source of truth from which all
#' DTP wait summaries are derived; users who want metrics not covered by
#' [dtp_wait_summary()] can compute them directly from this table.
#'
#' Replicates with zero waits contribute zero rows (but are still reflected
#' in `dtp_wait_summary()` with explicit `0`s). Replicates from non-DTP
#' designs contribute zero rows and are flagged by the design carrying no
#' event log at all, not by reporting zero waits.
#'
#' @param x A `simulations` or `simulations_collection` object.
#' @param ... Unused.
#' @return A tibble with columns `replicate` (and `design` for a
#'   `simulations_collection`), `cohort_idx`, `dose_before`, `time_in`,
#'   `time_out`, `wait_duration`, `projected_dose`, `dose_after`,
#'   `ended_by` (factor: `wait_end`/`dtp_decision`/`stopped`), `dose_delta`,
#'   `effective`, `num_extensions`, `n_held`, `n_missed`.
#'
#'   `wait_duration` is the total elapsed time between when the wait was
#'   triggered and when it ended, including any extensions caused by
#'   in-wait DLT events: each such event prompts a re-fit, and if the
#'   updated fit still requests a wait the original `wait_end` is bumped
#'   out. `num_extensions` records how many times that happened (0 means
#'   the wait ran exactly to its initial projection); as a result
#'   `wait_duration` can exceed the design's `t_max`.
#'
#'   `n_held` is the number of patients held through the wait and dosed at
#'   `wait_end`, **including** the cohort opener (the patient the wait is for);
#'   it is at least 1 for any wait that doses and ranges up to the cohort size.
#'   `n_missed` is the number who arrived while the wait was running to a full
#'   waiting room and were never enrolled. See [set_dtp_queue_size()] for the
#'   room capacity that bounds `n_held` and produces `n_missed`.
#' @export
dtp_wait_events <- function(x, ...) UseMethod("dtp_wait_events")

#' @export
dtp_wait_events.simulations <- function(x, ...) {
  per_rep <- lapply(seq_along(x$fits), function(i) {
    last <- x$fits[[i]][[length(x$fits[[i]])]]
    ev <- last$dtp_wait_events
    if (is.null(ev)) return(NULL)
    if (nrow(ev) == 0L) return(cbind(replicate = integer(0), ev))
    cbind(replicate = rep(i, nrow(ev)), ev)
  })
  out <- dplyr::bind_rows(per_rep)
  if (nrow(out) == 0L) {
    out <- cbind(replicate = integer(0), .wait_events_empty())
  }
  tibble::as_tibble(out)
}

#' @export
dtp_wait_events.simulations_collection <- function(x, ...) {
  dplyr::bind_rows(lapply(x, dtp_wait_events), .id = "design")
}

# -- dtp_wait_summary --------------------------------------------------

#' Per-trial DTP wait summary from simulated trials.
#'
#' Aggregates the raw event log returned by [dtp_wait_events()] into a
#' curated per-replicate (or per-replicate × starting dose) summary.
#'
#' @details
#' With the default `by_dose = FALSE`, one row per replicate with columns:
#'
#' - `num_waits` — how often DTP fired in this trial
#' - `total_wait_time` — sum of wait durations
#' - `wait_fraction` — `total_wait_time / trial_duration`, a single
#'   comparable "DTP overhead" number. The denominator is the full trial
#'   calendar time as returned by [escalation::trial_duration()], which
#'   includes the post-loop "advance to full follow-up" tail. Use the raw
#'   event log if you need a decision-phase-only denominator.
#' - `max_wait_time` — longest single wait (for slot / IRB conversations).
#'   Can exceed `t_max` when an in-wait DLT triggers an extension; see
#'   [dtp_wait_events()].
#' - `num_effective_waits` — number of waits where the post-wait
#'   recommendation was strictly higher than the pre-wait recommendation
#'   (`dose_after > dose_before`). Captures both outright escalations and
#'   waits that recovered from an otherwise unnecessary de-escalation.
#' - `num_missed` — per-trial sum of the event-log `n_missed`: patients who
#'   presented during a wait to a full waiting room and were never enrolled.
#'   This is DTP's accrual cost, and is zero whenever no wait ran longer than
#'   the room capacity allows.
#'
#' With `by_dose = TRUE`, one row per (replicate × `dose_before`).
#' `wait_fraction` (needs a whole-trial denominator) is dropped and a
#' `mean_wait_time` column is added; `max_wait_time` and
#' `num_effective_waits` are retained per dose level.
#'
#' Non-DTP designs — those whose simulation did not attach a wait event log
#' at all — return all-`NA` metric columns so that users can distinguish
#' "design did not run under DTP" from "ran under DTP and never waited"
#' (`num_waits = 0`).
#'
#' On a `simulations_collection` (e.g. the return value of
#' `simulate_compare`), both accessors prepend a `design` column naming
#' the arm.
#'
#' @param x A `simulations` or `simulations_collection` object.
#' @param by_dose If `TRUE`, roll up by `(replicate, dose_before)`;
#'   otherwise one row per replicate. Default `FALSE`.
#' @param ... Unused.
#' @return A tibble. See Details for columns.
#' @export
dtp_wait_summary <- function(x, by_dose = FALSE, ...) {
  UseMethod("dtp_wait_summary")
}

#' @export
dtp_wait_summary.simulations <- function(x, by_dose = FALSE, ...) {
  n_rep <- length(x$fits)
  has_log <- vapply(seq_len(n_rep), function(i) {
    last <- x$fits[[i]][[length(x$fits[[i]])]]
    !is.null(last$dtp_wait_events)
  }, logical(1))

  if (by_dose) {
    cols <- c("replicate", "dose_before", "num_waits", "total_wait_time",
              "mean_wait_time", "max_wait_time", "num_effective_waits",
              "num_missed")
  } else {
    cols <- c("replicate", "num_waits", "total_wait_time", "wait_fraction",
              "max_wait_time", "num_effective_waits", "num_missed")
  }

  # Non-DTP design: all-NA metrics, one row per replicate
  if (!any(has_log)) {
    na_tbl <- tibble::tibble(replicate = seq_len(n_rep))
    for (c in setdiff(cols, "replicate")) na_tbl[[c]] <- NA_real_
    return(na_tbl[, cols, drop = FALSE])
  }

  ev <- dtp_wait_events(x)

  if (by_dose) {
    # One row per (replicate × dose_before) actually observed in the log.
    if (nrow(ev) == 0L) {
      out <- tibble::tibble(
        replicate           = integer(0),
        dose_before         = integer(0),
        num_waits           = integer(0),
        total_wait_time     = numeric(0),
        mean_wait_time      = numeric(0),
        max_wait_time       = numeric(0),
        num_effective_waits = integer(0),
        num_missed          = integer(0)
      )
      return(out)
    }
    out <- dplyr::summarise(
      dplyr::group_by(ev, replicate, dose_before),
      num_waits           = dplyr::n(),
      total_wait_time     = sum(wait_duration),
      mean_wait_time      = mean(wait_duration),
      max_wait_time       = max(wait_duration),
      num_effective_waits = sum(effective),
      num_missed          = sum(n_missed),
      .groups = "drop"
    )
    return(tibble::as_tibble(out[, cols, drop = FALSE]))
  }

  # Default: per-replicate roll-up, with explicit 0 rows for replicates
  # whose event log is empty (DTP ran but never waited).
  if (nrow(ev) > 0L) {
    agg <- dplyr::summarise(
      dplyr::group_by(ev, replicate),
      num_waits           = dplyr::n(),
      total_wait_time     = sum(wait_duration),
      max_wait_time       = max(wait_duration),
      num_effective_waits = sum(effective),
      num_missed          = sum(n_missed),
      .groups = "drop"
    )
  } else {
    agg <- tibble::tibble(
      replicate           = integer(0),
      num_waits           = integer(0),
      total_wait_time     = numeric(0),
      max_wait_time       = numeric(0),
      num_effective_waits = integer(0),
      num_missed          = integer(0)
    )
  }

  base <- tibble::tibble(replicate = seq_len(n_rep))
  out <- dplyr::left_join(base, agg, by = "replicate")
  out$num_waits[is.na(out$num_waits)] <- 0L
  out$total_wait_time[is.na(out$total_wait_time)] <- 0
  out$max_wait_time[is.na(out$max_wait_time)] <- 0
  out$num_effective_waits[is.na(out$num_effective_waits)] <- 0L
  out$num_missed[is.na(out$num_missed)] <- 0L

  # wait_fraction = total_wait_time / per-replicate trial duration. Trust
  # escalation::trial_duration() to return one value per replicate; if it
  # ever doesn't, that's an upstream bug worth surfacing.
  dur <- trial_duration(x)
  out$wait_fraction <- ifelse(dur > 0, out$total_wait_time / dur, NA_real_)

  out[, cols, drop = FALSE]
}

#' @export
dtp_wait_summary.simulations_collection <- function(x, by_dose = FALSE, ...) {
  dplyr::bind_rows(
    lapply(x, dtp_wait_summary, by_dose = by_dose),
    .id = "design"
  )
}

# Silence R CMD check NOTE: undefined globals used via NSE in dplyr calls
utils::globalVariables(c(
  "replicate", "dose_before", "wait_duration", "effective"
))
