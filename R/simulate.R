# DTP-aware simulation for TITE dose-finding designs
#
# Provides a custom simulation function that implements DTP's "wait" mechanism
# so that Monte Carlo simulations capture the trial-duration impact of waiting
# for pending patients' follow-up before dosing the next cohort.
#
# Key additions over escalation's phase1_tite_sim:
# 1. Queue-aware cohort assembly (patients arriving during wait)
# 2. DTP wait block with event-driven DLT handling
# 3. simulation_function dispatch for dtp_selector_factory

# -- Exported: simulation_function dispatch --------------------------------

#' @importFrom escalation simulation_function
#' @export
simulation_function.dtp_selector_factory <- function(selector_factory) {
  phase1_dtp_tite_sim
}

# -- Exported: tite_patient_samples ----------------------------------------

#' Create patient samples with realistic TITE DLT timing.
#'
#' Generates a list of [escalation::PatientSample] objects whose
#' `tox_time` values are drawn from `U(0, max_time)` instead of the
#' default `U(0, 1)`. Pass the result to
#' [escalation::simulate_compare()] or [escalation::simulate_trials()]
#' via `patient_samples`.
#'
#' @param num_sims Number of simulation replicates.
#' @param max_time DLT observation window (e.g. 56 days). Must match the
#'   `max_time` argument passed to `simulate_compare` / `simulate_trials`.
#' @param num_patients Number of latent patients per sample (default 100).
#' @return A list of length `num_sims`, each element a `PatientSample`.
#' @export
tite_patient_samples <- function(num_sims, max_time, num_patients = 100) {
  stopifnot(is.numeric(num_sims), length(num_sims) == 1, num_sims >= 1)
  stopifnot(is.numeric(max_time), length(max_time) == 1, max_time > 0)
  stopifnot(is.numeric(num_patients), length(num_patients) == 1,
            num_patients >= 1)

  lapply(seq_len(num_sims), function(i) {
    escalation::PatientSample$new(
      num_patients = num_patients,
      time_to_tox_func = function() stats::runif(1, 0, max_time)
    )
  })
}

# -- Exported: set_dtp_queue_size ------------------------------------------

#' Set the queue size for DTP simulation.
#'
#' During a DTP wait, patients who arrive are queued (up to `queue_size`)
#' and dosed at the post-wait recommendation. Queued patients replace
#' (part of) the next cohort — they are not additional enrolments.
#'
#' @param dtp_selector_factory A `dtp_selector_factory` created by
#'   [apply_dtp()].
#' @param queue_size Integer >= 0. Default (when not set) equals the
#'   cohort size (inferred at simulation time from the first
#'   `sample_patient_arrivals` call).
#' @return The modified `dtp_selector_factory` (invisibly).
#' @export
set_dtp_queue_size <- function(dtp_selector_factory, queue_size) {
  fac <- .find_dtp_factory(dtp_selector_factory)
  if (is.null(fac)) {
    stop("No dtp_selector_factory found in the factory chain.")
  }
  stopifnot(is.numeric(queue_size), length(queue_size) == 1, queue_size >= 0)
  fac$sim_settings$queue_size <- as.integer(queue_size)
  invisible(dtp_selector_factory)
}

# -- Internal: find DTP factory in chain -----------------------------------

.find_dtp_factory <- function(factory) {
  if (inherits(factory, "dtp_selector_factory")) return(factory)
  if (!is.null(factory$parent)) return(.find_dtp_factory(factory$parent))
  NULL
}

# -- Internal: consume arrivals during wait --------------------------------

.consume_arrivals <- function(time_now, wait_end, queue_size,
                              sample_patient_arrivals, all_data) {
  if (queue_size == 0) return(numeric(0))

  queued_times <- numeric(0)
  cursor <- time_now

  while (length(queued_times) < queue_size && cursor < wait_end) {
    new_pts <- sample_patient_arrivals(all_data)
    arrival_deltas <- cumsum(new_pts$time_delta)
    abs_times <- cursor + arrival_deltas

    for (t in abs_times) {
      if (t > wait_end) break
      if (length(queued_times) >= queue_size) break
      queued_times <- c(queued_times, t)
    }
    cursor <- cursor + sum(new_pts$time_delta)
  }

  queued_times
}

# -- Internal: extend queued arrivals over incremental wait window ----------

.extend_queued_arrivals <- function(queued_arrival_times,
                                    old_wait_end,
                                    new_wait_end,
                                    queue_size,
                                    sample_patient_arrivals,
                                    all_data) {
  if (new_wait_end <= old_wait_end) return(queued_arrival_times)
  remaining_capacity <- queue_size - length(queued_arrival_times)
  if (remaining_capacity <= 0) return(queued_arrival_times)

  added <- .consume_arrivals(
    time_now = old_wait_end,
    wait_end = new_wait_end,
    queue_size = remaining_capacity,
    sample_patient_arrivals = sample_patient_arrivals,
    all_data = all_data
  )
  c(queued_arrival_times, added)
}

# -- Internal: event-driven wait loop --------------------------------------

.advance_wait_with_events <- function(time_now, wait_end,
                                      tite_dose, tite_time, tite_cohort,
                                      base_df, patient_sample,
                                      true_prob_tox, selector_factory,
                                      max_time, get_weight,
                                      queued_arrival_times,
                                      queue_size,
                                      sample_patient_arrivals,
                                      all_data) {
  # Count how many times an in-wait DLT-driven re-fit caused wait_end to be
  # bumped further out. 0 means the wait ran on its initially-projected
  # schedule; >0 means the wait was chained by mid-wait events.
  num_extensions <- 0L
  while (time_now < wait_end) {
    # Find earliest unobserved DLT from enrolled patients
    n_enrolled <- length(tite_dose)
    if (n_enrolled == 0) {
      time_now <- wait_end
      break
    }

    # For each enrolled patient, check if they'd have a DLT by wait_end
    follow_up_at_wait_end <- wait_end - tite_time
    tox_at_end <- patient_sample$get_patient_tox(
      i = seq_len(n_enrolled),
      prob_tox = true_prob_tox[tite_dose],
      time = follow_up_at_wait_end
    )

    # Check current tox status
    follow_up_now <- time_now - tite_time
    tox_now <- patient_sample$get_patient_tox(
      i = seq_len(n_enrolled),
      prob_tox = true_prob_tox[tite_dose],
      time = pmax(follow_up_now, 0)
    )

    # Find patients who will have a NEW DLT between now and wait_end
    new_dlt <- which(tox_at_end == 1 & tox_now == 0)

    if (length(new_dlt) == 0) {
      time_now <- wait_end
      break
    }

    # Find the earliest new DLT time
    # DLT occurs at: tite_time[i] + tox_time[i] from PatientSample
    # We need to find the absolute time of each new DLT
    dlt_abs_times <- numeric(length(new_dlt))
    for (k in seq_along(new_dlt)) {
      idx <- new_dlt[k]
      # Binary search for the exact DLT time
      # We know DLT happens between time_now and wait_end for this patient
      lo <- pmax(time_now - tite_time[idx], 0)
      hi <- wait_end - tite_time[idx]
      # Find minimum follow-up at which DLT is observed
      while (hi - lo > 0.5) {
        mid <- (lo + hi) / 2
        tox_mid <- patient_sample$get_patient_tox(
          i = idx,
          prob_tox = true_prob_tox[tite_dose[idx]],
          time = mid
        )
        if (tox_mid == 1) {
          hi <- mid
        } else {
          lo <- mid
        }
      }
      dlt_abs_times[k] <- tite_time[idx] + hi
    }

    earliest_dlt_time <- min(dlt_abs_times)
    if (earliest_dlt_time >= wait_end) {
      time_now <- wait_end
      break
    }

    # Advance to the DLT event
    time_now <- earliest_dlt_time

    # Re-compute tox/weights and re-fit
    tite_data <- data.frame(
      cohort = tite_cohort,
      dose = tite_dose,
      time = tite_time
    )
    tite_data$tox <- patient_sample$get_patient_tox(
      i = seq_len(n_enrolled),
      prob_tox = true_prob_tox[tite_dose],
      time = time_now - tite_time
    )
    tite_data$weight <- get_weight(
      now_time = time_now,
      recruited_time = tite_time,
      tox = tite_data$tox,
      max_time = max_time
    )

    all_data <- dplyr::bind_rows(base_df, tite_data)
    all_data$patient <- seq_len(nrow(all_data))

    fit <- selector_factory |> fit(all_data)
    next_dose <- recommended_dose(fit)

    # Re-evaluate DTP decision
    if (!dtp_should_wait(fit) || !continue(fit) || is.na(next_dose)) {
      ended_by <- if (!continue(fit) || is.na(next_dose)) {
        "stopped"
      } else {
        "dtp_decision"
      }
      return(list(
        time_now = time_now, fit = fit, next_dose = next_dose,
        all_data = all_data, tite_data = tite_data,
        ended_by = ended_by,
        num_extensions = num_extensions,
        queued_arrival_times = queued_arrival_times
      ))
    }

    # Update wait_end based on new DTP assessment. This is an extension —
    # the wait will run past its previous wait_end.
    old_wait_end <- wait_end
    wait_end <- time_now + dtp_wait_time(fit)
    queued_arrival_times <- .extend_queued_arrivals(
      queued_arrival_times = queued_arrival_times,
      old_wait_end = old_wait_end,
      new_wait_end = wait_end,
      queue_size = queue_size,
      sample_patient_arrivals = sample_patient_arrivals,
      all_data = all_data
    )
    num_extensions <- num_extensions + 1L
  }

  # Reached wait_end without interruption — do final re-fit
  n_enrolled <- length(tite_dose)
  tite_data <- data.frame(
    cohort = tite_cohort,
    dose = tite_dose,
    time = tite_time
  )
  if (n_enrolled > 0) {
    tite_data$tox <- patient_sample$get_patient_tox(
      i = seq_len(n_enrolled),
      prob_tox = true_prob_tox[tite_dose],
      time = time_now - tite_time
    )
  } else {
    tite_data$tox <- integer(0)
  }
  tite_data$weight <- get_weight(
    now_time = time_now,
    recruited_time = tite_time,
    tox = tite_data$tox,
    max_time = max_time
  )

  all_data <- dplyr::bind_rows(base_df, tite_data)
  all_data$patient <- seq_len(nrow(all_data))

  fit <- selector_factory |> fit(all_data)
  next_dose <- recommended_dose(fit)

  list(
    time_now = time_now, fit = fit, next_dose = next_dose,
    all_data = all_data, tite_data = tite_data,
    ended_by = "wait_end",
    num_extensions = num_extensions,
    queued_arrival_times = queued_arrival_times
  )
}

# -- Internal: main simulation function ------------------------------------
#
# `min_fup_time` is the minimum follow-up enforced between a cohort being
# fully enrolled and the model being (re)fit for its dose decision. It is
# passed through by `simulate_trials(..., min_fup_time = )`.
#
# Implication of the default `min_fup_time = 0`: when a cohort is dosed
# straight from the DTP queue at the wait-end instant, all of its patients
# have zero follow-up (weight 0) at the immediately-following fit, so they
# contribute no information and that cohort's dose decision (including whether
# DTP waits again) is a foregone conclusion determined by the earlier
# patients. To have each decision rest on real follow-up, set `min_fup_time`
# to roughly the fixed accrual interval (the per-patient inter-arrival time),
# which delays each fit until the freshest patients have accrued weight.
# Note this is a global timing change that shifts operating characteristics.

phase1_dtp_tite_sim <- function(selector_factory, true_prob_tox,
                                patient_sample = NULL,
                                sample_patient_arrivals = function(df) {
                                  escalation::cohorts_of_n(n = 1,
                                                           mean_time_delta = 1)
                                },
                                previous_outcomes = "",
                                next_dose = NULL,
                                time_now = NULL,
                                max_time,
                                min_fup_time = 0,
                                get_weight = escalation::linear_follow_up_weight,
                                i_like_big_trials = FALSE,
                                return_all_fits = FALSE) {

  if (length(max_time) > 1 | max_time <= 0) {
    stop("max_time should be a strictly positive scalar.")
  }

  if (is.null(patient_sample)) {
    patient_sample <- escalation::PatientSample$new(
      time_to_tox_func = function() stats::runif(1, 0, max_time)
    )
  }

  # -- Parse previous outcomes (matches phase1_tite_sim) -------------------
  if (is.character(previous_outcomes)) {
    base_df <- escalation::parse_phase1_outcomes(previous_outcomes,
                                                  as_list = FALSE)
  } else if (is.data.frame(previous_outcomes)) {
    base_df <- previous_outcomes
  } else {
    base_df <- escalation::parse_phase1_outcomes("", as_list = FALSE)
  }

  all_data <- base_df

  if (nrow(base_df) > 0) {
    base_df$weight <- 1
    if (is.null(time_now)) {
      if ("time" %in% colnames(base_df)) {
        time_now <- max(base_df$time)
      }
    }
  }

  if (is.null(time_now)) time_now <- 0

  cohort <- base_df$cohort
  next_cohort <- ifelse(length(cohort) > 0, max(cohort) + 1, 1)

  i <- 1
  max_i <- 30

  fit_obj <- selector_factory |> fit(base_df)
  if (is.null(next_dose)) next_dose <- recommended_dose(fit_obj)

  fits <- list()
  fits[[1]] <- list(.depth = i, time = time_now, fit = fit_obj)

  tite_dose <- integer(0)
  tite_cohort <- integer(0)
  tite_time <- numeric(0)

  # -- DTP-specific setup --------------------------------------------------
  dtp_fac <- .find_dtp_factory(selector_factory)
  t_max <- if (!is.null(dtp_fac)) dtp_fac$t_max else 0
  obswin <- if (!is.null(dtp_fac)) dtp_fac$obswin else max_time

  # Queue: filled during DTP waits, consumed at next cohort assembly
  queued_arrival_times <- numeric(0)
  cohort_size <- NULL  # inferred from first sample_patient_arrivals call
  queue_size <- if (!is.null(dtp_fac)) dtp_fac$sim_settings$queue_size else NULL
  # queue_size NULL means "default to cohort_size" — resolved after first call

  # DTP wait event log: one record per wait fired in this trial. Populated
  # only when a dtp_selector_factory is in the chain; stays empty otherwise.
  dtp_wait_events <- list()

  # -- Main loop -----------------------------------------------------------
  while (continue(fit_obj) & !is.na(next_dose) &
         (i_like_big_trials | i < max_i)) {

    # ── ASSEMBLE NEXT COHORT ──────────────────────────────────────────────
    # Produces: recruit_abs_times (absolute recruitment times for new patients)
    #           time_span (time from first to last patient arrival, for clock advance)
    if (!is.null(cohort_size) && length(queued_arrival_times) >= cohort_size) {
      # Queue already has a full cohort — they're already here, dose them now
      queued_arrival_times <- queued_arrival_times[-seq_len(cohort_size)]
      n_new_pts <- cohort_size
      recruit_abs_times <- rep(time_now, cohort_size)
      time_span <- 0

    } else if (length(queued_arrival_times) > 0) {
      # Partial queue — queued patients dosed now, fill remaining from arrivals
      n_queued <- length(queued_arrival_times)
      remaining <- cohort_size - n_queued
      fresh <- sample_patient_arrivals(all_data)
      n_fresh <- min(remaining, nrow(fresh))
      fresh_cumdeltas <- cumsum(fresh$time_delta[seq_len(n_fresh)])

      recruit_abs_times <- c(
        rep(time_now, n_queued),
        time_now + fresh_cumdeltas
      )
      time_span <- fresh_cumdeltas[n_fresh]
      n_new_pts <- n_queued + n_fresh
      queued_arrival_times <- numeric(0)

    } else {
      # No queue — standard fresh cohort
      new_pts <- sample_patient_arrivals(all_data)
      n_new_pts <- nrow(new_pts)
      cum_deltas <- cumsum(new_pts$time_delta)
      recruit_abs_times <- time_now + cum_deltas
      time_span <- sum(new_pts$time_delta)

      # Infer cohort_size from first call
      if (is.null(cohort_size)) {
        cohort_size <- n_new_pts
        if (is.null(queue_size)) queue_size <- cohort_size
        # Validate queue_size
        if (!is.null(dtp_fac) &&
            !is.null(dtp_fac$sim_settings$queue_size)) {
          if (dtp_fac$sim_settings$queue_size > cohort_size) {
            stop("queue_size (", dtp_fac$sim_settings$queue_size,
                 ") must be <= cohort_size (", cohort_size, ")")
          }
          queue_size <- dtp_fac$sim_settings$queue_size
        }
      }
    }

    # ── ENROLL + FIT ──────────────────────────────────────────────────────
    new_dose <- rep(next_dose, n_new_pts)
    new_cohort <- rep(next_cohort, n_new_pts)

    tite_cohort <- c(tite_cohort, new_cohort)
    tite_dose <- c(tite_dose, new_dose)
    tite_time <- c(tite_time, recruit_abs_times)

    tite_data <- data.frame(
      cohort = tite_cohort,
      dose = tite_dose,
      time = tite_time
    )

    time_now <- time_now + time_span + min_fup_time

    if (nrow(tite_data) > 0) {
      tite_data$tox <- patient_sample$get_patient_tox(
        i = seq_along(tite_dose),
        prob_tox = true_prob_tox[tite_dose],
        time = time_now - tite_time
      )
    } else {
      tite_data$tox <- integer(0)
    }

    tite_data$weight <- get_weight(
      now_time = time_now,
      recruited_time = tite_data$time,
      tox = tite_data$tox,
      max_time = max_time
    )

    all_data <- dplyr::bind_rows(base_df, tite_data)
    all_data$patient <- seq_len(nrow(all_data))

    fit_obj <- selector_factory |> fit(all_data)
    next_dose <- recommended_dose(fit_obj)

    i <- i + 1
    next_cohort <- next_cohort + 1
    fits[[i]] <- list(.depth = i, time = time_now, fit = fit_obj)

    # ── DTP DECISION ──────────────────────────────────────────────────────
    if (t_max > 0 && dtp_should_wait(fit_obj) &&
        continue(fit_obj) && !is.na(next_dose)) {

      wait_time_val <- dtp_wait_time(fit_obj)
      wait_end <- time_now + wait_time_val

      # Capture pre-wait state for the event log
      time_before <- time_now
      dose_before <- next_dose
      fit_before  <- fit_obj

      # A. Consume arrivals during wait window
      queued_arrival_times <- .consume_arrivals(
        time_now = time_now,
        wait_end = wait_end,
        queue_size = queue_size,
        sample_patient_arrivals = sample_patient_arrivals,
        all_data = all_data
      )

      # B. Event-driven wait (DLTs from enrolled patients)
      wait_result <- .advance_wait_with_events(
        time_now = time_now,
        wait_end = wait_end,
        tite_dose = tite_dose,
        tite_time = tite_time,
        tite_cohort = tite_cohort,
        base_df = base_df,
        patient_sample = patient_sample,
        true_prob_tox = true_prob_tox,
        selector_factory = selector_factory,
        max_time = max_time,
        get_weight = get_weight,
        queued_arrival_times = queued_arrival_times,
        queue_size = queue_size,
        sample_patient_arrivals = sample_patient_arrivals,
        all_data = all_data
      )

      time_now <- wait_result$time_now
      fit_obj <- wait_result$fit
      next_dose <- wait_result$next_dose
      all_data <- wait_result$all_data
      queued_arrival_times <- wait_result$queued_arrival_times
      # An in-wait DLT can end the wait before its projected wait_end. The
      # queue was pre-filled (.consume_arrivals / .extend_queued_arrivals)
      # against that projected window, so drop any queued arrival dated after
      # the actual wait end (time_now) — those patients never really arrived.
      queued_arrival_times <-
        queued_arrival_times[queued_arrival_times <= time_now]
      queue_this_wait <- length(queued_arrival_times)

      dose_delta <- if (is.na(next_dose) || is.na(dose_before)) {
        NA_integer_
      } else {
        as.integer(next_dose - dose_before)
      }
      projected_dose <- as.integer(dtp_projected_dose(fit_before))
      # `effective`: the wait moved the recommendation strictly upward
      # compared to the pre-wait recommendation. With DTP's invariant
      # `dose_before <= dose_current`, this captures both outright
      # escalations and recovered de-escalations. (`projection_realized`
      # was previously a separate flag here; dropped because for monotone
      # TITE-CRM/BOIN with integer dose levels and projected_dose ≈
      # dose_before + 1, it is row-for-row identical to `effective`. The
      # raw `projected_dose` column is retained for audit.)
      effective <- !is.na(dose_before) && !is.na(next_dose) &&
        as.integer(next_dose) > as.integer(dose_before)

      dtp_wait_events[[length(dtp_wait_events) + 1L]] <- list(
        cohort_idx     = as.integer(next_cohort - 1L),
        dose_before    = as.integer(dose_before),
        time_in        = time_before,
        time_out       = time_now,
        wait_duration  = time_now - time_before,
        projected_dose = projected_dose,
        dose_after     = as.integer(next_dose),
        ended_by       = wait_result$ended_by,
        dose_delta     = dose_delta,
        effective      = isTRUE(effective),
        num_extensions = as.integer(wait_result$num_extensions),
        queue_size     = as.integer(queue_this_wait)
      )
    }
  }

  # -- Safety valve --------------------------------------------------------
  if (!i_like_big_trials & i >= max_i) {
    warning(paste("Simulation stopped because max depth reached.",
                  "Set 'i_like_big_trials = TRUE' to avoid this constraint. "))
  }

  # -- Final analysis (advance to full follow-up) --------------------------
  if (!is.na(next_dose)) {
    time_now <- time_now + max_time - min_fup_time
    tite_data <- data.frame(
      cohort = tite_cohort,
      dose = tite_dose,
      time = tite_time
    )
    if (nrow(tite_data) > 0) {
      tite_data$tox <- patient_sample$get_patient_tox(
        i = seq_along(tite_dose),
        prob_tox = true_prob_tox[tite_dose],
        time = time_now - tite_time
      )
    } else {
      tite_data$tox <- integer(0)
    }
    tite_data$weight <- get_weight(
      now_time = time_now,
      recruited_time = tite_data$time,
      tox = tite_data$tox,
      max_time = max_time
    )
    all_data <- dplyr::bind_rows(base_df, tite_data)
    all_data$patient <- seq_len(nrow(all_data))
    fit_obj <- selector_factory |> fit(all_data)
    i <- i + 1
    fits[[i]] <- list(.depth = i, time = time_now, fit = fit_obj)
  }

  # Attach the DTP wait event log to the final fit entry. Escalation's
  # built-in extractors only inspect `$fit` / `$time` on tail(.x, 1)[[1]]
  # so adding a named field is non-breaking. Empty tibble (right columns,
  # zero rows) when no waits fired — downstream accessors rely on this.
  last_idx <- length(fits)
  fits[[last_idx]]$dtp_wait_events <- .build_wait_events_tibble(dtp_wait_events)

  if (return_all_fits) {
    return(fits)
  } else {
    return(fits[length(fits)])
  }
}

# -- Internal: assemble wait event log tibble ------------------------------

.wait_events_cols <- c("cohort_idx", "dose_before", "time_in", "time_out",
                       "wait_duration", "projected_dose", "dose_after",
                       "ended_by", "dose_delta", "effective",
                       "num_extensions", "queue_size")

.wait_events_empty <- function() {
  out <- data.frame(
    cohort_idx     = integer(0),
    dose_before    = integer(0),
    time_in        = numeric(0),
    time_out       = numeric(0),
    wait_duration  = numeric(0),
    projected_dose = integer(0),
    dose_after     = integer(0),
    ended_by       = factor(character(0),
                            levels = c("wait_end", "dtp_decision", "stopped")),
    dose_delta     = integer(0),
    effective      = logical(0),
    num_extensions = integer(0),
    queue_size     = integer(0),
    stringsAsFactors = FALSE
  )
  out
}

.build_wait_events_tibble <- function(records) {
  if (length(records) == 0L) return(.wait_events_empty())
  df <- dplyr::bind_rows(lapply(records, function(r) {
    r$ended_by <- as.character(r$ended_by)
    as.data.frame(r, stringsAsFactors = FALSE)
  }))
  df$ended_by <- factor(df$ended_by,
                         levels = c("wait_end", "dtp_decision", "stopped"))
  df[, .wait_events_cols, drop = FALSE]
}
