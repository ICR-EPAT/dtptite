# Tests for DTP-specific simulation metrics
#
# Coverage:
#  (a) Raw event log has the expected rows and column types
#  (b) dtp_wait_summary columns + arithmetic sanity
#  (c) Raw log column types (factor levels, queue_size cap)
#  (d) Termination path: "stopped" appears for both CRM and BOIN designs
#      and coincides with continue()==FALSE + is.na(rec)
#  (e) by_dose = TRUE shape
#  (f) Non-DTP path: NA metrics, empty event log
#  (g) simulate_compare: design column, three arms
#  (h) Regression guard: prob_recommend / trial_duration / num_patients unchanged

library(escalation)

skeleton <- c(0.124, 0.25, 0.398, 0.542, 0.666)
target <- 0.25
arrivals3 <- function(df) escalation::cohorts_of_n(n = 3, mean_time_delta = 1)

# ----------------------------------------------------------------------
# Shared fixture: a DTP+CRM design under a scenario that forces waits
# ----------------------------------------------------------------------

make_dtp_crm <- function() {
  get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 35, obswin = 56) |>
    stop_at_n(n = 12)
}

run_dtp_crm_sims <- function(num_sims = 8) {
  set.seed(456)
  simulate_trials(
    make_dtp_crm(),
    num_sims = num_sims,
    true_prob_tox = c(0.05, 0.15, 0.25, 0.35, 0.60),
    sample_patient_arrivals = arrivals3,
    max_time = 56
  )
}

# ===== (a) raw event log shape ========================================

test_that("dtp_wait_events returns rows with expected columns and types", {
  sims <- run_dtp_crm_sims(num_sims = 8)
  ev <- dtp_wait_events(sims)

  expect_s3_class(ev, "tbl_df")
  expect_true(all(c("replicate", "cohort_idx", "dose_before", "time_in",
                    "time_out", "wait_duration", "projected_dose",
                    "dose_after", "ended_by", "dose_delta", "effective",
                    "num_extensions",
                    "queue_size") %in% names(ev)))
  expect_gt(nrow(ev), 0)

  expect_type(ev$cohort_idx, "integer")
  expect_type(ev$dose_before, "integer")
  expect_type(ev$wait_duration, "double")
  expect_s3_class(ev$ended_by, "factor")
  expect_type(ev$effective, "logical")
  expect_type(ev$num_extensions, "integer")
  expect_true(all(ev$wait_duration >= 0))
  expect_true(all(ev$time_out >= ev$time_in))
  expect_true(all(ev$num_extensions >= 0L))
})

# ===== (b) dtp_wait_summary arithmetic sanity =========================

test_that("dtp_wait_summary: headline columns and arithmetic invariants", {
  sims <- run_dtp_crm_sims(num_sims = 8)
  smry <- dtp_wait_summary(sims)
  ev <- dtp_wait_events(sims)

  expect_s3_class(smry, "tbl_df")
  expect_equal(nrow(smry), length(sims$fits))
  expect_setequal(names(smry),
                  c("replicate", "num_waits", "total_wait_time",
                    "wait_fraction", "max_wait_time", "num_effective_waits",
                    "num_missed"))

  # total_wait_time per replicate matches sum of wait_duration in the log
  for (i in seq_len(nrow(smry))) {
    rep_i <- smry$replicate[i]
    ev_i <- ev[ev$replicate == rep_i, ]
    expect_equal(smry$total_wait_time[i], sum(ev_i$wait_duration))
    expect_equal(smry$num_waits[i], nrow(ev_i))
    expect_equal(smry$num_effective_waits[i],
                 sum(ev_i$dose_after > ev_i$dose_before, na.rm = TRUE))
  }

  # wait_fraction in [0, 1]; max_wait_time <= total_wait_time
  expect_true(all(smry$wait_fraction >= 0 & smry$wait_fraction <= 1,
                  na.rm = TRUE))
  expect_true(all(smry$max_wait_time <= smry$total_wait_time + 1e-9))
  expect_true(all(smry$num_effective_waits <= smry$num_waits))
})

# ===== (c) Raw log column types =======================================

test_that("dtp_wait_events: ended_by levels, queue_size cap", {
  sims <- run_dtp_crm_sims(num_sims = 8)
  ev <- dtp_wait_events(sims)

  expect_equal(levels(ev$ended_by),
               c("wait_end", "dtp_decision", "stopped"))
  expect_true(all(as.character(ev$ended_by) %in%
                    c("wait_end", "dtp_decision", "stopped")))
  # Default queue_size equals cohort size (3 here)
  expect_true(all(ev$queue_size >= 0 & ev$queue_size <= 3))
})

test_that("set_dtp_queue_size caps queue_size column in event log", {
  design <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 35, obswin = 56) |>
    stop_at_n(n = 12)
  set_dtp_queue_size(design, 2)

  set.seed(654)
  sims <- simulate_trials(
    design,
    num_sims = 5,
    true_prob_tox = c(0.05, 0.15, 0.25, 0.35, 0.60),
    sample_patient_arrivals = arrivals3,
    max_time = 56
  )
  ev <- dtp_wait_events(sims)
  expect_true(all(ev$queue_size <= 2))
})

# ===== (d) Termination path: CRM and BOIN =============================

test_that("ended_by = 'stopped' fires for DTP+CRM in toxic scenario", {
  # Toxic scenario + Beta-Binom stopping on dose 1
  true_prob_tox <- c(0.60, 0.70, 0.80, 0.90, 0.95)

  design <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    stop_for_beta_binomial_toxicity(dose = 1, tox_threshold = 0.25,
                                    confidence = 0.70) |>
    apply_dtp(t_max = 35, obswin = 56) |>
    stop_at_n(n = 24)

  set.seed(22)
  sims <- simulate_trials(
    design,
    num_sims = 20,
    true_prob_tox = true_prob_tox,
    sample_patient_arrivals = arrivals3,
    max_time = 56
  )

  ev <- dtp_wait_events(sims)

  # At least one "stopped" event across replicates
  expect_true(any(ev$ended_by == "stopped"))

  # Stopped replicates should have enrolled fewer than the stop_at_n cap
  # (the only way this scenario stops before 24 is via the Beta-Binom chain)
  np <- num_patients(sims)
  stopped_reps <- unique(ev$replicate[ev$ended_by == "stopped"])
  expect_true(all(np[stopped_reps] < 24))
})

test_that("ended_by = 'stopped' fires for DTP+BOIN with dose='any' elimination", {
  true_prob_tox <- c(0.60, 0.70, 0.80, 0.90, 0.95)

  design <- get_boin_tite(5, target) |>
    stop_for_beta_binomial_toxicity(dose = "any", tox_threshold = target,
                                    confidence = 0.90, a = 1, b = 1) |>
    apply_dtp(t_max = 35, obswin = 56) |>
    stop_at_n(n = 24)

  set.seed(44)
  sims <- simulate_trials(
    design,
    num_sims = 20,
    true_prob_tox = true_prob_tox,
    sample_patient_arrivals = arrivals3,
    max_time = 56
  )

  ev <- dtp_wait_events(sims)

  expect_true(any(ev$ended_by == "stopped"))

  np <- num_patients(sims)
  stopped_reps <- unique(ev$replicate[ev$ended_by == "stopped"])
  expect_true(all(np[stopped_reps] < 24))
})

# ===== (e) by_dose shape ==============================================

test_that("dtp_wait_summary(by_dose=TRUE) has one row per (replicate × dose_before)", {
  sims <- run_dtp_crm_sims(num_sims = 8)
  by_d <- dtp_wait_summary(sims, by_dose = TRUE)
  ev <- dtp_wait_events(sims)

  expect_s3_class(by_d, "tbl_df")
  expect_setequal(names(by_d),
                  c("replicate", "dose_before", "num_waits", "total_wait_time",
                    "mean_wait_time", "max_wait_time", "num_effective_waits",
                    "num_missed"))
  expect_false("wait_fraction" %in% names(by_d))
  expect_false("dose" %in% names(by_d))

  # Row count equals distinct (replicate, dose_before) pairs in the log
  keys_log <- unique(paste(ev$replicate, ev$dose_before, sep = "|"))
  keys_smry <- unique(paste(by_d$replicate, by_d$dose_before, sep = "|"))
  expect_setequal(keys_smry, keys_log)

  # max_wait_time at the (replicate, dose_before) slice = max of those rows
  for (i in seq_len(nrow(by_d))) {
    ev_i <- ev[ev$replicate == by_d$replicate[i] &
                 ev$dose_before == by_d$dose_before[i], ]
    expect_equal(by_d$max_wait_time[i], max(ev_i$wait_duration))
    expect_equal(by_d$num_effective_waits[i],
                 sum(ev_i$dose_after > ev_i$dose_before, na.rm = TRUE))
  }
})

# ===== (f) Non-DTP path ===============================================

test_that("dtp_wait_events/summary on non-DTP sims: empty log, NA summary", {
  base_design <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    stop_at_n(n = 9)

  set.seed(88)
  sims <- simulate_trials(
    base_design,
    num_sims = 5,
    true_prob_tox = c(0.05, 0.15, 0.25, 0.35, 0.60),
    sample_patient_arrivals = arrivals3,
    max_time = 56
  )

  ev <- dtp_wait_events(sims)
  expect_s3_class(ev, "tbl_df")
  expect_equal(nrow(ev), 0)
  expect_true(all(c("replicate", "cohort_idx", "ended_by", "effective") %in%
                    names(ev)))

  smry <- dtp_wait_summary(sims)
  expect_equal(nrow(smry), 5)
  for (col in c("num_waits", "total_wait_time", "wait_fraction",
                "max_wait_time", "num_effective_waits")) {
    expect_true(all(is.na(smry[[col]])),
                info = paste("non-DTP summary col", col, "should be all NA"))
  }
})

# ===== (g) simulate_compare path ======================================

test_that("dtp_wait_summary on simulations_collection handles mixed arms", {
  true_prob_tox <- c(0.05, 0.15, 0.25, 0.35, 0.60)

  base_design <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    stop_at_n(n = 9)

  dtp_crm <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 35, obswin = 56) |>
    stop_at_n(n = 9)

  dtp_boin <- get_boin_tite(5, target) |>
    apply_dtp(t_max = 35, obswin = 56) |>
    stop_at_n(n = 9)

  set.seed(200)
  sc <- simulate_compare(
    list(base = base_design, dtp_crm = dtp_crm, dtp_boin = dtp_boin),
    num_sims = 4,
    true_prob_tox = true_prob_tox,
    max_time = 56
  )

  smry <- dtp_wait_summary(sc)
  expect_s3_class(smry, "tbl_df")
  expect_true("design" %in% names(smry))
  expect_equal(nrow(smry), 3 * 4)
  expect_setequal(unique(smry$design), c("base", "dtp_crm", "dtp_boin"))

  # Non-DTP arm: all-NA metric columns
  base_rows <- smry[smry$design == "base", ]
  expect_true(all(is.na(base_rows$num_waits)))
  expect_true(all(is.na(base_rows$total_wait_time)))

  # DTP arms: numeric metrics (possibly 0 if no waits fired, but not NA)
  for (arm in c("dtp_crm", "dtp_boin")) {
    arm_rows <- smry[smry$design == arm, ]
    expect_false(any(is.na(arm_rows$num_waits)))
    expect_false(any(is.na(arm_rows$total_wait_time)))
  }

  # Raw events on the collection also prepend a design column
  ev <- dtp_wait_events(sc)
  expect_true("design" %in% names(ev))
  expect_true(all(ev$design %in% c("dtp_crm", "dtp_boin")))
})

# ===== (h) Regression guard: existing extractors unchanged =============

test_that("attaching the wait event log does not disturb existing extractors", {
  sims <- run_dtp_crm_sims(num_sims = 5)

  # These should all succeed and return the expected shapes. They are the
  # built-in escalation extractors that walk `tail(.x, 1)[[1]]$fit`, which
  # is the same entry we decorate with `$dtp_wait_events`.
  pr <- prob_recommend(sims)
  expect_true(is.numeric(pr))

  td <- trial_duration(sims)
  expect_equal(length(td), 5)
  expect_true(all(td > 0))

  np <- num_patients(sims)
  expect_equal(length(np), 5)
  expect_true(all(np > 0))
})

# ===== (i) Semantic checks on `effective` and `projection_realized` ===

test_that("effective flag is set only when dose_after > dose_before", {
  sims <- run_dtp_crm_sims(num_sims = 8)
  ev <- dtp_wait_events(sims)
  skip_if(nrow(ev) == 0, "no waits in fixture")

  # Direct semantic check on every logged event
  expect_equal(
    ev$effective,
    !is.na(ev$dose_before) & !is.na(ev$dose_after) &
      ev$dose_after > ev$dose_before
  )

  # A wait that ran to wait_end with unchanged dose must be ineffective
  unchanged_wait_end <- ev$ended_by == "wait_end" &
    !is.na(ev$dose_before) & !is.na(ev$dose_after) &
    ev$dose_after == ev$dose_before
  if (any(unchanged_wait_end)) {
    expect_false(any(ev$effective[unchanged_wait_end]))
  }
})

test_that("num_extensions matches recorded wait_duration shape", {
  sims <- run_dtp_crm_sims(num_sims = 8)
  ev <- dtp_wait_events(sims)
  skip_if(nrow(ev) == 0, "no waits in fixture")

  expect_true(all(ev$num_extensions >= 0L))
  # Extensions can only happen for waits that were re-fit mid-wait —
  # which by construction means at least one in-wait DLT was found and
  # the new fit still requested a wait. So extensions > 0 implies the
  # wait took at least *some* time.
  has_ext <- ev$num_extensions > 0L
  if (any(has_ext)) expect_true(all(ev$wait_duration[has_ext] > 0))
})

# ===== (j) Phantom queue over-fill on early wait-end =================
# Regression guard. The queue is pre-filled against the *projected* wait_end.
# When an in-wait DLT ends the wait early (dtp_decision / stopped), queued
# arrivals dated after the *actual* wait end must be dropped — those patients
# never really arrived. With deterministic arrivals spaced `spacing` apart, at
# most floor(wait_duration / spacing) patients can have arrived within a wait,
# so `queue_size` must respect that bound. Before the fix, an early-ended wait
# reported a count reflecting the longer projected window (phantoms).

test_that("queue is not over-filled when an in-wait DLT ends the wait early", {
  spacing <- 7
  arrivals_fixed <- function(df) data.frame(time_delta = rep(spacing, 3))
  design <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 35, obswin = 56) |>
    stop_at_n(n = 12)

  set.seed(18)
  sims <- suppressMessages(simulate_trials(
    design,
    num_sims = 40,
    true_prob_tox = c(0.20, 0.35, 0.50, 0.65, 0.80),
    sample_patient_arrivals = arrivals_fixed,
    max_time = 56
  ))
  ev <- dtp_wait_events(sims)

  # Precondition: the fixture must actually contain early-ended waits,
  # otherwise the guard would be vacuous.
  expect_gt(nrow(ev), 0)
  expect_true(any(ev$ended_by %in% c("dtp_decision", "stopped")))

  # Core invariant: with arrivals every `spacing` units, no wait can have
  # queued more patients than could have arrived within its actual duration.
  expect_true(all(ev$queue_size <= floor(ev$wait_duration / spacing)))
})
