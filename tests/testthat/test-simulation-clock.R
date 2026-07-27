# Timing semantics of phase1_dtp_tite_sim (issue #32).
#
# The pool of available patients and the trial clock are separate. A cohort
# opens when its first patient arrives; the model is fit then, before anyone in
# that cohort is dosed; DTP moves the earliest dosing time, never the update
# time. These tests pin that down with arithmetic that can be checked by hand.

skeleton <- c(0.140, 0.250, 0.376, 0.502, 0.615)
target <- 0.25
obswin <- 56

# Deterministic design/patient fixtures ------------------------------------

design_at <- function(t_max, n) {
  d <- escalation::get_dfcrm_tite(skeleton = skeleton, target = target,
                                  model = "empiric", scale = 0.759) |>
    escalation::dont_skip_doses(when_deescalating = TRUE) |>
    escalation::stop_at_n(n = n)
  apply_dtp(d, t_max = t_max, obswin = obswin, verbose = FALSE)
}

# No patient ever has a DLT.
sample_no_dlt <- function(n = 60) {
  escalation::PatientSample$new(num_patients = n,
                                time_to_tox_func = function() 1e9)
}

# Exactly the listed patients have a DLT, always 21 days after dosing.
# Toxicity is dose-independent so the fixture is fully controlled by tox_u.
sample_dlt_at <- function(who, n = 60) {
  p <- escalation::PatientSample$new(num_patients = n,
                                     time_to_tox_func = function() 21)
  p$tox_u <- rep(0.99, n)
  p$tox_u[who] <- 0.01
  p
}

arrivals_every <- function(gap, cohort_size) {
  force(gap); force(cohort_size)
  function(df) data.frame(time_delta = rep(gap, cohort_size))
}

update_times <- function(fits) vapply(fits, function(x) x$time, numeric(1))

dosing_frame <- function(fit) {
  while (!is.null(fit)) {
    d <- fit$df
    if (!is.null(d) && all(c("dose", "time") %in% names(d))) {
      return(as.data.frame(d))
    }
    fit <- fit$parent
  }
  NULL
}

run <- function(t_max, n, cohort_size, ps, min_fup_time = 0, gap = 14,
                tp = c(0.05, 0.10, 0.15, 0.20, 0.25)) {
  phase1_dtp_tite_sim(
    design_at(t_max, n), true_prob_tox = tp, patient_sample = ps,
    sample_patient_arrivals = arrivals_every(gap, cohort_size),
    max_time = obswin, min_fup_time = min_fup_time, next_dose = 1,
    return_all_fits = TRUE
  )
}

# --- t_max = 0 reduces to accrual-driven TITE-CRM --------------------------

test_that("cohort size 1: the model updates at every arrival", {
  f <- run(t_max = 0, n = 6, cohort_size = 1, ps = sample_no_dlt())

  # t = 0 is escalation's starting recommendation; then one update per arrival.
  expect_equal(head(update_times(f), 7), c(0, 14, 28, 42, 56, 70, 84))
  expect_equal(dosing_frame(f[[length(f)]]$fit)$time,
               c(14, 28, 42, 56, 70, 84))
})

test_that("cohort size 3: one update per cohort, at its first arrival", {
  f <- run(t_max = 0, n = 9, cohort_size = 3, ps = sample_no_dlt())

  expect_equal(head(update_times(f), 4), c(0, 14, 56, 98))
  # Members are dosed as they arrive, not batched.
  expect_equal(dosing_frame(f[[length(f)]]$fit)$time,
               c(14, 28, 42, 56, 70, 84, 98, 112, 126))
})

test_that("the update decides the dose of the cohort it opens", {
  f <- run(t_max = 0, n = 6, cohort_size = 1, ps = sample_no_dlt())
  d <- dosing_frame(f[[length(f)]]$fit)

  # Patient 1 is dosed at level 1 on no data. Patient 2 is dosed at the level
  # recommended by the t = 28 update -- previously they were locked to level 1
  # by an update made when nothing was known.
  expect_equal(d$dose[1], 1)
  expect_gt(d$dose[2], d$dose[1])
})

# --- min_fup_time is a stagger floor, not an accrual delay -----------------

test_that("min_fup_time does not stretch the accrual rhythm", {
  ps <- sample_no_dlt()
  a <- dosing_frame(run(0, 6, 1, ps, min_fup_time = 0)[[7]]$fit)$time
  b <- dosing_frame(run(0, 6, 1, ps, min_fup_time = 14)[[7]]$fit)$time

  # 14-day accrual with a 14-day floor is non-binding: both give a 14-day gap.
  expect_equal(a, b)
  expect_equal(diff(a), rep(14, length(a) - 1))
})

test_that("a binding min_fup_time delays the update without moving arrivals", {
  f <- run(t_max = 0, n = 4, cohort_size = 1, ps = sample_no_dlt(),
           min_fup_time = 28)

  # Arrivals stay on the 14-day pool rhythm, but each cohort is held until the
  # last dosed patient has 28 days of follow-up: 14, 42, 70, ...
  expect_equal(dosing_frame(f[[length(f)]]$fit)$time, c(14, 42, 70, 98))
})

# --- DTP moves dosing, never the update ------------------------------------

test_that("a DTP wait leaves the pool untouched", {
  # Same seed, same patients: the arrival stream must be identical whether or
  # not DTP waits. Only which arrivals enrol may differ.
  seen <- list()
  spy <- function(tag) {
    force(tag)
    function(df) {
      seen[[tag]] <<- c(seen[[tag]], 14)
      data.frame(time_delta = 14)
    }
  }
  ps <- sample_dlt_at(2)
  tp <- rep(0.30, 5)

  phase1_dtp_tite_sim(design_at(0, 8), tp, patient_sample = ps,
                      sample_patient_arrivals = spy("base"), max_time = obswin,
                      min_fup_time = 0, next_dose = 1, return_all_fits = TRUE)
  phase1_dtp_tite_sim(design_at(35, 8), tp, patient_sample = ps,
                      sample_patient_arrivals = spy("dtp"), max_time = obswin,
                      min_fup_time = 0, next_dose = 1, return_all_fits = TRUE)

  # The stream is a renewal process with a fixed gap; a wait must not reset it.
  expect_true(length(seen$dtp) >= length(seen$base))
  expect_equal(unique(seen$base), 14)
  expect_equal(unique(seen$dtp), 14)
})

test_that("no two trajectory entries share a timestamp", {
  # The waiting-room cap makes same-instant refits structurally impossible:
  # while a wait runs the room is full, so the patient who opens the next
  # cohort arrives strictly after wait_end.
  f <- run(t_max = 35, n = 8, cohort_size = 1, ps = sample_dlt_at(2),
           tp = rep(0.30, 5))
  expect_equal(anyDuplicated(update_times(f)), 0L)
})

test_that("patients held through a wait are dosed at its end, the rest on arrival", {
  f <- run(t_max = 35, n = 12, cohort_size = 3, ps = sample_dlt_at(2),
           tp = rep(0.30, 5))
  d <- dosing_frame(f[[length(f)]]$fit)
  ev <- f[[length(f)]]$dtp_wait_events

  expect_gt(nrow(ev), 0)
  w <- ev[1, ]
  held <- d$time[abs(d$time - w$time_out) < 1e-9]
  expect_gte(length(held), 1)
  # A cohort spanning a wait does not share one dosing time.
  cohort_of_held <- d$cohort[which(abs(d$time - w$time_out) < 1e-9)[1]]
  expect_gt(length(unique(d$time[d$cohort == cohort_of_held])), 1)
})

test_that("a mid-wait DLT ends the wait early and de-escalates the waiting cohort", {
  # Patient 2 is dosed at t = 28 and their DLT surfaces at t = 49.125, inside a
  # wait opened at t = 42 and planned to run to 54. Regression for #30 under
  # the new loop: the cohort that was waiting is dosed at the de-escalated
  # dose, not at the projection.
  f <- run(t_max = 35, n = 8, cohort_size = 1, ps = sample_dlt_at(2),
           tp = rep(0.30, 5))
  ev <- f[[length(f)]]$dtp_wait_events

  early <- ev[ev$ended_by == "dtp_decision", ]
  expect_equal(nrow(early), 1)
  expect_lt(early$time_out, early$time_in + 14)
  expect_lt(early$dose_after, early$dose_before)

  d <- dosing_frame(f[[length(f)]]$fit)
  expect_equal(d$dose[3], early$dose_after)
})

# --- the waiting room and its cost ----------------------------------------

test_that("n_held counts the opener, so a cohort-1 wait reports 1 (not 0)", {
  # Regression: n_held used to subtract the opener, so it was always 0 at
  # cohort size 1 even though exactly one patient was held through the wait.
  f <- run(t_max = 35, n = 8, cohort_size = 1, ps = sample_dlt_at(2),
           tp = rep(0.30, 5))
  ev <- f[[length(f)]]$dtp_wait_events
  dosing <- ev[ev$ended_by != "stopped", ]
  expect_gt(nrow(dosing), 0)
  expect_true(all(dosing$n_held == 1))
})

test_that("patients are missed exactly when the wait outruns the room", {
  # A patient is missed iff more arrive during the wait than the room can hold.
  f1 <- run(t_max = 35, n = 8, cohort_size = 1, ps = sample_dlt_at(2),
            tp = rep(0.30, 5))
  ev1 <- f1[[length(f1)]]$dtp_wait_events
  expect_true("n_missed" %in% names(ev1))
  # held-during-wait = n_held - 1 (the opener is always held); everything else
  # that arrived within the wait window is missed.
  expect_equal(ev1$n_missed,
               pmax(0, floor(ev1$wait_duration / 14) - (ev1$n_held - 1)))

  # Cohort 3 has a room three times as large, so short waits cost nothing.
  f3 <- run(t_max = 35, n = 12, cohort_size = 3, ps = sample_dlt_at(2),
            tp = rep(0.30, 5))
  ev3 <- f3[[length(f3)]]$dtp_wait_events
  expect_true(all(ev3$wait_duration[ev3$n_missed == 0] <= 3 * 14))
})

test_that("n_held and n_missed split the arrivals during a wait", {
  f <- run(t_max = 35, n = 8, cohort_size = 1, ps = sample_dlt_at(2),
           tp = rep(0.30, 5))
  ev <- f[[length(f)]]$dtp_wait_events

  # Everyone who presents during a wait is either held or turned away. n_held
  # also counts the opener (arrived before the wait), so subtract it here.
  expect_equal((ev$n_held - 1) + ev$n_missed, floor(ev$wait_duration / 14))
})

test_that("the empty wait-event log carries the n_missed column", {
  f <- run(t_max = 0, n = 6, cohort_size = 1, ps = sample_no_dlt())
  ev <- f[[length(f)]]$dtp_wait_events

  expect_equal(nrow(ev), 0)
  expect_true(all(c("n_held", "n_missed") %in% names(ev)))
})
