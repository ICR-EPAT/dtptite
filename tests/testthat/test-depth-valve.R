# Regression: the depth safety valve must not truncate wait-heavy DTP trials.
#
# `phase1_dtp_tite_sim` records a trajectory entry per cohort AND per DTP wait,
# but the safety valve counts trajectory entries (`i < max_i`, max_i = 30). A
# wait-heavy cohort-size-1 trial therefore used to hit the cap and terminate
# with a valid recommendation, one or more patients short of its planned sample.
# The simulator now defaults `i_like_big_trials = TRUE` so the valve is off;
# termination is guaranteed by the design's own stopping rules.

library(escalation)

test_that("a wait-heavy cohort-1 DTP trial enrols its full planned sample", {
  skeleton <- c(0.140, 0.250, 0.376, 0.502, 0.615)
  n_max <- 24
  design <- get_dfcrm_tite(skeleton = skeleton, target = 0.25,
                           model = "empiric", scale = 0.759) |>
    stop_for_beta_binomial_toxicity(a = 1, b = 1, dose = 1,
                                    tox_threshold = 0.25, confidence = 0.88) |>
    stop_at_n(n = n_max) |>
    dont_skip_doses(when_escalating = TRUE, when_deescalating = TRUE) |>
    apply_dtp(t_max = 35, obswin = 56)

  arrivals <- function(df) data.frame(time_delta = 14)
  true_prob_tox <- c(0.25, 0.38, 0.50, 0.62, 0.72)   # lowest dose at target

  # Seed 94 produces 7 DTP waits at cohort size 1 — enough to push the trajectory
  # index past max_i = 30 before all 24 patients enrol.
  set.seed(94)
  ps <- PatientSample$new(num_patients = n_max,
                          time_to_tox_func = function() runif(1, 0, 56))
  res <- phase1_dtp_tite_sim(
    design, true_prob_tox, patient_sample = ps,
    sample_patient_arrivals = arrivals, max_time = 56,
    return_all_fits = TRUE
  )                                     # default i_like_big_trials = TRUE
  last <- res[[length(res)]]

  # The trial waits repeatedly and still enrols its full sample.
  expect_gte(nrow(last$dtp_wait_events), 6L)
  expect_equal(num_patients(last$fit), n_max)

  # With the valve re-enabled the same trial is truncated early and warns — this
  # is the bug the default guards against.
  set.seed(94)
  ps2 <- PatientSample$new(num_patients = n_max,
                           time_to_tox_func = function() runif(1, 0, 56))
  expect_warning(
    truncated <- phase1_dtp_tite_sim(
      design, true_prob_tox, patient_sample = ps2,
      sample_patient_arrivals = arrivals, max_time = 56,
      return_all_fits = TRUE, i_like_big_trials = FALSE
    ),
    "max depth reached"
  )
  expect_lt(num_patients(truncated[[length(truncated)]]$fit), n_max)
})
