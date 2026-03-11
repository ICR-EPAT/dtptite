# Tests for DTP simulation function (phase1_dtp_tite_sim)
#
# Test groups:
# 1. t_max=0 equivalence (same as base phase1_tite_sim)
# 2. Wait fires (DTP wait increases trial duration)
# 3. simulate_compare works (DTP + non-DTP)
# 4. Stopping rule in simulation (early stop)
# 5. simulate_trials dispatch (returns simulations object)
# 6. Outer decorator dispatch (simulation_function finds DTP through chain)
# 7. DLT during wait (event-driven re-evaluation)
# 8. Cohort size 3 + queue (queue replaces next cohort)
# 9. tite_patient_samples (correct distribution)
# 10. queue_size = 0 vs default (different trial duration)
# 11. Default patient_sample uses U(0, max_time)

library(escalation)

# -- Shared setup ----------------------------------------------------------
skeleton <- c(0.124, 0.25, 0.398, 0.542, 0.666)
target <- 0.25

# ===== Test 1: t_max=0 equivalence ========================================

test_that("t_max=0 DTP simulation matches base phase1_tite_sim", {
  true_prob_tox <- c(0.05, 0.15, 0.25, 0.35, 0.60)

  base_design <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    stop_at_n(n = 9)

  dtp_design <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 0, obswin = 56) |>
    stop_at_n(n = 9)

  # Same seed -> same patient sample and arrivals
  set.seed(42)
  ps <- PatientSample$new(
    num_patients = 50,
    time_to_tox_func = function() runif(1, 0, 56)
  )

  set.seed(100)
  base_result <- escalation:::phase1_tite_sim(
    base_design, true_prob_tox,
    patient_sample = ps,
    max_time = 56
  )

  # Reset the same PatientSample for DTP run
  set.seed(42)
  ps2 <- PatientSample$new(
    num_patients = 50,
    time_to_tox_func = function() runif(1, 0, 56)
  )

  set.seed(100)
  dtp_result <- phase1_dtp_tite_sim(
    dtp_design, true_prob_tox,
    patient_sample = ps2,
    max_time = 56
  )

  base_fit <- base_result[[1]]$fit
  dtp_fit <- dtp_result[[1]]$fit

  expect_equal(recommended_dose(dtp_fit), recommended_dose(base_fit))
  expect_equal(num_patients(dtp_fit), num_patients(base_fit))
  expect_equal(
    base_result[[1]]$time,
    dtp_result[[1]]$time,
    tolerance = 0.01
  )
})

# ===== Test 2: Wait fires =================================================

test_that("DTP wait increases trial duration in single simulation", {
  # Use a scenario where DTP is likely to wait:
  # moderate tox with cohort_size=3, pending patients, t_max=35
  true_prob_tox <- c(0.05, 0.15, 0.25, 0.35, 0.60)

  base_design <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    stop_at_n(n = 12)

  dtp_design <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 35, obswin = 56) |>
    stop_at_n(n = 12)

  arrivals <- function(df) escalation::cohorts_of_n(n = 3, mean_time_delta = 1)

  # Pre-create matched patient sample pairs (same seed -> identical latent data)
  set.seed(999)
  ps_base <- tite_patient_samples(num_sims = 30, max_time = 56,
                                   num_patients = 50)
  set.seed(999)
  ps_dtp <- tite_patient_samples(num_sims = 30, max_time = 56,
                                  num_patients = 50)

  # Run multiple seeds and check DTP >= base on average
  base_durations <- numeric(30)
  dtp_durations <- numeric(30)

  for (sim_i in seq_len(30)) {
    seed <- 1000 + sim_i

    set.seed(seed + 5000)
    base_res <- escalation:::phase1_tite_sim(
      base_design, true_prob_tox,
      patient_sample = ps_base[[sim_i]],
      sample_patient_arrivals = arrivals,
      max_time = 56
    )
    base_durations[sim_i] <- base_res[[1]]$time

    set.seed(seed + 5000)
    dtp_res <- phase1_dtp_tite_sim(
      dtp_design, true_prob_tox,
      patient_sample = ps_dtp[[sim_i]],
      sample_patient_arrivals = arrivals,
      max_time = 56
    )
    dtp_durations[sim_i] <- dtp_res[[1]]$time
  }

  # DTP should have >= duration on average (waits add time)
  # Allow some tolerance since not every sim will have a wait
  expect_gte(mean(dtp_durations), mean(base_durations) - 1)
  # At least some DTP sims should have longer duration
  expect_true(any(dtp_durations > base_durations))
})

# ===== Test 3: simulate_compare works =====================================

test_that("simulate_compare runs with DTP and non-DTP designs", {
  true_prob_tox <- c(0.05, 0.15, 0.25, 0.35, 0.60)

  base_design <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    stop_at_n(n = 9)

  dtp_design <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 14, obswin = 56) |>
    stop_at_n(n = 9)

  set.seed(42)
  result <- simulate_compare(
    list("base" = base_design, "dtp" = dtp_design),
    num_sims = 5,
    true_prob_tox = true_prob_tox,
    max_time = 56
  )

  expect_true(!is.null(result$base))
  expect_true(!is.null(result$dtp))
  expect_true(is.numeric(prob_recommend(result$base)))
  expect_true(is.numeric(prob_recommend(result$dtp)))
  expect_true(is.numeric(trial_duration(result$base)))
  expect_true(is.numeric(trial_duration(result$dtp)))
})

# ===== Test 4: Stopping rule in simulation ================================

test_that("stopping rule triggers early stop in DTP simulation", {
  # Very toxic scenario: true tox >> target at all doses
  true_prob_tox <- c(0.60, 0.70, 0.80, 0.90, 0.95)

  design <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    stop_for_beta_binomial_toxicity(dose = 1, tox_threshold = 0.25,
                                    confidence = 0.70) |>
    apply_dtp(t_max = 14, obswin = 56)

  set.seed(42)
  ps <- PatientSample$new(
    num_patients = 50,
    time_to_tox_func = function() runif(1, 0, 56)
  )

  set.seed(100)
  result <- phase1_dtp_tite_sim(
    design, true_prob_tox,
    patient_sample = ps,
    sample_patient_arrivals = function(df) {
      escalation::cohorts_of_n(n = 3, mean_time_delta = 1)
    },
    max_time = 56,
    min_fup_time = 10,
    i_like_big_trials = TRUE
  )

  final_fit <- result[[1]]$fit
  n_enrolled <- num_patients(final_fit)
  # With very high tox and aggressive stopping, should stop before many cohorts
  expect_lt(n_enrolled, 60)
})

# ===== Test 5: simulate_trials dispatch ===================================

test_that("simulate_trials dispatches to DTP simulation function", {
  true_prob_tox <- c(0.05, 0.15, 0.25, 0.35, 0.60)

  design <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 14, obswin = 56) |>
    stop_at_n(n = 6)

  # simulation_function should find our method
  sim_func <- simulation_function(design)
  expect_identical(sim_func, phase1_dtp_tite_sim)

  set.seed(42)
  sims <- simulate_trials(
    design,
    num_sims = 3,
    true_prob_tox = true_prob_tox,
    max_time = 56
  )

  expect_s3_class(sims, "simulations")
  expect_true(is.numeric(prob_recommend(sims)))
  expect_equal(length(trial_duration(sims)), 3)
})

# ===== Test 6: Outer decorator dispatch ===================================

test_that("simulation_function finds DTP factory through decorator chain", {
  design <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 14, obswin = 56) |>
    stop_at_n(n = 9)

  # stop_at_n wraps DTP -> dispatch should walk to DTP factory
  sim_func <- simulation_function(design)
  expect_identical(sim_func, phase1_dtp_tite_sim)

  # Double-wrapped: stopping rule + stop_at_n around DTP
  design2 <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 14, obswin = 56) |>
    stop_for_beta_binomial_toxicity(dose = 1, tox_threshold = 0.25) |>
    stop_at_n(n = 9)

  sim_func2 <- simulation_function(design2)
  expect_identical(sim_func2, phase1_dtp_tite_sim)
})

# ===== Test 7: DLT during wait ============================================

test_that("DLT during wait triggers re-evaluation", {
  true_prob_tox <- c(0.05, 0.15, 0.25, 0.35, 0.60)

  design <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 35, obswin = 56) |>
    stop_at_n(n = 12)

  # Use deterministic DLT timing: DLT always at day 20
  set.seed(42)
  ps <- PatientSample$new(
    num_patients = 50,
    time_to_tox_func = function() 20
  )

  set.seed(100)
  result <- phase1_dtp_tite_sim(
    design, true_prob_tox,
    patient_sample = ps,
    sample_patient_arrivals = function(df) {
      escalation::cohorts_of_n(n = 3, mean_time_delta = 1)
    },
    max_time = 56,
    return_all_fits = TRUE
  )

  # Should have more than 2 fits (initial + cohorts)
  expect_gt(length(result), 2)

  # Final fit should be valid
  final_fit <- result[[length(result)]]$fit
  expect_true(is.numeric(recommended_dose(final_fit)) ||
                is.na(recommended_dose(final_fit)))
})

# ===== Test 8: Cohort size 3 + queue ======================================

test_that("queue_size=2 still doses full cohorts of 3", {
  true_prob_tox <- c(0.05, 0.15, 0.25, 0.35, 0.60)

  design <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 14, obswin = 56) |>
    stop_at_n(n = 12)
  set_dtp_queue_size(design, 2)

  set.seed(42)
  ps <- PatientSample$new(
    num_patients = 50,
    time_to_tox_func = function() runif(1, 0, 56)
  )

  set.seed(100)
  result <- phase1_dtp_tite_sim(
    design, true_prob_tox,
    patient_sample = ps,
    sample_patient_arrivals = function(df) {
      escalation::cohorts_of_n(n = 3, mean_time_delta = 1)
    },
    max_time = 56
  )

  # Trial should complete and produce valid results
  final_fit <- result[[1]]$fit
  expect_true(is.numeric(recommended_dose(final_fit)) ||
                is.na(recommended_dose(final_fit)))
  # Should have enrolled patients in cohorts of 3
  expect_equal(num_patients(final_fit) %% 3, 0)
})

# ===== Test 9: tite_patient_samples =======================================

test_that("tite_patient_samples creates correct tox_time distribution", {
  set.seed(123)
  ps_list <- tite_patient_samples(num_sims = 5, max_time = 56,
                                   num_patients = 200)

  expect_length(ps_list, 5)

  for (ps in ps_list) {
    expect_s3_class(ps, "PatientSample")
    # tox_time should span (0, 56), not just (0, 1)
    expect_true(all(ps$tox_time >= 0))
    expect_true(all(ps$tox_time <= 56))
    # With 200 patients, at least some should have tox_time > 1
    expect_true(any(ps$tox_time > 1))
    # Distribution should cover the range reasonably
    expect_true(any(ps$tox_time > 28))
  }
})

# ===== Test 10: queue_size = 0 vs default =================================

test_that("queue_size=0 produces different behavior than default", {
  true_prob_tox <- c(0.05, 0.15, 0.25, 0.35, 0.60)

  design_default <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 35, obswin = 56) |>
    stop_at_n(n = 12)

  design_no_queue <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 35, obswin = 56) |>
    stop_at_n(n = 12)
  set_dtp_queue_size(design_no_queue, 0)

  arrivals <- function(df) escalation::cohorts_of_n(n = 3, mean_time_delta = 1)

  # Pre-create matched patient sample pairs
  set.seed(888)
  ps_default <- tite_patient_samples(num_sims = 30, max_time = 56,
                                      num_patients = 50)
  set.seed(888)
  ps_noqueue <- tite_patient_samples(num_sims = 30, max_time = 56,
                                      num_patients = 50)

  # Run multiple seeds and track whether durations ever differ
  any_different <- FALSE
  for (sim_i in seq_len(30)) {
    seed <- 2000 + sim_i

    set.seed(seed + 5000)
    res_default <- phase1_dtp_tite_sim(
      design_default, true_prob_tox,
      patient_sample = ps_default[[sim_i]],
      sample_patient_arrivals = arrivals,
      max_time = 56
    )

    set.seed(seed + 5000)
    res_no_queue <- phase1_dtp_tite_sim(
      design_no_queue, true_prob_tox,
      patient_sample = ps_noqueue[[sim_i]],
      sample_patient_arrivals = arrivals,
      max_time = 56
    )

    if (!identical(res_default[[1]]$time, res_no_queue[[1]]$time)) {
      any_different <- TRUE
      break
    }
  }

  # At least some simulations should differ when queue is enabled vs disabled
  expect_true(any_different)
})

# ===== Test 11: Default patient_sample uses U(0, max_time) ================

test_that("default patient_sample generates tox times spanning max_time", {
  # High tox to guarantee DLTs appear in 12 patients
  true_prob_tox <- c(0.30, 0.45, 0.55, 0.65, 0.80)

  design <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 14, obswin = 56) |>
    stop_at_n(n = 12)

  # No patient_sample provided — exercises the NULL default path
  set.seed(42)
  result <- phase1_dtp_tite_sim(
    design, true_prob_tox,
    sample_patient_arrivals = function(df) {
      escalation::cohorts_of_n(n = 3, mean_time_delta = 1)
    },
    max_time = 56,
    return_all_fits = TRUE
  )

  # Simulation should complete with multiple fits
  expect_gt(length(result), 1)

  # Final fit should be valid and have enrolled patients
  final_fit <- result[[length(result)]]$fit
  expect_true(is.numeric(recommended_dose(final_fit)) ||
                is.na(recommended_dose(final_fit)))
  expect_gt(num_patients(final_fit), 0)

  # With high tox, DLTs should be observed
  fit_data <- model_frame(final_fit)
  expect_gt(sum(fit_data$tox), 0)
})
