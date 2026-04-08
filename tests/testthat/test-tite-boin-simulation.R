# Tests for TITE-BOIN simulation integration (Step 7)
#
# Mirrors DTP-TITE-CRM simulation tests (test-dtp-simulation.R) using
# TITE-BOIN as the base selector, plus BOIN-specific tests.
#
# Test groups:
# 1.  t_max=0 equivalence (same as base phase1_tite_sim)
# 2.  Wait fires (DTP wait increases trial duration)
# 3.  simulate_compare works (DTP + non-DTP BOIN)
# 4.  Stopping rule in simulation (early stop)
# 5.  simulate_trials dispatch (returns simulations object)
# 6.  Outer decorator dispatch (simulation_function finds DTP through chain)
# 7.  DLT during wait (event-driven re-evaluation)
# 8.  Cohort size 3 + queue (queue replaces next cohort)
# 9.  queue_size = 0 vs default (different trial duration)
# 10. TITE-BOIN standalone simulation (no DTP)
# 11. select_boin_mtd in simulation
# 12. All-toxic scenario with TITE-BOIN
# 13. Elimination + DTP in simulation
# 14. Elimination without DTP in simulation
# 15. Cross-design comparison (all 4 designs)
# 16. Operating characteristics match TITEgBOIN reference

library(escalation)

# -- Shared setup ----------------------------------------------------------
target <- 0.25

# ===== Test 1: t_max=0 equivalence ========================================

test_that("t_max=0 DTP-BOIN simulation matches base phase1_tite_sim", {
  true_prob_tox <- c(0.05, 0.15, 0.25, 0.35, 0.60)

  base_design <- get_boin_tite(5, target) |>
    stop_at_n(n = 9)

  dtp_design <- get_boin_tite(5, target) |>
    apply_dtp(t_max = 0, obswin = 56) |>
    stop_at_n(n = 9)

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

test_that("DTP wait increases trial duration for TITE-BOIN", {
  # BOIN's discrete lambda decisions need: (1) moderate tox so phat lands
  # in the "stay" zone, and (2) enough pending patients with low weights
  # so that waiting could push phat below lambda_e.
  # Use single-patient accrual (cohort_size=1) with fast arrivals so
  # multiple pending patients accumulate at the same dose.
  true_prob_tox <- c(0.15, 0.25, 0.35, 0.50, 0.65)

  base_design <- get_boin_tite(5, target) |>
    stop_at_n(n = 24)

  dtp_design <- get_boin_tite(5, target) |>
    apply_dtp(t_max = 35, obswin = 56) |>
    stop_at_n(n = 24)

  arrivals <- function(df) data.frame(time_delta = 7)

  set.seed(999)
  ps_base <- tite_patient_samples(num_sims = 30, max_time = 56,
                                   num_patients = 100)
  set.seed(999)
  ps_dtp <- tite_patient_samples(num_sims = 30, max_time = 56,
                                  num_patients = 100)

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

  expect_gte(mean(dtp_durations), mean(base_durations) - 1)
  expect_true(any(dtp_durations > base_durations))
})

# ===== Test 3: simulate_compare works =====================================

test_that("simulate_compare runs with DTP and non-DTP BOIN designs", {
  true_prob_tox <- c(0.05, 0.15, 0.25, 0.35, 0.60)

  base_design <- get_boin_tite(5, target) |>
    stop_at_n(n = 9)

  dtp_design <- get_boin_tite(5, target) |>
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

test_that("stopping rule triggers early stop in DTP-BOIN simulation", {
  true_prob_tox <- c(0.60, 0.70, 0.80, 0.90, 0.95)

  design <- get_boin_tite(5, target) |>
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
  # With 60%+ tox at all doses, stopping should fire before exhausting
  # the 50 available patients
  expect_lt(n_enrolled, 50)
})

# ===== Test 5: simulate_trials dispatch ===================================

test_that("simulate_trials dispatches to DTP simulation for BOIN", {
  true_prob_tox <- c(0.05, 0.15, 0.25, 0.35, 0.60)

  design <- get_boin_tite(5, target) |>
    apply_dtp(t_max = 14, obswin = 56) |>
    stop_at_n(n = 6)

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

test_that("simulation_function finds DTP factory through BOIN decorator chain", {
  design <- get_boin_tite(5, target) |>
    apply_dtp(t_max = 14, obswin = 56) |>
    stop_at_n(n = 9)

  sim_func <- simulation_function(design)
  expect_identical(sim_func, phase1_dtp_tite_sim)

  design2 <- get_boin_tite(5, target) |>
    apply_dtp(t_max = 14, obswin = 56) |>
    stop_for_beta_binomial_toxicity(dose = 1, tox_threshold = 0.25) |>
    stop_at_n(n = 9)

  sim_func2 <- simulation_function(design2)
  expect_identical(sim_func2, phase1_dtp_tite_sim)
})

# ===== Test 7: DLT during wait ============================================

test_that("DLT during wait triggers re-evaluation for BOIN", {
  true_prob_tox <- c(0.05, 0.15, 0.25, 0.35, 0.60)

  design <- get_boin_tite(5, target) |>
    apply_dtp(t_max = 35, obswin = 56) |>
    stop_at_n(n = 12)

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

  expect_gt(length(result), 2)

  final_fit <- result[[length(result)]]$fit
  expect_true(is.numeric(recommended_dose(final_fit)) ||
                is.na(recommended_dose(final_fit)))
})

# ===== Test 8: Cohort size 3 + queue ======================================

test_that("queue_size=2 still doses full cohorts of 3 for BOIN", {
  true_prob_tox <- c(0.05, 0.15, 0.25, 0.35, 0.60)

  design <- get_boin_tite(5, target) |>
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

  final_fit <- result[[1]]$fit
  expect_true(is.numeric(recommended_dose(final_fit)) ||
                is.na(recommended_dose(final_fit)))
  expect_equal(num_patients(final_fit) %% 3, 0)
})

# ===== Test 9: queue_size = 0 vs default ==================================

test_that("queue_size=0 produces different behavior than default for BOIN", {
  true_prob_tox <- c(0.15, 0.25, 0.35, 0.50, 0.65)

  design_default <- get_boin_tite(5, target) |>
    apply_dtp(t_max = 35, obswin = 56) |>
    stop_at_n(n = 24)

  design_no_queue <- get_boin_tite(5, target) |>
    apply_dtp(t_max = 35, obswin = 56) |>
    stop_at_n(n = 24)
  set_dtp_queue_size(design_no_queue, 0)

  arrivals <- function(df) data.frame(time_delta = 7)

  set.seed(888)
  ps_default <- tite_patient_samples(num_sims = 30, max_time = 56,
                                      num_patients = 50)
  set.seed(888)
  ps_noqueue <- tite_patient_samples(num_sims = 30, max_time = 56,
                                      num_patients = 50)

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

  expect_true(any_different)
})

# ===== Test 10: TITE-BOIN standalone simulation ===========================

test_that("TITE-BOIN standalone simulation works via simulate_trials", {
  true_prob_tox <- c(0.05, 0.15, 0.25, 0.35, 0.60)

  design <- get_boin_tite(5, target) |>
    stop_at_n(n = 12)

  sim_func <- simulation_function(design)
  expect_identical(sim_func, escalation:::phase1_tite_sim)

  set.seed(42)
  sims <- simulate_trials(
    design,
    num_sims = 5,
    true_prob_tox = true_prob_tox,
    max_time = 56
  )

  expect_s3_class(sims, "simulations")
  expect_length(prob_recommend(sims), 5 + 1)  # num_doses + "no dose"
  expect_true(all(trial_duration(sims) > 0))
  expect_true(all(num_patients(sims) > 0))
  expect_true(all(num_patients(sims) <= 12))
})

# ===== Test 11: select_boin_mtd in simulation =============================

test_that("select_boin_mtd works in TITE-BOIN simulation", {
  true_prob_tox <- c(0.05, 0.15, 0.25, 0.35, 0.60)

  design <- get_boin_tite(5, target) |>
    stop_at_n(n = 12) |>
    select_boin_mtd()

  set.seed(42)
  sims <- simulate_trials(
    design,
    num_sims = 5,
    true_prob_tox = true_prob_tox,
    max_time = 56
  )

  pr <- prob_recommend(sims)
  expect_length(pr, 5 + 1)
  expect_equal(sum(pr), 1, tolerance = 0.01)
})

# ===== Test 12: All-toxic scenario ========================================

test_that("all-toxic scenario triggers early stop and recommends no dose", {
  true_prob_tox <- c(0.60, 0.70, 0.80, 0.90, 0.95)

  design <- get_boin_tite(5, target) |>
    stop_for_beta_binomial_toxicity(dose = 1, tox_threshold = 0.25,
                                    confidence = 0.95) |>
    stop_at_n(n = 24) |>
    select_boin_mtd()

  set.seed(42)
  sims <- simulate_trials(
    design,
    num_sims = 10,
    true_prob_tox = true_prob_tox,
    max_time = 56
  )

  pr <- prob_recommend(sims)
  # "NoDose" is the first element
  expect_gt(pr["NoDose"], 0.3)
  # Early stopping should reduce mean patients below max
  expect_lt(mean(num_patients(sims)), 24)
})

# ===== Test 13: Elimination + DTP in simulation ===========================

test_that("elimination works correctly in DTP-BOIN simulation", {
  # Dose 1 safe, doses 2+ very toxic -> elimination should fire for dose 2+
  true_prob_tox <- c(0.10, 0.50, 0.60, 0.70, 0.80)

  design <- get_boin_tite(5, target) |>
    stop_for_beta_binomial_toxicity(dose = 1, tox_threshold = target,
                                    confidence = 0.95, a = 1, b = 1) |>
    apply_dtp(t_max = 35, obswin = 56) |>
    stop_at_n(n = 24)

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
    return_all_fits = TRUE
  )

  # Simulation should complete
  expect_gt(length(result), 1)

  # Check that elimination fired at some point: at least one intermediate fit
  # should have dose_admissible() with FALSE entries
  any_eliminated <- FALSE
  for (i in seq_along(result)) {
    fit_i <- result[[i]]$fit
    da <- dose_admissible(fit_i)
    if (any(!da)) {
      any_eliminated <- TRUE
      # Eliminated dose should not be recommended
      rec <- recommended_dose(fit_i)
      if (!is.na(rec)) {
        expect_true(da[rec],
          info = paste("Recommended dose", rec, "is not admissible at fit", i))
      }
      break
    }
  }
  expect_true(any_eliminated,
    info = "Expected at least one fit with eliminated doses")
})

# ===== Test 14: Elimination without DTP in simulation =====================

test_that("elimination works in standalone TITE-BOIN simulation", {
  true_prob_tox <- c(0.10, 0.50, 0.60, 0.70, 0.80)

  design <- get_boin_tite(5, target) |>
    stop_for_beta_binomial_toxicity(dose = 1, tox_threshold = target,
                                    confidence = 0.95, a = 1, b = 1) |>
    stop_at_n(n = 24) |>
    select_boin_mtd()

  set.seed(42)
  sims <- simulate_trials(
    design,
    num_sims = 10,
    true_prob_tox = true_prob_tox,
    max_time = 56
  )

  pr <- prob_recommend(sims)
  # Dose 1 is safe (0.10), doses 2-5 are very toxic -> should concentrate
  # on dose 1 or no-dose
  expect_gt(pr["NoDose"] + pr["1"], 0.5)
})

# ===== Test 15: Cross-design comparison (all 4 designs) ===================

test_that("simulate_compare works with all 4 designs", {
  skip_if_not_installed("dfcrm")

  true_prob_tox <- c(0.05, 0.15, 0.25, 0.35, 0.60)
  skeleton <- c(0.124, 0.25, 0.398, 0.542, 0.666)

  titecrm <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    stop_at_n(n = 9)

  dtp_titecrm <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 14, obswin = 56) |>
    stop_at_n(n = 9)

  boin <- get_boin_tite(5, target) |>
    stop_at_n(n = 9)

  dtp_boin <- get_boin_tite(5, target) |>
    apply_dtp(t_max = 14, obswin = 56) |>
    stop_at_n(n = 9)

  designs <- list(
    "TITE-CRM"     = titecrm,
    "DTP-TITE-CRM" = dtp_titecrm,
    "TITE-BOIN"    = boin,
    "DTP-TITE-BOIN" = dtp_boin
  )

  set.seed(42)
  result <- simulate_compare(
    designs,
    num_sims = 5,
    true_prob_tox = true_prob_tox,
    max_time = 56
  )

  for (name in names(designs)) {
    expect_true(!is.null(result[[name]]),
      info = paste("Missing results for", name))
    expect_true(is.numeric(prob_recommend(result[[name]])),
      info = paste("prob_recommend not numeric for", name))
    expect_true(is.numeric(trial_duration(result[[name]])),
      info = paste("trial_duration not numeric for", name))
  }
})

# ===== Test 16: Operating characteristics match TITEgBOIN ==================

test_that("TITE-BOIN selection probabilities match TITEgBOIN reference", {
  skip_if_not_installed("TITEgBOIN")

  # -- Shared parameters -----------------------------------------------------
  true_prob_tox <- c(0.02, 0.08, 0.12, 0.25, 0.38)
  tgt <- 0.25
  num_doses <- 5
  n_patients <- 24
  obswin <- 56
  accrual_interval <- 14  # 1 patient every 14 days
  p.saf <- 0.6 * tgt
  p.tox <- 1.4 * tgt
  cutoff.eli <- 0.88
  n_sims <- 1000

  # -- TITEgBOIN reference ---------------------------------------------------
  # accrual = patients per observation window = obswin / accrual_interval
  tg_result <- TITEgBOIN::get_oc_TITE_QuasiBOIN(
    target     = tgt,
    prob       = true_prob_tox,
    score      = NA,
    TITE       = TRUE,
    ncohort    = n_patients,
    cohortsize = 1,
    maxt       = obswin,
    accrual    = obswin / accrual_interval,
    maxpen     = 0.5,
    alpha1     = 0.5,     # uniform DLT timing within window
    alpha2     = 0.5,
    n.earlystop = 100,
    Neli       = 3,
    startdose  = 1,
    p.saf      = p.saf,
    p.tox      = p.tox,
    cutoff.eli = cutoff.eli,
    ntrial     = n_sims,
    seed       = 42
  )

  # Extract selection percentages (list element $selpercent, in %)
  tg_sel_prop <- tg_result$selpercent / 100

  # -- dtptite simulation ----------------------------------------------------
  design <- get_boin_tite(num_doses, tgt, p.saf = p.saf, p.tox = p.tox) |>
    stop_for_beta_binomial_toxicity(
      dose = 1, tox_threshold = tgt,
      confidence = cutoff.eli, a = 1, b = 1
    ) |>
    dont_skip_doses() |>
    stop_at_n(n = n_patients) |>
    select_boin_mtd()

  # Explicitly create patient samples with U(0, obswin) DLT timing.
  # Use simulate_compare (not simulate_trials) because it correctly
  # iterates patient_samples[[i]] per sim replicate.
  set.seed(42)
  ps_list <- tite_patient_samples(
    num_sims = n_sims, max_time = obswin, num_patients = 100
  )

  set.seed(42)
  result <- simulate_compare(
    list("boin" = design),
    num_sims = n_sims,
    true_prob_tox = true_prob_tox,
    max_time = obswin,
    sample_patient_arrivals = function(df) data.frame(time_delta = accrual_interval),
    patient_samples = ps_list
  )

  our_pr <- prob_recommend(result$boin)
  our_no_dose <- our_pr["NoDose"]
  our_sel <- as.numeric(our_pr[as.character(1:num_doses)])

  # -- Comparison -------------------------------------------------------------
  # Dose with highest selection % must match
  tg_best <- which.max(tg_sel_prop)
  our_best <- which.max(our_sel)
  expect_equal(our_best, tg_best,
    info = sprintf(
      "Best dose mismatch: dtptite=%d (%.1f%%), TITEgBOIN=%d (%.1f%%)",
      our_best, our_sel[our_best] * 100,
      tg_best, tg_sel_prop[tg_best] * 100
    )
  )

  # All selection probabilities within 10% absolute difference
  abs_diff <- abs(our_sel - tg_sel_prop)
  for (d in seq_len(num_doses)) {
    expect_lt(abs_diff[d], 0.10,
      label = sprintf(
        "Dose %d diff (dtptite=%.1f%%, TITEgBOIN=%.1f%%, diff=%.1f%%)",
        d, our_sel[d] * 100, tg_sel_prop[d] * 100, abs_diff[d] * 100
      )
    )
  }
})
