# Baseline property tests for escalation's TITE-CRM
#
# These tests verify properties that dtptite depends on. They do NOT re-test
# what escalation already covers (see escalation's test_dfcrm.R,
# test_simulate_trials.R, test_patient_sample.R, test_linear_follow_up_weight.R).
#
# Test groups:
# 1. TITE-CRM = CRM equivalence (all weights = 1)
# 2. Weight effects on dose recommendation
# 3. Edge cases relevant to DTP
# 4. simulate_compare() with shared PatientSample
# 5. time_to_tox_func and DLT timing
# 6. Custom weight function hook

library(escalation)

# -- Shared setup -------------------------------------------------------
skeleton <- c(0.124, 0.25, 0.398, 0.542, 0.666)
target <- 0.25

# ===== Test Group 1: TITE-CRM = CRM Equivalence ========================

test_that("TITE-CRM with all weights=1 matches CRM recommended_dose", {
  tite_model <- get_dfcrm_tite(skeleton = skeleton, target = target)
  crm_model  <- get_dfcrm(skeleton = skeleton, target = target)

  # Outcomes with all weights = 1 (fully observed)
  tite_outcomes <- data.frame(
    dose   = c(1, 1, 1, 2, 2, 2),
    tox    = c(0, 0, 0, 0, 1, 0),
    weight = c(1, 1, 1, 1, 1, 1),
    cohort = 1:6
  )
  crm_outcomes <- "1NNN 2NTN"

  tite_fit <- tite_model %>% fit(tite_outcomes)
  crm_fit  <- crm_model %>% fit(crm_outcomes)

  expect_equal(recommended_dose(tite_fit), recommended_dose(crm_fit))
})

test_that("TITE-CRM with all weights=1 matches CRM mean_prob_tox", {
  tite_model <- get_dfcrm_tite(skeleton = skeleton, target = target)
  crm_model  <- get_dfcrm(skeleton = skeleton, target = target)

  tite_outcomes <- data.frame(
    dose   = c(1, 1, 1, 2, 2, 2),
    tox    = c(0, 0, 0, 0, 1, 0),
    weight = c(1, 1, 1, 1, 1, 1),
    cohort = 1:6
  )
  crm_outcomes <- "1NNN 2NTN"

  tite_fit <- tite_model %>% fit(tite_outcomes)
  crm_fit  <- crm_model %>% fit(crm_outcomes)

  expect_equal(
    mean_prob_tox(tite_fit),
    mean_prob_tox(crm_fit),
    tolerance = 1e-6
  )
})

# ===== Test Group 2: Weight Effects on Dose Recommendation ==============

test_that("mean_prob_tox decreases monotonically as non-tox weight increases", {
  model <- get_dfcrm_tite(skeleton = skeleton, target = target)

  # 3 complete patients at dose 2, no tox; 1 pending patient at dose 2
  base_data <- data.frame(
    dose   = c(2, 2, 2, 2),
    tox    = c(0, 0, 0, 0),
    weight = c(1, 1, 1, NA),  # last weight varies
    cohort = 1:4
  )

  weights_to_test <- c(0.25, 0.5, 0.75, 1.0)
  prob_at_dose2 <- numeric(length(weights_to_test))

  for (i in seq_along(weights_to_test)) {
    base_data$weight[4] <- weights_to_test[i]
    fit_obj <- model %>% fit(base_data)
    prob_at_dose2[i] <- mean_prob_tox(fit_obj)[3]  # dose 2 is index 3 (NoDose + 5 doses)
  }

  # More weight with no tox → lower estimated toxicity: monotonically decreasing
  for (i in seq_along(prob_at_dose2)[-1]) {
    expect_lt(prob_at_dose2[i], prob_at_dose2[i - 1])
  }
})

test_that("increasing non-tox weights can cause dose escalation", {
  model <- get_dfcrm_tite(skeleton = skeleton, target = target)

  # Start with some data where pending patients hold back escalation
  # 3 complete patients at dose 1 (no tox), 3 pending patients at dose 2
  outcomes_low_weight <- data.frame(
    dose   = c(1, 1, 1, 2, 2, 2),
    tox    = c(0, 0, 0, 0, 0, 0),
    weight = c(1, 1, 1, 0.1, 0.1, 0.1),
    cohort = 1:6
  )

  outcomes_high_weight <- data.frame(
    dose   = c(1, 1, 1, 2, 2, 2),
    tox    = c(0, 0, 0, 0, 0, 0),
    weight = c(1, 1, 1, 1, 1, 1),
    cohort = 1:6
  )

  fit_low  <- model %>% fit(outcomes_low_weight)
  fit_high <- model %>% fit(outcomes_high_weight)

  # With more follow-up (higher weight) and no tox, recommended dose should

  # be >= the low-weight recommendation
  expect_gte(recommended_dose(fit_high), recommended_dose(fit_low))
})

test_that("partial vs complete weights produce different mean_prob_tox", {
  model <- get_dfcrm_tite(skeleton = skeleton, target = target)

  outcomes_partial <- data.frame(
    dose   = c(1, 1, 1, 2, 2, 2),
    tox    = c(0, 0, 0, 0, 1, 0),
    weight = c(1, 1, 1, 0.5, 1, 0.5),
    cohort = 1:6
  )

  outcomes_complete <- data.frame(
    dose   = c(1, 1, 1, 2, 2, 2),
    tox    = c(0, 0, 0, 0, 1, 0),
    weight = c(1, 1, 1, 1, 1, 1),
    cohort = 1:6
  )

  fit_partial  <- model %>% fit(outcomes_partial)
  fit_complete <- model %>% fit(outcomes_complete)

  expect_false(
    all(abs(mean_prob_tox(fit_partial) - mean_prob_tox(fit_complete)) < 1e-10)
  )
})

# ===== Test Group 3: Edge Cases Relevant to DTP =========================

test_that("TITE-CRM fits with all patients pending (no tox)", {
  model <- get_dfcrm_tite(skeleton = skeleton, target = target)

  outcomes <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(0, 0, 0),
    weight = c(0.3, 0.5, 0.7),
    cohort = 1:3
  )

  fit_obj <- model %>% fit(outcomes)
  expect_true(is.numeric(recommended_dose(fit_obj)))
  expect_true(!is.na(recommended_dose(fit_obj)))
})

test_that("TITE-CRM fits with mixed complete tox and pending patients", {
  model <- get_dfcrm_tite(skeleton = skeleton, target = target)

  outcomes <- data.frame(
    dose   = c(1, 1, 2, 2, 2),
    tox    = c(0, 1, 0, 0, 0),
    weight = c(1, 1, 0.5, 0.3, 0.8),
    cohort = 1:5
  )

  fit_obj <- model %>% fit(outcomes)
  expect_true(is.numeric(recommended_dose(fit_obj)))
  expect_true(!is.na(recommended_dose(fit_obj)))
  expect_true(all(is.numeric(mean_prob_tox(fit_obj))))
})

test_that("base CRM has no built-in stopping (continue() returns TRUE)", {
  model <- get_dfcrm_tite(skeleton = skeleton, target = target)

  outcomes <- data.frame(
    dose   = c(1, 1, 1, 1, 1, 1),
    tox    = c(1, 1, 1, 0, 0, 0),
    weight = c(1, 1, 1, 1, 1, 1),
    cohort = 1:6
  )

  fit_obj <- model %>% fit(outcomes)
  expect_true(continue(fit_obj))
})

# ===== Test Group 4: simulate_compare() with Shared Patients ============

test_that("simulate_compare with identical TITE-CRM designs gives same results", {
  true_prob_tox <- c(0.05, 0.15, 0.25, 0.35, 0.60)

  design <- get_dfcrm_tite(skeleton = skeleton, target = target) %>%
    stop_at_n(n = 9)

  designs <- list(
    "Design A" = design,
    "Design B" = design
  )

  set.seed(42)
  result <- simulate_compare(
    designs,
    num_sims = 5,
    true_prob_tox = true_prob_tox,
    max_time = 56
  )

  # Both designs should give identical recommended dose distributions
  expect_equal(
    prob_recommend(result$`Design A`),
    prob_recommend(result$`Design B`)
  )
})

test_that("simulate_compare returns simulations_collection with per-design access", {
  true_prob_tox <- c(0.05, 0.15, 0.25, 0.35, 0.60)

  design <- get_dfcrm_tite(skeleton = skeleton, target = target) %>%
    stop_at_n(n = 9)

  designs <- list(
    "Design A" = design,
    "Design B" = design
  )

  set.seed(42)
  result <- simulate_compare(
    designs,
    num_sims = 3,
    true_prob_tox = true_prob_tox,
    max_time = 56
  )

  # Per-design sims objects should exist and have standard accessors
  sims_a <- result$`Design A`
  expect_true(!is.null(sims_a))
  expect_true(is.numeric(prob_recommend(sims_a)))
  expect_true(is.numeric(trial_duration(sims_a)))
  expect_equal(length(trial_duration(sims_a)), 3)
})

# ===== Test Group 5: time_to_tox_func and DLT Timing ====================

test_that("default PatientSample: DLT times are in (0,1)", {
  set.seed(123)
  ps <- PatientSample$new(num_patients = 100)

  # tox_time should be in (0, 1) from default runif(1)
  expect_true(all(ps$tox_time >= 0))
  expect_true(all(ps$tox_time <= 1))
})

test_that("custom time_to_tox_func produces DLT times in (0, obswin)", {
  obswin <- 56
  set.seed(123)
  ps <- PatientSample$new(
    num_patients = 100,
    time_to_tox_func = function() runif(1, 0, obswin)
  )

  expect_true(all(ps$tox_time >= 0))
  expect_true(all(ps$tox_time <= obswin))
  # Should NOT all be < 1 (unlike default)
  expect_true(any(ps$tox_time > 1))
})

test_that("get_patient_tox respects time parameter for DLT observation", {
  set.seed(42)
  # Create patient with known tox_time via set_eff_and_tox
  ps <- PatientSample$new()
  ps$set_eff_and_tox(
    tox_u = c(0.05),      # low tox_u means tox will occur at high prob_tox
    eff_u = c(0.5),
    tox_time = c(28),      # DLT manifests at day 28
    eff_time = c(0)
  )

  prob_tox <- 0.5  # tox_u=0.05 < 0.5, so tox occurs

  # At time < 28, DLT not yet observed
  expect_equal(ps$get_patient_tox(1, prob_tox, time = 14), 0)

  # At time >= 28, DLT is observed
  expect_equal(ps$get_patient_tox(1, prob_tox, time = 28), 1)
  expect_equal(ps$get_patient_tox(1, prob_tox, time = 56), 1)
})

test_that("custom time_to_tox_func changes simulation behaviour", {
  true_prob_tox <- c(0.05, 0.15, 0.25, 0.35, 0.60)

  design <- get_dfcrm_tite(skeleton = skeleton, target = target) %>%
    stop_at_n(n = 12)

  num_sims <- 30

  # Create patient samples with SAME tox_u (who gets tox) but DIFFERENT

  # tox_time (when tox manifests). We do this by manually setting tox_u/tox_time.
  set.seed(100)
  ps_narrow <- lapply(1:num_sims, function(i) {
    ps <- PatientSample$new(num_patients = 50)
    # Override: tox times in (0, 1) — DLTs always visible at first interim
    ps_tox_u <- ps$tox_u
    ps_eff_u <- ps$eff_u
    ps$set_eff_and_tox(
      tox_u = ps_tox_u, eff_u = ps_eff_u,
      tox_time = runif(50, 0, 1), eff_time = rep(0, 50)
    )
    ps
  })

  set.seed(100)
  ps_wide <- lapply(1:num_sims, function(i) {
    ps <- PatientSample$new(num_patients = 50)
    # Override: tox times in (0, 56) — DLTs may be hidden at interims
    ps_tox_u <- ps$tox_u
    ps_eff_u <- ps$eff_u
    ps$set_eff_and_tox(
      tox_u = ps_tox_u, eff_u = ps_eff_u,
      tox_time = runif(50, 0, 56), eff_time = rep(0, 50)
    )
    ps
  })

  # Same seed for simulation mechanics (arrival times etc.)
  set.seed(200)
  sims_narrow <- simulate_compare(
    list("design" = design),
    num_sims = num_sims,
    true_prob_tox = true_prob_tox,
    max_time = 56,
    patient_samples = ps_narrow
  )

  set.seed(200)
  sims_wide <- simulate_compare(
    list("design" = design),
    num_sims = num_sims,
    true_prob_tox = true_prob_tox,
    max_time = 56,
    patient_samples = ps_wide
  )

  # With narrow timing, DLTs are observed immediately → different dose decisions
  # than wide timing where DLTs may be hidden during interim analyses
  expect_false(identical(
    prob_recommend(sims_narrow$design),
    prob_recommend(sims_wide$design)
  ))
})

# ===== Test Group 6: Custom Weight Function Hook ========================

test_that("custom get_weight function works in simulate_trials", {
  true_prob_tox <- c(0.05, 0.15, 0.25, 0.35, 0.60)

  design <- get_dfcrm_tite(skeleton = skeleton, target = target) %>%
    stop_at_n(n = 9)

  # Step weight: 0 for first 14 days, then jumps to 1
  step_weight <- function(now_time, recruited_time, tox, max_time,
                          tox_has_weight_1 = TRUE) {
    follow_up <- now_time - recruited_time
    w <- ifelse(follow_up >= 14, 1, 0)
    if (tox_has_weight_1) {
      w <- ifelse(tox == 1, 1, w)
    }
    w
  }

  set.seed(200)
  sims_step <- design %>%
    simulate_trials(
      num_sims = 20,
      true_prob_tox = true_prob_tox,
      max_time = 56,
      get_weight = step_weight
    )

  set.seed(200)
  sims_linear <- design %>%
    simulate_trials(
      num_sims = 20,
      true_prob_tox = true_prob_tox,
      max_time = 56
    )

  # Different weight functions should produce different recommendation probs
  expect_false(
    identical(prob_recommend(sims_step), prob_recommend(sims_linear))
  )
})

test_that("get_dfcrm_tite with pre-computed weights matches dfcrm::titecrm", {
  # Verify our entry point faithfully wraps dfcrm with sarcoma trial parameters
  outcomes <- data.frame(
    dose   = c(1, 2, 2, 2),
    tox    = c(0, 0, 1, 0),
    weight = c(1, 0.5, 1, 0.75),
    cohort = 1:4
  )

  esc_fit <- get_dfcrm_tite(
    skeleton = skeleton, target = target, scale = sqrt(1.34)
  ) %>% fit(outcomes)

  dfcrm_fit <- dfcrm::titecrm(
    prior   = skeleton,
    target  = target,
    tox     = c(0, 0, 1, 0),
    level   = c(1, 2, 2, 2),
    weights = c(1, 0.5, 1, 0.75),
    model   = "empiric",
    method  = "bayes",
    scale   = sqrt(1.34)
  )

  expect_equal(recommended_dose(esc_fit), dfcrm_fit$mtd)
  expect_equal(mean_prob_tox(esc_fit), dfcrm_fit$ptox, tolerance = 1e-6)
})

# ===== Test Group 7: Credible Intervals ==================================

test_that("TITE-CRM credible intervals match underlying dfcrm fit", {
  model <- get_dfcrm_tite(skeleton = skeleton, target = target)

  outcomes <- data.frame(
    dose   = c(1, 1, 2, 2, 2),
    tox    = c(0, 0, 0, 1, 0),
    weight = c(1, 1, 0.5, 1, 0.75),
    cohort = 1:5
  )

  fit_obj <- model %>% fit(outcomes)

  # Access credible intervals from underlying dfcrm fit (default conf.level = 0.9)
  dfcrm_lower <- fit_obj$dfcrm_fit$ptoxL
  dfcrm_upper <- fit_obj$dfcrm_fit$ptoxU

  # escalation's prob_tox_quantile at matching quantiles (5th/95th for 90% CI)
  lower <- prob_tox_quantile(fit_obj, p = 0.05)
  upper <- prob_tox_quantile(fit_obj, p = 0.95)

  expect_equal(lower, dfcrm_lower, tolerance = 1e-6)
  expect_equal(upper, dfcrm_upper, tolerance = 1e-6)

  # Lower < mean < upper for all doses
  mean_pt <- mean_prob_tox(fit_obj)
  expect_true(all(lower < mean_pt))
  expect_true(all(mean_pt < upper))

  # prob_tox_exceeds returns probabilities in [0, 1]
  exceed <- prob_tox_exceeds(fit_obj, threshold = target)
  expect_true(is.numeric(exceed))
  expect_length(exceed, length(skeleton))
  expect_true(all(exceed >= 0 & exceed <= 1))
})
