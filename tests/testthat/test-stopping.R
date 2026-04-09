# Beta-Binomial conjugate stopping rule tests

library(escalation)

# -- Shared setup -------------------------------------------------------
skeleton <- c(0.124, 0.25, 0.398, 0.542, 0.666)
target <- 0.25

# Helper: create a fitted model with stopping rule
fit_with_stop <- function(outcomes, dose = 1, tox_threshold = target,
                          confidence = 0.88, a = NULL, b = NULL) {
  model <- get_dfcrm_tite(skeleton = skeleton, target = target) %>%
    stop_for_beta_binomial_toxicity(
      dose = dose, tox_threshold = tox_threshold,
      confidence = confidence, a = a, b = b
    )
  model %>% fit(outcomes)
}

# ===== Test 1: 3/3 DLTs at dose 1 triggers stop ==========================

test_that("3/3 DLTs at dose 1 triggers stop", {
  outcomes <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(1, 1, 1),
    weight = c(1, 1, 1),
    cohort = 1:3
  )
  fit_obj <- fit_with_stop(outcomes)

  expect_false(continue(fit_obj))
  expect_true(is.na(recommended_dose(fit_obj)))
})

# ===== Test 2: 2/3 DLTs does NOT trigger stop ============================

test_that("1/3 DLTs at dose 1 does not trigger stop", {
  outcomes <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(1, 0, 0),
    weight = c(1, 1, 1),
    cohort = 1:3
  )
  fit_obj <- fit_with_stop(outcomes)

  # 1 - pbeta(0.25, 0.25+1, 0.75+2) ≈ 0.550 < 0.88
  expect_true(continue(fit_obj))
})

# ===== Test 3: 0 patients at monitored dose — no stop ====================

test_that("No patients at monitored dose does not trigger stop", {
  outcomes <- data.frame(
    dose   = c(2, 2, 2),
    tox    = c(1, 1, 1),
    weight = c(1, 1, 1),
    cohort = 1:3
  )
  # Monitor dose 1 only (default)
  fit_obj <- fit_with_stop(outcomes, dose = 1)

  expect_true(continue(fit_obj))
})

# ===== Test 4: Exact boundary math =======================================

test_that("Exceedance probability matches hand-computed value", {
  # 3/3 DLTs at dose 1 with default prior Beta(0.25, 0.75)
  # Posterior: Beta(0.25+3, 0.75+0) = Beta(3.25, 0.75)
  # Pr(p > 0.25) = 1 - pbeta(0.25, 3.25, 0.75)
  expected <- 1 - pbeta(0.25, 3.25, 0.75)

  outcomes <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(1, 1, 1),
    weight = c(1, 1, 1),
    cohort = 1:3
  )
  fit_obj <- fit_with_stop(outcomes)
  exceedance <- dtptite:::.bb_exceedance_probs(fit_obj)

  expect_equal(exceedance[1], expected)
  expect_true(expected > 0.88)  # confirms stop
})

test_that("1/3 DLTs exceedance below confidence threshold", {
  # Posterior: Beta(0.25+1, 0.75+2) = Beta(1.25, 2.75)
  expected <- 1 - pbeta(0.25, 1.25, 2.75)

  outcomes <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(1, 0, 0),
    weight = c(1, 1, 1),
    cohort = 1:3
  )
  fit_obj <- fit_with_stop(outcomes)
  exceedance <- dtptite:::.bb_exceedance_probs(fit_obj)

  expect_equal(exceedance[1], expected)
  expect_true(expected < 0.88)  # confirms no stop
})

# ===== Test 5: Default prior from tox_target() ===========================

test_that("Default prior resolves from tox_target()", {
  outcomes <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(1, 1, 1),
    weight = c(1, 1, 1),
    cohort = 1:3
  )
  # Default: a=NULL, b=NULL -> Beta(0.25, 0.75)
  fit_obj <- fit_with_stop(outcomes)

  # Verify resolved prior
  expect_equal(fit_obj$a, 0.25)
  expect_equal(fit_obj$b, 0.75)
  expect_false(continue(fit_obj))
})

# ===== Test 6: Custom prior Beta(1,1) ====================================

test_that("Custom prior Beta(1,1) gives different boundary", {
  # 1/1 DLTs at dose 1
  outcomes <- data.frame(
    dose   = c(1),
    tox    = c(1),
    weight = c(1),
    cohort = 1
  )

  # Default prior Beta(0.25, 0.75): exceedance ≈ 0.868 < 0.88 -> no stop
  fit_default <- fit_with_stop(outcomes)
  expect_true(continue(fit_default))

  # Uniform prior Beta(1, 1): exceedance = 1 - pbeta(0.25, 2, 1) ≈ 0.9375 > 0.88
  fit_uniform <- fit_with_stop(outcomes, a = 1, b = 1)
  expect_false(continue(fit_uniform))
})

# ===== Test 7: Monitors specific dose only ================================

test_that("Monitoring specific dose ignores other doses", {
  outcomes <- data.frame(
    dose   = c(1, 1, 1, 2, 2, 2),
    tox    = c(0, 0, 0, 1, 1, 1),
    weight = c(1, 1, 1, 1, 1, 1),
    cohort = 1:6
  )

  # Monitor dose 1 only — no DLTs there, no stop
  fit_d1 <- fit_with_stop(outcomes, dose = 1)
  expect_true(continue(fit_d1))

  # Monitor dose 2 — all DLTs, should stop
  fit_d2 <- fit_with_stop(outcomes, dose = 2)
  expect_false(continue(fit_d2))
})

# ===== Test 8: dose = "recommended" ======================================

test_that("dose='recommended' checks the recommended dose", {
  # 3/3 DLTs at dose 1, model still recommends dose 1 (only dose with data)
  outcomes <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(1, 1, 1),
    weight = c(1, 1, 1),
    cohort = 1:3
  )
  fit_obj <- fit_with_stop(outcomes, dose = "recommended")

  # The parent's recommended_dose should be 1 (lowest dose, all toxic)
  parent_rec <- recommended_dose(fit_obj$parent)
  expect_equal(parent_rec, 1)

  # With recommended dose being dose 1 and 3/3 DLTs there: stop
  expect_false(continue(fit_obj))
})

# ===== Test 9: dose = "any" ==============================================

test_that("dose='any' eliminates doses but continues if admissible doses remain", {
  outcomes <- data.frame(
    dose   = c(1, 1, 1, 3, 3, 3),
    tox    = c(0, 0, 0, 1, 1, 1),
    weight = c(1, 1, 1, 1, 1, 1),
    cohort = 1:6
  )

  fit_obj <- fit_with_stop(outcomes, dose = "any")
  # Dose 3 eliminated (3/3 DLTs) but dose 1 is safe -> trial continues
  expect_true(continue(fit_obj))
  expect_false(dose_admissible(fit_obj)[3])
  # Recommended dose clamps below eliminated dose
  expect_lte(recommended_dose(fit_obj), 2)
})

test_that("dose='any' stops when all reachable doses are eliminated", {
  outcomes <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(1, 1, 1),
    weight = c(1, 1, 1),
    cohort = 1:3
  )

  fit_obj <- fit_with_stop(outcomes, dose = "any")
  # Dose 1 eliminated (3/3 DLTs), no lower dose -> trial stops
  expect_false(continue(fit_obj))
  expect_true(is.na(recommended_dose(fit_obj)))
})

# ===== Test 10: dose = vector ============================================

test_that("dose=c(1,2) stops when either monitored dose exceeds", {
  outcomes <- data.frame(
    dose   = c(1, 1, 1, 2, 2, 2),
    tox    = c(0, 0, 0, 1, 1, 1),
    weight = c(1, 1, 1, 1, 1, 1),
    cohort = 1:6
  )

  # Monitor doses 1 and 2; dose 2 has 3/3 DLTs -> stop
  fit_obj <- fit_with_stop(outcomes, dose = c(1, 2))
  expect_false(continue(fit_obj))

  # Monitor dose 1 only; dose 1 has 0/3 DLTs -> no stop
  fit_d1 <- fit_with_stop(outcomes, dose = 1)
  expect_true(continue(fit_d1))
})

# ===== Test 10b: Invalid string dose errors ===============================

test_that("String dose level like '1' is rejected", {
  model <- get_dfcrm_tite(skeleton = skeleton, target = target)
  expect_error(
    stop_for_beta_binomial_toxicity(
      model, dose = "1", tox_threshold = target, confidence = 0.88
    ),
    "\"recommended\", or \"any\""
  )
})

# ===== Test 10d: parent stops first =======================================

test_that("continue=FALSE when parent stops but beta-binom would not", {
  outcomes <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(0, 0, 0),
    weight = c(1, 1, 1),
    cohort = 1:3
  )

  # Parent stop_at_n(3) fires; beta-binom has 0 DLTs so would not stop
  model <- get_dfcrm_tite(skeleton = skeleton, target = target) %>%
    stop_at_n(n = 3) %>%
    stop_for_beta_binomial_toxicity(
      dose = 1, tox_threshold = target, confidence = 0.88
    )
  fit_obj <- model %>% fit(outcomes)

  expect_false(continue(fit_obj))
})

# ===== Test 10e: recommended_dose delegates when not stopping =============

test_that("recommended_dose delegates to parent when not stopping", {
  outcomes <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(0, 0, 0),
    weight = c(1, 1, 1),
    cohort = 1:3
  )
  fit_obj <- fit_with_stop(outcomes)

  # No DLTs -> no stop -> recommended_dose should match parent
  expect_equal(recommended_dose(fit_obj), recommended_dose(fit_obj$parent))
  expect_false(is.na(recommended_dose(fit_obj)))
})

# ===== Test 11: Chains with DTP ==========================================

test_that("Stopping rule chains with DTP decorator", {
  outcomes <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(1, 1, 1),
    weight = c(1, 1, 1),
    cohort = 1:3
  )

  model <- get_dfcrm_tite(skeleton = skeleton, target = target) %>%
    stop_for_beta_binomial_toxicity(
      dose = 1, tox_threshold = target, confidence = 0.88
    ) %>%
    apply_dtp(t_max = 56, obswin = 56)

  fit_obj <- model %>% fit(outcomes)

  # Stopping takes effect even through DTP
  expect_false(continue(fit_obj))
  expect_true(is.na(recommended_dose(fit_obj)))

  # DTP generics still accessible
  expect_false(dtp_should_wait(fit_obj))
})

# ===== Test 12: Chains with stop_at_n ====================================

test_that("Stopping rule chains with stop_at_n", {
  outcomes <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(0, 0, 0),
    weight = c(1, 1, 1),
    cohort = 1:3
  )

  model <- get_dfcrm_tite(skeleton = skeleton, target = target) %>%
    stop_for_beta_binomial_toxicity(
      dose = 1, tox_threshold = target, confidence = 0.88
    ) %>%
    stop_at_n(n = 3)

  fit_obj <- model %>% fit(outcomes)

  # stop_at_n should trigger (3 patients = n)
  expect_false(continue(fit_obj))
})

# ===== Test 13: dose_admissible marks toxic doses =========================

test_that("dose_admissible cascades elimination to higher doses", {
  # Dose 1 eliminated -> all doses above also eliminated (BOIN rule)
  outcomes_d1 <- data.frame(
    dose   = c(1, 1, 1, 2, 2, 2),
    tox    = c(1, 1, 1, 0, 0, 0),
    weight = c(1, 1, 1, 1, 1, 1),
    cohort = 1:6
  )

  fit_obj <- fit_with_stop(outcomes_d1, dose = "any")
  admissible <- dose_admissible(fit_obj)

  expect_false(admissible[1])
  # Cascade: dose 2+ also eliminated even though 0/3 DLTs
  expect_false(admissible[2])

  # Dose 3 eliminated -> doses 4-5 also eliminated, doses 1-2 safe
  outcomes_d3 <- data.frame(
    dose   = c(1,1,1, 2,2,2, 3,3,3),
    tox    = c(0,0,0, 0,0,0, 1,1,1),
    weight = rep(1, 9),
    cohort = 1:9
  )

  fit_d3 <- fit_with_stop(outcomes_d3, dose = "any")
  admissible_d3 <- dose_admissible(fit_d3)

  expect_true(admissible_d3[1])
  expect_true(admissible_d3[2])
  expect_false(admissible_d3[3])
  expect_false(admissible_d3[4])  # cascade
  expect_false(admissible_d3[5])  # cascade
})
