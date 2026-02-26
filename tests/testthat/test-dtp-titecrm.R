# DTP decorator tests on TITE-CRM
#
# Tests the 7 branches of the DTP decision flow:
# 1. Escalate immediately (no wait needed)
# 2. Stay at current dose, no wait benefit
# 3. De-escalate, no wait benefit
# 4. Wait leads to dose increase
# 5. Wait too long (t_max too short), proceed without waiting
# 6. Stop trial (continue() == FALSE)
# 7. t_max=0 equivalence (DTP has no effect)
# 8. t_max=Inf (unlimited waiting, effective bound)

library(escalation)

# -- Shared setup -------------------------------------------------------
skeleton <- c(0.124, 0.25, 0.398, 0.542, 0.666)
target <- 0.25

# ===== Test 1: Escalate immediately ======================================

test_that("DTP does not wait when model already recommends escalation", {
  # All patients complete at dose 1, no DLTs -> model recommends higher dose
  outcomes <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(0, 0, 0),
    weight = c(1, 1, 1),
    cohort = 1:3
  )

  dtp_model <- get_dfcrm_tite(skeleton = skeleton, target = target) %>%
    apply_dtp(t_max = 56, obswin = 56)

  fit_obj <- dtp_model %>% fit(outcomes)

  # Model recommends escalation (rec > current dose of 1)
  expect_gt(recommended_dose(fit_obj), 1)
  # No waiting needed

  expect_false(dtp_should_wait(fit_obj))
  expect_equal(dtp_wait_time(fit_obj), 0)
})

# ===== Test 2: Stay, no wait benefit ====================================

test_that("DTP does not wait when all complete and model says stay", {

  # 3 complete at D2 with 1 DLT + 1 complete at D2 no tox -> model says stay at D2
  # (all weights = 1, so no pending patients to project)
  outcomes <- data.frame(
    dose   = c(2, 2, 2, 2),
    tox    = c(0, 1, 0, 0),
    weight = c(1, 1, 1, 1),
    cohort = 1:4
  )

  dtp_model <- get_dfcrm_tite(skeleton = skeleton, target = target) %>%
    apply_dtp(t_max = 56, obswin = 56)

  fit_obj <- dtp_model %>% fit(outcomes)

  # Model recommends dose 2 (stay at current)
  expect_equal(recommended_dose(fit_obj), 2)
  # No waiting (all complete, projection engine finds no pending patients)
  expect_false(dtp_should_wait(fit_obj))
})

# ===== Test 3: De-escalate, no wait benefit ==============================

test_that("DTP does not wait when model says de-escalate and all complete", {
  # Heavy toxicity at dose 2 -> model de-escalates to dose 1
  outcomes <- data.frame(
    dose   = c(2, 2, 2),
    tox    = c(1, 1, 0),
    weight = c(1, 1, 1),
    cohort = 1:3
  )

  dtp_model <- get_dfcrm_tite(skeleton = skeleton, target = target) %>%
    apply_dtp(t_max = 56, obswin = 56)

  fit_obj <- dtp_model %>% fit(outcomes)

  # Model de-escalates
  expect_equal(recommended_dose(fit_obj), 1)
  # No waiting (all complete, no projection benefit)
  expect_false(dtp_should_wait(fit_obj))
})

# ===== Test 4: Wait leads to dose increase ===============================
#
# Data: 3@D2(w=1, 1 DLT) + 1@D2(w=0.3, no tox)
# At w=0.3: rec_dose=1 (de-escalate from current dose 2)
# At w≈0.51: rec_dose flips to 2
# With obswin=56: flip occurs at step 12 (w=0.3 + 12/56 ≈ 0.514)

test_that("DTP recommends waiting when projection shows dose increase", {
  outcomes <- data.frame(
    dose   = c(2, 2, 2, 2),
    tox    = c(0, 1, 0, 0),
    weight = c(1, 1, 1, 0.3),
    cohort = 1:4
  )

  dtp_model <- get_dfcrm_tite(skeleton = skeleton, target = target) %>%
    apply_dtp(t_max = 56, obswin = 56)

  fit_obj <- dtp_model %>% fit(outcomes)

  # Current recommendation is dose 1 (de-escalate)
  expect_equal(recommended_dose(fit_obj), 1)
  # DTP says wait
  expect_true(dtp_should_wait(fit_obj))
  # Wait time is positive
  expect_gt(dtp_wait_time(fit_obj), 0)
  # Projected dose is higher than current recommendation
  expect_gt(dtp_projected_dose(fit_obj), recommended_dose(fit_obj))
  # Weight-space reporting also works
  expect_gt(dtp_wait_time(fit_obj, type = "weight"), 0)
})

# ===== Test 5: Wait too long, proceed ====================================
#
# Same data as test 4 but t_max is too short for the flip to occur.
# Flip needs ~12 steps; use t_max = 5.

test_that("DTP does not wait when t_max is too short for dose increase", {
  outcomes <- data.frame(
    dose   = c(2, 2, 2, 2),
    tox    = c(0, 1, 0, 0),
    weight = c(1, 1, 1, 0.3),
    cohort = 1:4
  )

  dtp_model <- get_dfcrm_tite(skeleton = skeleton, target = target) %>%
    apply_dtp(t_max = 5, obswin = 56)

  fit_obj <- dtp_model %>% fit(outcomes)

  # Current recommendation is still dose 1
  expect_equal(recommended_dose(fit_obj), 1)
  # DTP says don't wait (t_max too short)
  expect_false(dtp_should_wait(fit_obj))
})

# ===== Test 6: Stop trial ================================================

test_that("DTP respects trial stopping from parent decorator", {
  # stop_at_n(n = 6) with 6 patients -> continue() == FALSE
  outcomes <- data.frame(
    dose   = c(1, 1, 1, 1, 1, 1),
    tox    = c(1, 1, 0, 0, 0, 0),
    weight = c(1, 1, 1, 1, 1, 1),
    cohort = 1:6
  )

  dtp_model <- get_dfcrm_tite(skeleton = skeleton, target = target) %>%
    stop_at_n(n = 6) %>%
    apply_dtp(t_max = 56, obswin = 56)

  fit_obj <- dtp_model %>% fit(outcomes)

  # Trial should stop
  expect_false(continue(fit_obj))
  # recommended_dose still returns a value (model recommendation)
  expect_true(!is.na(recommended_dose(fit_obj)))
})

# ===== Test 7: t_max=0 equivalence =======================================

test_that("apply_dtp with t_max=0 is equivalent to unwrapped model", {
  outcomes <- data.frame(
    dose   = c(2, 2, 2, 2),
    tox    = c(0, 1, 0, 0),
    weight = c(1, 1, 1, 0.3),
    cohort = 1:4
  )

  base_model <- get_dfcrm_tite(skeleton = skeleton, target = target)
  dtp_model  <- base_model %>% apply_dtp(t_max = 0, obswin = 56)

  fit_base <- base_model %>% fit(outcomes)
  fit_dtp  <- dtp_model  %>% fit(outcomes)

  expect_equal(recommended_dose(fit_dtp), recommended_dose(fit_base))
  expect_equal(mean_prob_tox(fit_dtp), mean_prob_tox(fit_base))
  expect_false(dtp_should_wait(fit_dtp))
})

# ===== Test 8: t_max=Inf (unlimited waiting) =============================
#
# Same data as test 4. With t_max=Inf, the effective bound is
# (1 - 0.3) * 56 = 39.2 -> effective max = 40 steps.
# The flip still occurs at step 12.

test_that("DTP with t_max=Inf terminates via effective bound", {
  outcomes <- data.frame(
    dose   = c(2, 2, 2, 2),
    tox    = c(0, 1, 0, 0),
    weight = c(1, 1, 1, 0.3),
    cohort = 1:4
  )

  dtp_model <- get_dfcrm_tite(skeleton = skeleton, target = target) %>%
    apply_dtp(t_max = Inf, obswin = 56)

  fit_obj <- dtp_model %>% fit(outcomes)

  # DTP should find the same flip as test 4
  expect_true(dtp_should_wait(fit_obj))
  expect_gt(dtp_wait_time(fit_obj), 0)
  # Weight-space reporting
  wt <- dtp_wait_time(fit_obj, type = "weight")
  expect_gt(wt, 0)
  expect_lte(wt, 1)  # weight gain cannot exceed 1
})
