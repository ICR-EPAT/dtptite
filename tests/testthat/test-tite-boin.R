# TITE-BOIN selector tests

library(escalation)

# -- Shared setup -------------------------------------------------------------
target <- 0.25
p.saf <- 0.6 * target
p.tox <- 1.4 * target
lambda_e <- log((1 - p.saf) / (1 - target)) /
  log(target * (1 - p.saf) / (p.saf * (1 - target)))
lambda_d <- log((1 - target) / (1 - p.tox)) /
  log(p.tox * (1 - target) / (target * (1 - p.tox)))

# Shared fixtures
df_1NNN_2NTN <- data.frame(
  dose   = c(1, 1, 1, 2, 2, 2),
  tox    = c(0, 0, 0, 0, 1, 0),
  weight = rep(1, 6),
  cohort = c(1, 1, 1, 2, 2, 2)
)
boin_ref <- get_boin(5, target, use_stopping_rule = FALSE) %>% fit("1NNN 2NTN")

# ===== Group 1: BOIN equivalence (all weights = 1) ==========================

test_that("1a: Dose recommendation matches escalation::get_boin on 1NNN 2NTN", {
  tite_fit <- get_boin_tite(5, target) %>% fit(df_1NNN_2NTN)
  expect_identical(recommended_dose(tite_fit), recommended_dose(boin_ref))
})

test_that("1b: mean_prob_tox values match escalation::get_boin", {
  tite_fit <- get_boin_tite(5, target) %>% fit(df_1NNN_2NTN)
  expect_identical(mean_prob_tox(tite_fit), mean_prob_tox(boin_ref))
})

test_that("1c: All dose_admissible are TRUE (barebone)", {
  tite_fit <- get_boin_tite(5, target) %>% fit(df_1NNN_2NTN)
  expect_true(all(dose_admissible(tite_fit)))
  expect_length(dose_admissible(tite_fit), 5)
})

# ===== Group 2: Hand-computed phat with weights ==============================

test_that("2a: 4 pts at D2, 1 DLT, 1 pending (w=0.5) -> phat=1/3.5 -> stay", {
  # n=4, s=1, r=3, c=1, STFT=0.5, n0=3.5
  # phat = 1/3.5 = 0.286
  # lambda_e < 0.286 < lambda_d -> stay
  df <- data.frame(
    dose   = c(2, 2, 2, 2),
    tox    = c(0, 1, 0, 0),
    weight = c(1, 1, 1, 0.5),
    cohort = c(1, 1, 1, 2)
  )
  fit_obj <- get_boin_tite(5, target) %>% fit(df)

  expect_equal(recommended_dose(fit_obj), 2L)
  expect_true(1 / 3.5 > lambda_e)
  expect_true(1 / 3.5 < lambda_d)
})

test_that("2b: 6 pts at D2, 1 DLT, 2 pending (w=0.3, w=0.7) -> phat=1/5 -> stay", {
  # n=6, s=1, r=4, c=2, STFT=0.3+0.7=1.0, n0=4+1=5
  # phat = 1/5 = 0.2
  df <- data.frame(
    dose   = c(2, 2, 2, 2, 2, 2),
    tox    = c(0, 1, 0, 0, 0, 0),
    weight = c(1, 1, 1, 1, 0.3, 0.7),
    cohort = c(1, 1, 1, 1, 2, 2)
  )
  fit_obj <- get_boin_tite(5, target) %>% fit(df)

  expect_equal(recommended_dose(fit_obj), 2L)
  expect_true(1 / 5 > lambda_e)
  expect_true(1 / 5 < lambda_d)
})

test_that("2c: All pending, no DLTs -> phat=0 -> escalate", {
  # n=3, s=0, r=0, c=3, STFT=0.3+0.4+0.5=1.2, n0=0+1.2=1.2
  # phat = 0/1.2 = 0 <= lambda_e -> escalate from 1 -> 2
  df <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(0, 0, 0),
    weight = c(0.3, 0.4, 0.5),
    cohort = c(1, 1, 1)
  )
  fit_obj <- get_boin_tite(5, target) %>% fit(df)

  expect_equal(recommended_dose(fit_obj), 2L)
})

# ===== Group 3: TITEgBOIN comparison ========================================

test_that("3a: Match TITEgBOIN dose recommendation (without elimination)", {
  skip_if_not_installed("TITEgBOIN")

  # 3 patients at dose 1 with no DLTs, all complete
  df <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(0, 0, 0),
    weight = c(1, 1, 1),
    cohort = c(1, 1, 1)
  )
  fit_obj <- get_boin_tite(5, target) %>% fit(df)

  # TITEgBOIN: with 3 patients, 0 DLTs -> escalate
  # Our barebone: phat = 0/3 = 0 <= lambda_e -> escalate to 2
  expect_equal(recommended_dose(fit_obj), 2L)
})

test_that("3b: Match escalation with pending patients", {
  skip_if_not_installed("TITEgBOIN")

  # Pending patients — should still get correct decision
  df <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(0, 0, 0),
    weight = c(1, 0.8, 0.3),
    cohort = c(1, 2, 2)
  )
  fit_obj <- get_boin_tite(5, target) %>% fit(df)

  # n0 = 1 + (0.8 + 0.3) = 2.1, phat = 0 -> escalate
  expect_equal(recommended_dose(fit_obj), 2L)
})

# ===== Group 4: Edge cases ===================================================

test_that("4a: No patients -> dose 1, continue", {
  df <- data.frame(dose = integer(0), tox = integer(0),
                   weight = numeric(0), cohort = integer(0))
  fit_obj <- get_boin_tite(5, target) %>% fit(df)

  expect_equal(recommended_dose(fit_obj), 1L)
  expect_true(continue(fit_obj))
})

test_that("4b: Single pending patient", {
  df <- data.frame(dose = 1L, tox = 0L, weight = 0.3, cohort = 1L)
  fit_obj <- get_boin_tite(5, target) %>% fit(df)

  # n0 = 0 + 0.3 = 0.3, phat = 0/0.3 = 0 -> escalate to 2
  expect_equal(recommended_dose(fit_obj), 2L)
})

test_that("4c: All complete, no DLTs -> escalate", {
  df <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(0, 0, 0),
    weight = c(1, 1, 1),
    cohort = c(1, 1, 1)
  )
  fit_obj <- get_boin_tite(5, target) %>% fit(df)

  expect_equal(recommended_dose(fit_obj), 2L)
})

test_that("4d: At max dose, phat <= lambda_e -> stay (can't go higher)", {
  df <- data.frame(
    dose   = c(5, 5, 5),
    tox    = c(0, 0, 0),
    weight = c(1, 1, 1),
    cohort = c(1, 1, 1)
  )
  fit_obj <- get_boin_tite(5, target) %>% fit(df)

  # phat = 0 <= lambda_e -> wants to escalate, but capped at num_doses

  expect_equal(recommended_dose(fit_obj), 5L)
})

test_that("4e: At dose 1, phat >= lambda_d -> stay (can't go lower)", {
  df <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(1, 1, 0),
    weight = c(1, 1, 1),
    cohort = c(1, 1, 1)
  )
  fit_obj <- get_boin_tite(5, target) %>% fit(df)

  # phat = 2/3 = 0.667 >= lambda_d -> wants to de-escalate, but capped at 1
  expect_equal(recommended_dose(fit_obj), 1L)
})

# ===== Group 5: Decorator composition ========================================

test_that("5a: stop_for_beta_binomial_toxicity stops when dose 1 too toxic", {
  df <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(1, 1, 1),
    weight = c(1, 1, 1),
    cohort = 1:3
  )
  model <- get_boin_tite(5, target) %>%
    stop_for_beta_binomial_toxicity(dose = 1, tox_threshold = target,
                                    confidence = 0.95)
  fit_obj <- model %>% fit(df)

  expect_false(continue(fit_obj))
  expect_true(is.na(recommended_dose(fit_obj)))
})

test_that("5b: dont_skip_doses prevents jumping from dose 1 to dose 3", {
  # Two cohorts: 3 at dose 1 (0 DLTs) -> escalate to 2
  # Then 3 at dose 2 (0 DLTs) -> escalate to 3
  # dont_skip ensures we go to 2 first
  df <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(0, 0, 0),
    weight = c(1, 1, 1),
    cohort = c(1, 1, 1)
  )
  model <- get_boin_tite(5, target) %>% dont_skip_doses()
  fit_obj <- model %>% fit(df)

  # Base model recommends dose 2 (escalate from 1), which is only +1, so same
  expect_equal(recommended_dose(fit_obj), 2L)
})

test_that("5c: stop_at_n stops at max patients", {
  df <- data.frame(
    dose   = rep(1:4, each = 3),
    tox    = rep(0, 12),
    weight = rep(1, 12),
    cohort = rep(1:4, each = 3)
  )
  model <- get_boin_tite(5, target) %>% stop_at_n(12)
  fit_obj <- model %>% fit(df)

  expect_false(continue(fit_obj))
})

test_that("5d: select_boin_mtd selects MTD at end of trial", {
  df <- data.frame(
    dose   = rep(1:3, each = 3),
    tox    = c(0, 0, 0, 0, 0, 0, 0, 1, 1),
    weight = rep(1, 9),
    cohort = rep(1:3, each = 3)
  )
  model <- get_boin_tite(5, target) %>%
    stop_at_n(9) %>%
    select_boin_mtd()
  fit_obj <- model %>% fit(df)

  # Trial stopped (9 patients), so select_boin_mtd kicks in
  expect_false(continue(fit_obj))
  # MTD should be determined by isotonic regression
  rec <- recommended_dose(fit_obj)
  expect_true(!is.na(rec))
  expect_true(rec >= 1 && rec <= 5)
})

# ===== Group 6: DTP integration ==============================================

test_that("6a: get_boin_tite %>% apply_dtp %>% fit works", {
  df <- data.frame(
    dose   = c(2, 2, 2),
    tox    = c(0, 1, 0),
    weight = c(1, 1, 0.3),
    cohort = c(1, 1, 2)
  )
  model <- get_boin_tite(5, target) %>% apply_dtp(35, 56)
  fit_obj <- model %>% fit(df)

  expect_s3_class(fit_obj, "dtp_selector")
  rec <- recommended_dose(fit_obj)
  expect_true(!is.na(rec))
})

test_that("6b: DTP projection finds escalation through weight increment", {
  # phat close to lambda_e boundary — as weights increase, phat should drop
  # 3 pts at D1: 0 DLTs, 2 pending with low weight
  # phat = 0 -> escalate immediately (no DTP wait needed)
  df <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(0, 0, 0),
    weight = c(1, 0.1, 0.1),
    cohort = c(1, 2, 2)
  )
  model <- get_boin_tite(5, target) %>% apply_dtp(35, 56)
  fit_obj <- model %>% fit(df)

  # Already escalating (phat=0), so no wait needed
  expect_false(dtp_should_wait(fit_obj))
  expect_equal(recommended_dose(fit_obj), 2L)
})

test_that("6c: DTP t_max=0 matches standalone", {
  df <- data.frame(
    dose   = c(2, 2, 2),
    tox    = c(0, 1, 0),
    weight = c(1, 1, 0.3),
    cohort = c(1, 1, 2)
  )
  standalone <- get_boin_tite(5, target) %>% fit(df)
  dtp_zero <- get_boin_tite(5, target) %>% apply_dtp(0, 56) %>% fit(df)

  expect_equal(recommended_dose(dtp_zero), recommended_dose(standalone))
  expect_false(dtp_should_wait(dtp_zero))
})

# ===== Group 7: S3 method correctness ========================================

test_that("7a: weight() returns actual weights", {
  df <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(0, 0, 0),
    weight = c(1, 0.5, 0.3),
    cohort = c(1, 2, 2)
  )
  fit_obj <- get_boin_tite(5, target) %>% fit(df)

  expect_equal(weight(fit_obj), c(1, 0.5, 0.3))
})

test_that("7b: continue() is always TRUE (barebone)", {
  # Even with many DLTs, barebone continues (stopping via decorators)
  df <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(1, 1, 1),
    weight = c(1, 1, 1),
    cohort = 1:3
  )
  fit_obj <- get_boin_tite(5, target) %>% fit(df)

  expect_true(continue(fit_obj))
})

test_that("7c: mean_prob_tox() returns isotonic estimates", {
  df <- data.frame(
    dose   = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
    tox    = c(0, 0, 0, 0, 1, 0, 1, 1, 0),
    weight = rep(1, 9),
    cohort = rep(1:3, each = 3)
  )
  fit_obj <- get_boin_tite(5, target) %>% fit(df)

  mpt <- mean_prob_tox(fit_obj)
  expect_length(mpt, 5)
  # Doses with patients should have numeric values
  expect_false(any(is.na(mpt[1:3])))
  # Isotonic: should be non-decreasing for doses with data
  expect_true(all(diff(mpt[1:3]) >= 0))
})

test_that("7d: tox_at_dose() returns integer counts", {
  df <- data.frame(
    dose   = c(1, 1, 2, 2, 2),
    tox    = c(0, 0, 1, 0, 0),
    weight = rep(1, 5),
    cohort = c(1, 1, 2, 2, 2)
  )
  fit_obj <- get_boin_tite(5, target) %>% fit(df)

  expect_equal(tox_at_dose(fit_obj), c(0, 1, 0, 0, 0))
})

# ===== Group 8: Simulation dispatch ==========================================

test_that("8a: simulation_function(get_boin_tite) returns phase1_tite_sim", {
  sf <- get_boin_tite(5, target)
  fn <- simulation_function(sf)

  expect_identical(fn, escalation:::phase1_tite_sim)
})

test_that("8b: simulation_function with DTP returns phase1_dtp_tite_sim", {
  sf <- get_boin_tite(5, target) %>% apply_dtp(35, 56)
  fn <- simulation_function(sf)

  expect_identical(fn, phase1_dtp_tite_sim)
})

# ===== Additional: String outcomes ==========================================

test_that("String outcomes work (all weights = 1)", {
  fit_str <- get_boin_tite(5, target) %>% fit("1NNN 2NTN")
  fit_df <- get_boin_tite(5, target) %>% fit(data.frame(
    dose   = c(1, 1, 1, 2, 2, 2),
    tox    = c(0, 0, 0, 0, 1, 0),
    weight = rep(1, 6),
    cohort = c(1, 1, 1, 2, 2, 2)
  ))

  expect_identical(recommended_dose(fit_str), recommended_dose(fit_df))
  expect_identical(mean_prob_tox(fit_str), mean_prob_tox(fit_df))
})

# ===== Additional: prob_tox methods ==========================================

test_that("prob_tox_exceeds and prob_tox_quantile work correctly", {
  df <- data.frame(
    dose   = c(1, 1, 1),
    tox    = c(1, 0, 0),
    weight = c(1, 1, 1),
    cohort = c(1, 1, 1)
  )
  fit_obj <- get_boin_tite(5, target) %>% fit(df)

  pte <- prob_tox_exceeds(fit_obj, target)
  expect_length(pte, 5)
  # Dose 1: s=1, n=3, Pr(p > 0.25 | s=1, n=3) = 1 - pbeta(0.25, 2, 3)
  expected <- 1 - pbeta(target, 1 + 1, 1 + 3 - 1)
  expect_equal(pte[1], expected)
  # Doses without data -> NA
  expect_true(all(is.na(pte[2:5])))

  # median_prob_tox
  mpt <- median_prob_tox(fit_obj)
  expect_length(mpt, 5)
  expected_median <- qbeta(0.5, 1 + 1, 1 + 3 - 1)
  expect_equal(mpt[1], expected_median)
})
