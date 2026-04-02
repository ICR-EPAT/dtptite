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

test_that("1d: Dose elimination via decorator matches TITEgBOIN", {
  skip_if_not_installed("TITEgBOIN")

  # BOIN elimination uses Beta(0.5, 0.5) Jeffreys prior with cutoff.eli = 0.95

  # --- Scenario A: dose 1 too toxic -> stop trial ---
  df_a <- data.frame(
    dose = c(1, 1, 1), tox = c(1, 1, 1),
    weight = c(1, 1, 1), cohort = 1:3
  )

  tg_a <- TITEgBOIN::next_TITE_QuasiBOIN(
    target = target, n = c(3, 0, 0, 0, 0), npend = c(0, 0, 0, 0, 0),
    y = c(3, 0, 0, 0, 0), ft = c(0, 0, 0, 0, 0), d = 1, maxt = 28,
    cutoff.eli = 0.95, Neli = 3, gdesign = FALSE, print_d = FALSE
  )
  expect_equal(tg_a$d, 0)
  expect_equal(tg_a$earlystop, 1)

  model_a <- get_boin_tite(5, target) %>%
    stop_for_beta_binomial_toxicity(dose = 1, tox_threshold = target,
                                    confidence = 0.95, a = 0.5, b = 0.5)
  fit_a <- model_a %>% fit(df_a)

  expect_false(continue(fit_a))
  expect_true(is.na(recommended_dose(fit_a)))

  # --- Scenario B: dose 2 toxic, dose 1 safe -> continue at dose 1 ---
  # Last cohort at dose 2; bare model deescalates to 1.
  df_b <- data.frame(
    dose = c(1, 1, 1, 2, 2, 2),
    tox = c(0, 0, 0, 1, 1, 1),
    weight = rep(1, 6),
    cohort = c(1, 1, 1, 2, 2, 2)
  )

  tg_b <- TITEgBOIN::next_TITE_QuasiBOIN(
    target = target, n = c(3, 3, 0, 0, 0), npend = c(0, 0, 0, 0, 0),
    y = c(0, 3, 0, 0, 0), ft = c(0, 0, 0, 0, 0), d = 2, maxt = 28,
    cutoff.eli = 0.95, Neli = 3, gdesign = FALSE, print_d = FALSE
  )
  expect_equal(tg_b$d, 1)

  model_b <- get_boin_tite(5, target) %>%
    stop_for_beta_binomial_toxicity(dose = 1, tox_threshold = target,
                                    confidence = 0.95, a = 0.5, b = 0.5)
  fit_b <- model_b %>% fit(df_b)

  expect_true(continue(fit_b))
  expect_equal(recommended_dose(fit_b), tg_b$d)  # Both recommend dose 1
  expect_false(dose_admissible(fit_b)[2])         # dose 2 eliminated

  # --- Scenario C: dose 2 toxic, last cohort at dose 1 -> still dose 1 ---
  # Bare model would escalate from dose 1 to inadmissible dose 2;
  # recommended_dose must clamp to the highest admissible dose.
  df_c <- data.frame(
    dose = c(2, 2, 2, 1, 1, 1),
    tox = c(1, 1, 1, 0, 0, 0),
    weight = rep(1, 6),
    cohort = c(1, 1, 1, 2, 2, 2)
  )

  tg_c <- TITEgBOIN::next_TITE_QuasiBOIN(
    target = target, n = c(3, 3, 0, 0, 0), npend = c(0, 0, 0, 0, 0),
    y = c(0, 3, 0, 0, 0), ft = c(0, 0, 0, 0, 0), d = 1, maxt = 28,
    cutoff.eli = 0.95, Neli = 3, gdesign = FALSE, print_d = FALSE
  )
  expect_equal(tg_c$d, 1)

  model_c <- get_boin_tite(5, target) %>%
    stop_for_beta_binomial_toxicity(dose = 1, tox_threshold = target,
                                    confidence = 0.95, a = 0.5, b = 0.5)
  fit_c <- model_c %>% fit(df_c)

  expect_true(continue(fit_c))
  expect_equal(recommended_dose(fit_c), tg_c$d)  # Clamped to dose 1, not 2
  expect_false(dose_admissible(fit_c)[2])         # dose 2 eliminated
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

test_that("3a: Escalation matches TITEgBOIN with pending patients", {
  skip_if_not_installed("TITEgBOIN")
  maxt <- 28

  # Dose 1: 3 complete, 0 DLTs. Dose 2: 2 complete + 1 pending (w=0.3), 0 DLTs.
  # Current dose = 2. phat = 0/2.3 = 0 <= lambda_e -> escalate to 3.
  df <- data.frame(
    dose   = c(1, 1, 1, 2, 2, 2),
    tox    = c(0, 0, 0, 0, 0, 0),
    weight = c(1, 1, 1, 1, 1, 0.3),
    cohort = c(1, 1, 1, 2, 2, 3)
  )
  our_rec <- recommended_dose(get_boin_tite(5, target) %>% fit(df))

  tg <- TITEgBOIN::next_TITE_QuasiBOIN(
    target = target,
    n = c(3, 3, 0, 0, 0), npend = c(0, 1, 0, 0, 0),
    y = c(0, 0, 0, 0, 0), ft = c(0, 0.3 * maxt, 0, 0, 0),
    d = 2, maxt = maxt, maxpen = 1,
    gdesign = FALSE, print_d = FALSE
  )

  expect_equal(our_rec, tg$d)
})

test_that("3b: Stay decision matches TITEgBOIN with pending patients", {
  skip_if_not_installed("TITEgBOIN")
  maxt <- 28

  # Dose 2: 3 complete + 1 pending (w=0.5), 1 DLT.
  # n0 = 3 + 0.5 = 3.5, phat = 1/3.5 = 0.286
  # lambda_e < 0.286 < lambda_d -> stay at dose 2.
  df <- data.frame(
    dose   = c(2, 2, 2, 2),
    tox    = c(0, 1, 0, 0),
    weight = c(1, 1, 1, 0.5),
    cohort = c(1, 1, 1, 2)
  )
  our_rec <- recommended_dose(get_boin_tite(5, target) %>% fit(df))

  tg <- TITEgBOIN::next_TITE_QuasiBOIN(
    target = target,
    n = c(0, 4, 0, 0, 0), npend = c(0, 1, 0, 0, 0),
    y = c(0, 1, 0, 0, 0), ft = c(0, 0.5 * maxt, 0, 0, 0),
    d = 2, maxt = maxt, maxpen = 1,
    gdesign = FALSE, print_d = FALSE
  )

  expect_equal(our_rec, tg$d)
})

test_that("3c: De-escalation matches TITEgBOIN with pending patients", {
  skip_if_not_installed("TITEgBOIN")
  maxt <- 28

  # Dose 1: 3 complete, 0 DLTs. Dose 2: 3 complete + 1 pending (w=0.6), 2 DLTs.
  # n0 = 3 + 0.6 = 3.6, phat = 2/3.6 = 0.556 >= lambda_d -> deescalate to 1.
  df <- data.frame(
    dose   = c(1, 1, 1, 2, 2, 2, 2),
    tox    = c(0, 0, 0, 1, 1, 0, 0),
    weight = c(1, 1, 1, 1, 1, 1, 0.6),
    cohort = c(1, 1, 1, 2, 2, 2, 3)
  )
  our_rec <- recommended_dose(get_boin_tite(5, target) %>% fit(df))

  tg <- TITEgBOIN::next_TITE_QuasiBOIN(
    target = target,
    n = c(3, 4, 0, 0, 0), npend = c(0, 1, 0, 0, 0),
    y = c(0, 2, 0, 0, 0), ft = c(0, 0.6 * maxt, 0, 0, 0),
    d = 2, maxt = maxt, maxpen = 1,
    gdesign = FALSE, print_d = FALSE
  )

  expect_equal(our_rec, tg$d)
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

test_that("5b: dont_skip_doses prevents jumping past untried dose 2", {
  # Cohort 1 at dose 3 (0/3 DLTs), Cohort 2 at dose 1 (0/3 DLTs).
  # Dose 2 is untried. Last dose given = 1.
  df <- data.frame(
    dose   = c(3, 3, 3, 1, 1, 1),
    tox    = c(0, 0, 0, 0, 0, 0),
    weight = rep(1, 6),
    cohort = c(1, 1, 1, 2, 2, 2)
  )

  # Without dont_skip: select_boin_mtd(when="always") recommends isotonic MTD = 3
  model_no_skip <- get_boin_tite(5, target) %>%
    select_boin_mtd(when = "always")
  fit_no_skip <- model_no_skip %>% fit(df)
  expect_equal(recommended_dose(fit_no_skip), 3L)

  # With dont_skip: constrained to last_d + 1 = 1 + 1 = 2
  model_skip <- get_boin_tite(5, target) %>%
    select_boin_mtd(when = "always") %>%
    dont_skip_doses()
  fit_skip <- model_skip %>% fit(df)
  expect_equal(recommended_dose(fit_skip), 2L)
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

test_that("5d: select_boin_mtd overrides next-dose with isotonic MTD at trial end", {
  # Dose 1: 0/3, Dose 2: 0/3, Dose 3: 1/3
  # Bare model: last_dose = 3, phat = 1/3 = 0.333 >= lambda_d -> deescalate to 2
  # Isotonic MTD: dose 3 (closest to target = 0.25 from below)
  df <- data.frame(
    dose   = rep(1:3, each = 3),
    tox    = c(0, 0, 0, 0, 0, 0, 0, 0, 1),
    weight = rep(1, 9),
    cohort = rep(1:3, each = 3)
  )

  # Bare model recommends dose 2 (deescalate from 3)
  bare <- get_boin_tite(5, target) %>% fit(df)
  expect_equal(recommended_dose(bare), 2L)

  # With stop_at_n(9) + select_boin_mtd: isotonic MTD overrides to dose 3
  model <- get_boin_tite(5, target) %>%
    stop_at_n(9) %>%
    select_boin_mtd()
  fit_obj <- model %>% fit(df)

  expect_false(continue(fit_obj))
  mtd <- recommended_dose(fit_obj)
  expect_equal(mtd, 3L)
  # MTD differs from the next-dose recommendation
  expect_false(identical(mtd, recommended_dose(bare)))
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
