# Tests for the time-to-first-MTD simulation metric
#
# Coverage:
#  (a) .patient_dose_times() finds dose/time across every design class, and
#      fails loudly when it cannot — the escalation-internals regression guard
#  (b) scenario_mtd() band rule, no simulation needed
#  (c) first-MTD event logic on synthetic frames, no simulation needed
#  (d) Censoring arithmetic: denominator, reached fraction, undefined median
#  (e) No-MTD scenarios yield NA rather than FALSE
#  (g) simulations_collection shape
#  (h) plot() axes and smoke test

library(escalation)

skeleton <- c(0.124, 0.25, 0.398, 0.542, 0.666)
target <- 0.25
arrivals3 <- function(df) escalation::cohorts_of_n(n = 3, mean_time_delta = 1)
true_prob_tox <- c(0.05, 0.15, 0.25, 0.35, 0.60)   # MTD = dose 3
all_toxic <- c(0.60, 0.70, 0.80, 0.90, 0.95)       # no MTD

# ===== (a) patient dose/time frame across design classes ==============

test_that(".patient_dose_times finds dose and time for every design class", {
  designs <- list(
    `dtp+tite-crm` = get_dfcrm_tite(skeleton = skeleton, target = target) |>
      apply_dtp(t_max = 20, obswin = 56) |> stop_at_n(n = 9),
    `tite-crm` = get_dfcrm_tite(skeleton = skeleton, target = target) |>
      stop_at_n(n = 9),
    `dtp+tite-boin` = get_boin_tite(num_doses = 5, target = target) |>
      apply_dtp(t_max = 20, obswin = 56) |> stop_at_n(n = 9)
  )

  for (nm in names(designs)) {
    set.seed(456)
    sims <- suppressWarnings(simulate_trials(
      designs[[nm]], num_sims = 2, true_prob_tox = true_prob_tox,
      sample_patient_arrivals = arrivals3, max_time = 56
    ))
    fit <- sims$fits[[1]][[length(sims$fits[[1]])]]$fit
    pdt <- dtptite:::.patient_dose_times(fit)
    expect_true(all(c("dose", "time") %in% names(pdt)), info = nm)
    expect_gt(nrow(pdt), 0)
    expect_false(any(is.na(pdt$time)), info = nm)
  }
})

test_that(".patient_dose_times covers non-TITE designs too", {
  for (design in list(get_dfcrm(skeleton = skeleton, target = target) |>
                        stop_at_n(n = 9),
                      get_boin(num_doses = 5, target = target) |>
                        stop_at_n(n = 9))) {
    set.seed(456)
    sims <- simulate_trials(design, num_sims = 2,
                            true_prob_tox = true_prob_tox)
    fit <- sims$fits[[1]][[length(sims$fits[[1]])]]$fit
    pdt <- dtptite:::.patient_dose_times(fit)
    expect_true(all(c("dose", "time") %in% names(pdt)))
  }
})

test_that(".patient_dose_times errors loudly when no frame is present", {
  expect_error(dtptite:::.patient_dose_times(list(df = data.frame(dose = 1))),
               "Could not find a patient-level data frame")
  expect_error(dtptite:::.patient_dose_times(NULL),
               "Could not find a patient-level data frame")
})

# ===== (b) scenario_mtd ================================================

test_that("scenario_mtd picks the dose within tolerance", {
  expect_equal(scenario_mtd(true_prob_tox, 0.25), 3L)
  expect_type(scenario_mtd(true_prob_tox, 0.25), "integer")
})

test_that("scenario_mtd returns the lowest of several qualifying doses", {
  expect_equal(scenario_mtd(c(0.05, 0.22, 0.27, 0.50), 0.25), 2L)
})

test_that("scenario_mtd treats the tolerance boundary as inclusive", {
  expect_equal(scenario_mtd(c(0.05, 0.30, 0.60), 0.25, tol = 0.05), 2L)
  expect_warning(scenario_mtd(c(0.05, 0.3001, 0.60), 0.25, tol = 0.05),
                 "no dose lies within tol")
})

test_that("scenario_mtd reports no MTD when every dose is too toxic", {
  expect_warning(res <- scenario_mtd(all_toxic, 0.25), "exceeds target")
  expect_true(is.na(res))
})

test_that("scenario_mtd reports no MTD when every dose is too safe", {
  expect_warning(res <- scenario_mtd(c(0.01, 0.02, 0.03), 0.25),
                 "falls below target")
  expect_true(is.na(res))
})

test_that("scenario_mtd handles non-monotone scenarios by taking the lowest", {
  expect_equal(scenario_mtd(c(0.26, 0.10, 0.24), 0.25), 1L)
})

test_that("scenario_mtd validates its arguments", {
  expect_error(scenario_mtd(numeric(0), 0.25), "non-empty numeric")
  expect_error(scenario_mtd(true_prob_tox, c(0.2, 0.3)), "single finite")
  expect_error(scenario_mtd(true_prob_tox, 0.25, tol = -1), "non-negative")
})

test_that("scenario_mtd rejects missing or infinite toxicity probabilities", {
  expect_error(scenario_mtd(c(0.05, NA, 0.25), 0.25), "missing or infinite")
  expect_error(scenario_mtd(c(0.05, Inf, 0.25), 0.25), "missing or infinite")
})

test_that("scenario_mtd rejects toxicity probabilities outside 0 and 1", {
  expect_error(scenario_mtd(c(0.05, 1.4), 0.25), "between 0 and 1")
  expect_error(scenario_mtd(c(-0.1, 0.25), 0.25), "between 0 and 1")
})

test_that("scenario_mtd rejects an out-of-range or non-finite target", {
  expect_error(scenario_mtd(true_prob_tox, NA), "single finite")
  expect_error(scenario_mtd(true_prob_tox, Inf), "single finite")
  expect_error(scenario_mtd(true_prob_tox, 1.5), "between 0 and 1")
})

test_that("scenario_mtd rejects a non-finite tol", {
  expect_error(scenario_mtd(true_prob_tox, 0.25, tol = Inf), "finite")
})

# ===== (c) first-MTD event on synthetic frames =========================

test_that("first-MTD event picks the earliest patient at the MTD", {
  p <- data.frame(dose = c(1, 1, 3, 2, 3), time = c(1, 2, 5, 3, 9))
  ev <- dtptite:::.first_mtd_event(p, 3L)
  expect_true(ev$reached)
  expect_equal(ev$time, 5)
  expect_equal(ev$n_before, 3L)   # times 1, 2, 3 precede time 5
})

test_that("first-MTD event returns NA when the MTD was never given", {
  p <- data.frame(dose = c(1, 1, 2, 2), time = c(1, 2, 3, 4))
  ev <- dtptite:::.first_mtd_event(p, 3L)
  expect_false(ev$reached)
  expect_true(is.na(ev$time))
  expect_true(is.na(ev$n_before))
})

test_that("first-MTD event is unaffected by ties at a shared timestamp", {
  # DTP doses a queued cohort at one instant
  p <- data.frame(dose = c(1, 1, 1, 3, 3, 3), time = c(1, 2, 3, 8, 8, 8))
  ev <- dtptite:::.first_mtd_event(p, 3L)
  expect_equal(ev$time, 8)
  expect_equal(ev$n_before, 3L)
})

test_that("first-MTD event returns NA throughout when there is no MTD", {
  p <- data.frame(dose = c(1, 2), time = c(1, 2))
  ev <- dtptite:::.first_mtd_event(p, NA_integer_)
  expect_true(is.na(ev$reached))
  expect_true(is.na(ev$time))
})

# ===== (d) censoring arithmetic ========================================

run_no_dtp <- function(num_sims = 60) {
  set.seed(2026)
  simulate_trials(
    get_dfcrm_tite(skeleton = skeleton, target = target) |> stop_at_n(n = 18),
    num_sims = num_sims, true_prob_tox = true_prob_tox,
    sample_patient_arrivals = arrivals3, max_time = 56
  )
}

test_that("every replicate is retained, reached or not", {
  r <- time_to_first_mtd(run_no_dtp())
  expect_s3_class(r, "time_to_first_mtd")
  expect_equal(nrow(r), 60)
  expect_equal(attr(r, "num_sims"), 60L)
  expect_true(any(!r$reached))          # some never reach
  expect_equal(r$reached, !is.na(r$time_to_mtd))
})

test_that("reached_fraction uses all replicates as the denominator", {
  r <- time_to_first_mtd(run_no_dtp())
  s <- summary(r)
  expect_equal(s$num_sims, 60L)
  expect_equal(s$num_reached, sum(r$reached))
  expect_equal(s$reached_fraction, sum(r$reached) / 60)
  expect_lt(s$reached_fraction, 0.5)    # fixture premise for the next test
})

test_that("the median is undefined when fewer than half of trials reach", {
  r <- time_to_first_mtd(run_no_dtp())
  s <- summary(r)
  expect_true(is.na(s$median_time_to_mtd))
  # but the conditional median over reached replicates does exist, which is
  # precisely the quantity this must not report
  expect_false(is.na(median(r$time_to_mtd, na.rm = TRUE)))
})

test_that("quantiles below the reached fraction are estimable", {
  r <- time_to_first_mtd(run_no_dtp())
  s <- summary(r, probs = c(0.1, 0.9))
  expect_false(is.na(s$q10_time_to_mtd))
  expect_true(is.na(s$q90_time_to_mtd))
  expect_true(all(c("q10_time_to_mtd", "q90_time_to_mtd") %in% names(s)))
})

test_that("summary warns when the object has been subset", {
  r <- time_to_first_mtd(run_no_dtp())
  expect_warning(summary(r[1:10, ]), "appears to have been subset")
})

test_that("summary validates probs", {
  r <- time_to_first_mtd(run_no_dtp())
  expect_error(summary(r, probs = 1.5), "between 0 and 1")
  expect_error(summary(r, probs = c(0.1, NA)), "finite numbers")
})

test_that("summary rejects duplicate probs", {
  r <- time_to_first_mtd(run_no_dtp())
  expect_error(summary(r, probs = c(0.1, 0.5, 0.1)), "duplicate values")
})

test_that("fractional percentages get distinct quantile columns", {
  r <- time_to_first_mtd(run_no_dtp())
  s <- summary(r, probs = c(0.125, 0.126))
  expect_true(all(c("q12.5_time_to_mtd", "q12.6_time_to_mtd") %in% names(s)))
  # integer percentages keep their existing, unrounded-looking names
  expect_equal(names(summary(r, probs = 0.9)) |> tail(1), "q90_time_to_mtd")
})

# ===== (e) scenarios with no MTD =======================================

test_that("a scenario with no MTD yields NA, not FALSE", {
  set.seed(11)
  sims <- simulate_trials(
    get_dfcrm_tite(skeleton = skeleton, target = target) |> stop_at_n(n = 12),
    num_sims = 5, true_prob_tox = all_toxic,
    sample_patient_arrivals = arrivals3, max_time = 56
  )
  r <- suppressWarnings(time_to_first_mtd(sims))
  expect_true(all(is.na(r$mtd)))
  expect_true(all(is.na(r$reached)))
  expect_equal(nrow(r), 5)

  s <- summary(r)
  expect_true(is.na(s$reached_fraction))
  expect_true(is.na(s$median_time_to_mtd))
  expect_equal(s$num_sims, 5L)
})

test_that("an explicit mtd overrides derivation and is validated", {
  sims <- run_no_dtp(num_sims = 4)
  expect_equal(unique(time_to_first_mtd(sims, mtd = 2)$mtd), 2L)
  expect_true(all(is.na(time_to_first_mtd(sims, mtd = NA)$mtd)))
  expect_error(time_to_first_mtd(sims, mtd = 99), "between 1 and 5")
  expect_error(time_to_first_mtd(sims, mtd = c(1, 2)), "single value")
})

test_that("a non-whole mtd is rejected rather than truncated", {
  sims <- run_no_dtp(num_sims = 4)
  expect_error(time_to_first_mtd(sims, mtd = 2.9), "whole number")
  expect_error(time_to_first_mtd(sims, mtd = "2"), "must be a number")
  expect_error(time_to_first_mtd(sims, mtd = Inf), "must be finite")
  # whole numbers supplied as doubles still work
  expect_equal(unique(time_to_first_mtd(sims, mtd = 2.0)$mtd), 2L)
})

# ===== (g) simulations_collection ======================================

test_that("a collection gains a design column and per-design summaries", {
  set.seed(2026)
  sc <- simulate_compare(
    list(
      DTP = get_dfcrm_tite(skeleton = skeleton, target = target) |>
        apply_dtp(t_max = 20, obswin = 56) |> stop_at_n(n = 18),
      plain = get_dfcrm_tite(skeleton = skeleton, target = target) |>
        stop_at_n(n = 18)
    ),
    num_sims = 20, true_prob_tox = true_prob_tox,
    sample_patient_arrivals = arrivals3, max_time = 56
  )
  r <- suppressWarnings(time_to_first_mtd(sc))
  expect_true("design" %in% names(r))
  expect_setequal(unique(r$design), c("DTP", "plain"))
  expect_equal(nrow(r), 40)

  s <- summary(r)
  expect_equal(nrow(s), 2)
  expect_true("design" %in% names(s))
  expect_true(all(s$num_sims == 20L))
})

# ===== (h) plot ========================================================

test_that("plot returns a ggplot on both axes", {
  skip_if_not_installed("ggplot2")
  r <- time_to_first_mtd(run_no_dtp(num_sims = 20))
  expect_s3_class(plot(r), "ggplot")
  expect_s3_class(plot(r, axis = "patients"), "ggplot")
  expect_s3_class(plot(r, ci = TRUE), "ggplot")
  expect_error(plot(r, axis = "nonsense"))
  expect_error(plot(r, ci = "yes"), "TRUE or FALSE")
  expect_error(plot(r, ci = NA), "TRUE or FALSE")
})
