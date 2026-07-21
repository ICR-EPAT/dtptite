# Tests for the DTP projection search strategy (Stage 1 runtime optimisation)
#
# `.dtp_project()` locates the minimum wait time at which the parent model's
# recommendation rises above the current one. `"linear"` scans every candidate
# day; `"binary"` screens the top of the range then bisects. The two must agree
# everywhere, so `"linear"` acts as the verification oracle for `"binary"`.
#
# Test groups:
# 1. Argument validation
# 2. binary == linear on fitted selectors (TITE-CRM and TITE-BOIN)
# 3. binary == linear through a full simulation
# 4. NA-at-top-of-range fallback
# 5. Admissibility invariant (never wait toward an eliminated dose)
# 6. verbose

library(escalation)

skeleton <- c(0.124, 0.25, 0.398, 0.542, 0.666)
target <- 0.25

# Random-but-reproducible outcome sets with pending follow-up, mirroring the
# shapes seen mid-trial: sorted doses, occasional DLTs, partial weights.
.random_outcomes <- function(seed, n_sets = 40) {
  set.seed(seed)
  out <- list()
  while (length(out) < n_sets) {
    m <- sample(6:18, 1)
    d <- data.frame(
      cohort = rep(seq_len(ceiling(m / 3)), each = 3)[seq_len(m)],
      dose = sort(sample(1:5, m, replace = TRUE)),
      tox = stats::rbinom(m, 1, stats::runif(1, 0.05, 0.45)),
      weight = stats::runif(m)
    )
    d$weight[d$tox == 1] <- 1
    if (any(d$tox == 0 & d$weight < 1)) out[[length(out) + 1L]] <- d
  }
  out
}

.projection_of <- function(model, outcomes) {
  f <- suppressMessages(fit(model, outcomes))
  list(
    wait = dtp_wait_time(f),
    gain = dtp_wait_time(f, type = "weight"),
    dose = dtp_projected_dose(f)
  )
}

# ===== Test 1: argument validation ========================================

test_that("apply_dtp validates projection_search and verbose", {
  base <- get_dfcrm_tite(skeleton = skeleton, target = target)

  expect_error(apply_dtp(base, t_max = 10, obswin = 56,
                         projection_search = "bogus"))
  expect_error(apply_dtp(base, t_max = 10, obswin = 56, verbose = NA))
  expect_error(apply_dtp(base, t_max = 10, obswin = 56, verbose = c(TRUE, TRUE)))

  expect_equal(apply_dtp(base, t_max = 10, obswin = 56)$projection_search,
               "binary")
  expect_true(apply_dtp(base, t_max = 10, obswin = 56)$verbose)
})

# ===== Test 2: binary == linear on fitted selectors =======================

test_that("binary and linear projections agree for TITE-CRM", {
  bin <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 35, obswin = 56, projection_search = "binary")
  lin <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 35, obswin = 56, projection_search = "linear")

  for (o in .random_outcomes(seed = 4242)) {
    expect_equal(.projection_of(bin, o), .projection_of(lin, o))
  }
})

test_that("binary and linear projections agree for TITE-BOIN", {
  bin <- get_boin_tite(5, target) |>
    apply_dtp(t_max = 35, obswin = 56, projection_search = "binary")
  lin <- get_boin_tite(5, target) |>
    apply_dtp(t_max = 35, obswin = 56, projection_search = "linear")

  for (o in .random_outcomes(seed = 99)) {
    expect_equal(.projection_of(bin, o), .projection_of(lin, o))
  }
})

test_that("binary and linear agree at t_max boundaries", {
  o <- .random_outcomes(seed = 7, n_sets = 12)
  for (t_max in c(1, 2, 35, 56, Inf)) {
    bin <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
      apply_dtp(t_max = t_max, obswin = 56, projection_search = "binary",
                verbose = FALSE)
    lin <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
      apply_dtp(t_max = t_max, obswin = 56, projection_search = "linear",
                verbose = FALSE)
    for (d in o) {
      expect_equal(.projection_of(bin, d), .projection_of(lin, d),
                   info = paste("t_max =", t_max))
    }
  }
})

test_that("a fractional t_max cannot make the two strategies disagree", {
  # t_max_effective is floored, so both strategies probe the same integer days.
  o <- .random_outcomes(seed = 21, n_sets = 12)
  bin <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 20.7, obswin = 56, projection_search = "binary")
  lin <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 20.7, obswin = 56, projection_search = "linear")
  for (d in o) expect_equal(.projection_of(bin, d), .projection_of(lin, d))
})

# ===== Test 3: binary == linear through a full simulation =================

test_that("search strategy does not change simulation results", {
  true_prob_tox <- c(0.05, 0.15, 0.25, 0.35, 0.60)
  arrivals <- function(df) cohorts_of_n(n = 3, mean_time_delta = 7)
  n_sims <- 15

  run <- function(search) {
    design <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
      apply_dtp(t_max = 35, obswin = 56, projection_search = search) |>
      stop_at_n(n = 18)
    sim_func <- simulation_function(design)
    set.seed(1)
    ps <- tite_patient_samples(n_sims, max_time = 56, num_patients = 60)
    lapply(seq_len(n_sims), function(i) {
      set.seed(5000L + i)
      sim_func(design, true_prob_tox = true_prob_tox,
               patient_sample = ps[[i]], max_time = 56,
               sample_patient_arrivals = arrivals)
    })
  }

  r_bin <- run("binary")
  r_lin <- run("linear")

  doses <- function(r) vapply(r, function(x) {
    rd <- recommended_dose(x[[1]]$fit)
    if (is.na(rd)) -1L else as.integer(rd)
  }, integer(1))
  times <- function(r) vapply(r, function(x) x[[1]]$time, numeric(1))
  waits <- function(r) dplyr::bind_rows(lapply(r, function(x) {
    tail(x, 1)[[1]]$dtp_wait_events
  }))

  expect_identical(doses(r_bin), doses(r_lin))
  expect_equal(times(r_bin), times(r_lin))
  expect_equal(waits(r_bin), waits(r_lin))
})

# ===== Test 4: stopping rules inside the projection chain =================

# The projection re-fits apply_dtp()'s parent at every candidate step, so a
# stopping rule chained *before* apply_dtp is re-evaluated there and can
# withhold a recommendation. These cover the two monitoring modes the package
# actually uses: `dose = 1` for CRM (the default) and `dose = "any"` for BOIN
# elimination. Both are invariant to projected weight -- exceedance is computed
# from raw n_at_dose()/tox_at_dose() counts, which the projection never changes
# -- so the two search strategies must agree.

.stopper_chains <- list(
  `dose = 1` = function() {
    get_dfcrm_tite(skeleton = skeleton, target = target) |>
      stop_for_beta_binomial_toxicity(
        dose = 1, tox_threshold = target + 0.1, confidence = 0.7
      )
  },
  `dose = "any"` = function() {
    get_boin_tite(5, target) |>
      stop_for_beta_binomial_toxicity(
        dose = "any", tox_threshold = target + 0.1, confidence = 0.7
      )
  }
)

test_that("binary and linear agree with a stopping rule inside the chain", {
  for (nm in names(.stopper_chains)) {
    inner <- .stopper_chains[[nm]]
    bin <- inner() |> apply_dtp(t_max = 56, obswin = 56,
                                projection_search = "binary")
    lin <- inner() |> apply_dtp(t_max = 56, obswin = 56,
                                projection_search = "linear")
    for (o in .random_outcomes(seed = 4242, n_sets = 60)) {
      expect_equal(.projection_of(bin, o), .projection_of(lin, o), info = nm)
    }
  }
})

# ===== Test 5: admissibility invariant ====================================

test_that("DTP never projects a wait toward an inadmissible dose", {
  for (nm in names(.stopper_chains)) {
    inner <- .stopper_chains[[nm]]
    for (search in c("binary", "linear")) {
      model <- inner() |> apply_dtp(t_max = 56, obswin = 56,
                                     projection_search = search)
      for (o in .random_outcomes(seed = 4242, n_sets = 60)) {
        f <- suppressMessages(fit(model, o))
        if (dtp_wait_time(f) == 0) next

        projected <- dtp_projected_dose(f)
        expect_false(is.na(projected))

        # Re-fit the parent chain at the projected wait to confirm the dose it
        # lands on is one the design still permits.
        pending <- o$tox == 0 & o$weight < 1
        at_wait <- o
        at_wait$weight[pending] <- pmin(
          1, o$weight[pending] + dtp_wait_time(f) / 56
        )
        expect_true(
          dose_admissible(suppressMessages(fit(inner(), at_wait)))[projected],
          info = paste(nm, "| search =", search)
        )
      }
    }
  }
})

# ===== Test 6: verbose ====================================================

test_that("verbose controls the effective-wait message without changing results", {
  # Needs DLTs at the top dose so the recommendation sits at or below the
  # current dose -- otherwise fit.dtp_selector_factory skips the projection
  # entirely and there is nothing to message about.
  o <- data.frame(
    cohort = rep(1:3, each = 3),
    dose = rep(c(1, 2, 3), each = 3),
    tox = c(0, 0, 0, 0, 1, 0, 1, 0, 0),
    weight = c(1, 1, 1, 1, 1, 1, 1, 0.2, 0.2)
  )
  loud <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 56, obswin = 56, verbose = TRUE)
  quiet <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 56, obswin = 56, verbose = FALSE)

  expect_message(fit(loud, o), "covers all remaining follow-up")
  expect_no_message(fit(quiet, o))
  expect_equal(.projection_of(loud, o), .projection_of(quiet, o))
})

test_that("simulation silences the effective-wait message", {
  design <- get_dfcrm_tite(skeleton = skeleton, target = target) |>
    apply_dtp(t_max = 56, obswin = 56, verbose = TRUE) |>
    stop_at_n(n = 9)

  set.seed(1)
  ps <- tite_patient_samples(1, max_time = 56, num_patients = 30)

  set.seed(202)
  expect_no_message(
    phase1_dtp_tite_sim(
      design, true_prob_tox = c(0.05, 0.15, 0.25, 0.35, 0.60),
      patient_sample = ps[[1]], max_time = 56,
      sample_patient_arrivals = function(df) {
        cohorts_of_n(n = 3, mean_time_delta = 7)
      }
    )
  )

  # ...and restores the caller's setting afterwards.
  dtp_fac <- dtptite:::.find_dtp_factory(design)
  expect_null(dtp_fac$sim_settings$verbose)
})
