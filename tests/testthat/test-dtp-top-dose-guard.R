# Tests for the top-dose short-circuit in `fit.dtp_selector_factory()`
#
# `.dtp_project()` reports a wait only when the projected dose rises strictly
# above `dose_recommended`. When the recommendation is already the top dose no
# amount of extra follow-up can satisfy that, so the projection is skipped
# outright. The guard must not catch a de-escalation off the top dose, where
# recovering to the top dose is a genuine hit.

library(escalation)

target <- 0.25
num_doses <- 5

.boin_dtp <- function(t_max = 56) {
  get_boin_tite(num_doses, target, p.saf = 0.15, p.tox = 0.35) |>
    apply_dtp(t_max = t_max, obswin = 56, verbose = FALSE)
}

test_that("at the top dose with a stay recommendation, DTP returns no wait", {
  # Escalating through to dose 5; the recommendation clamps at `num_doses`, so
  # dose_recommended == dose_current == 5.
  df <- data.frame(
    dose   = c(1, 2, 3, 4, 5),
    tox    = c(0, 0, 0, 0, 0),
    weight = c(1, 1, 1, 0.5, 0.25),
    cohort = c(1, 2, 3, 4, 5)
  )
  fit_obj <- .boin_dtp() |> fit(df)

  expect_equal(recommended_dose(fit_obj), num_doses)
  expect_false(dtp_should_wait(fit_obj))
  expect_equal(dtp_wait_time(fit_obj), 0)
  expect_equal(dtp_wait_time(fit_obj, type = "weight"), 0)
  expect_true(is.na(dtp_projected_dose(fit_obj)))
})

test_that("the guard still fires when pending follow-up exists at the top dose", {
  # Same as above but every dose-5 patient is pending, so there is plenty of
  # weight left to project — the skip must be driven by the recommendation
  # being the top dose, not by an absence of pending patients.
  df <- data.frame(
    dose   = c(1, 1, 1, 5, 5, 5),
    tox    = c(0, 0, 0, 0, 0, 0),
    weight = c(1, 1, 1, 0.25, 0.25, 0.25),
    cohort = c(1, 1, 1, 2, 2, 2)
  )
  fit_obj <- .boin_dtp() |> fit(df)

  expect_equal(recommended_dose(fit_obj), num_doses)
  expect_false(dtp_should_wait(fit_obj))
  expect_equal(dtp_wait_time(fit_obj), 0)
})

test_that("a de-escalation off the top dose still projects a wait", {
  # Dose 5 holds 1 DLT and 3 pending patients at weight 0.5:
  #   n0 = 1 + 3 * 0.5 = 2.5, phat = 1 / 2.5 = 0.4 >= lambda_d -> de-escalate.
  # Projecting the pending patients to full follow-up gives n0 = 4 and
  # phat = 0.25, which sits between lambda_e and lambda_d -> stay at dose 5.
  # Recovering to dose 5 is strictly above the recommended dose 4, so the
  # projection hits and the guard must not suppress it.
  df <- data.frame(
    dose   = c(1, 2, 3, 4, 5, 5, 5, 5),
    tox    = c(0, 0, 0, 0, 1, 0, 0, 0),
    weight = c(1, 1, 1, 1, 1, 0.5, 0.5, 0.5),
    cohort = c(1, 2, 3, 4, 5, 5, 5, 5)
  )
  fit_obj <- .boin_dtp() |> fit(df)

  expect_equal(recommended_dose(fit_obj), 4L)
  expect_true(dtp_should_wait(fit_obj))
  expect_equal(dtp_wait_time(fit_obj), 16)
  expect_equal(dtp_projected_dose(fit_obj), num_doses)
})

test_that("the guard is inert below the top dose", {
  # Identical shape to the de-escalation case but one dose lower, so neither
  # the recommendation nor the current dose is the top dose.
  df <- data.frame(
    dose   = c(1, 1, 1, 2, 2, 2, 2),
    tox    = c(0, 0, 0, 1, 0, 0, 0),
    weight = c(1, 1, 1, 1, 0.5, 0.5, 0.5),
    cohort = c(1, 1, 1, 2, 2, 2, 3)
  )
  fit_obj <- .boin_dtp() |> fit(df)

  expect_equal(recommended_dose(fit_obj), 1L)
  expect_true(dtp_should_wait(fit_obj))
  expect_equal(dtp_wait_time(fit_obj), 16)
  expect_equal(dtp_projected_dose(fit_obj), 2L)
})
