# Beta-Binomial conjugate stopping rule for (DTP)-TITE-CRM
#
# Stops the trial when the posterior probability of excess toxicity
# at a monitored dose exceeds a confidence threshold.
#
# Prior:     p ~ Beta(a, b)
# Posterior: p | data ~ Beta(a + s, b + n - s)
# Stop when: Pr(p > tox_threshold | data) > confidence

# -- Factory -----------------------------------------------------------

#' Stop for excessive toxicity using a Beta-Binomial conjugate posterior.
#'
#' Creates a stopping decorator that monitors one or more doses and halts
#' the trial when the posterior probability of the true toxicity rate
#' exceeding `tox_threshold` is greater than `confidence`.
#'
#' @param parent_selector_factory A `selector_factory` object.
#' @param dose Which dose(s) to monitor. A numeric scalar (e.g. `1`),
#'   numeric vector (e.g. `c(1, 2)`), `"recommended"`, or `"any"`.
#'   Default: `1`.
#' @param tox_threshold Toxicity probability threshold. Required.
#' @param confidence Posterior exceedance cutoff (default: `0.88`).
#' @param a,b Beta prior hyperparameters. When `NULL` (default), resolved
#'   at fit time from `tox_target()` of the parent selector:
#'   `a = tox_target, b = 1 - tox_target`.
#' @return A `beta_binom_tox_selector_factory` object.
#' @export
stop_for_beta_binomial_toxicity <- function(parent_selector_factory,
                                            dose = 1,
                                            tox_threshold,
                                            confidence = 0.88,
                                            a = NULL,
                                            b = NULL) {
  if (is.character(dose)) {
    if (!identical(dose, "recommended") && !identical(dose, "any")) {
      stop("'dose' must be numeric, \"recommended\", or \"any\"")
    }
  } else {
    stopifnot(is.numeric(dose), length(dose) >= 1, all(dose >= 1))
  }
  stopifnot(is.numeric(tox_threshold), length(tox_threshold) == 1,
            tox_threshold > 0, tox_threshold < 1)
  stopifnot(is.numeric(confidence), length(confidence) == 1,
            confidence > 0, confidence < 1)
  if (!is.null(a)) stopifnot(is.numeric(a), length(a) == 1, a > 0)
  if (!is.null(b)) stopifnot(is.numeric(b), length(b) == 1, b > 0)

  x <- list(
    parent = parent_selector_factory,
    dose = dose,
    tox_threshold = tox_threshold,
    confidence = confidence,
    a = a,
    b = b
  )
  class(x) <- c("beta_binom_tox_selector_factory",
                 "derived_dose_selector_factory",
                 "selector_factory")
  x
}

# -- Fit ---------------------------------------------------------------

#' @export
fit.beta_binom_tox_selector_factory <- function(selector_factory,
                                                outcomes, ...) {
  parent_selector <- selector_factory$parent |> fit(outcomes, ...)

  a <- selector_factory$a
  b <- selector_factory$b
  if (is.null(a)) a <- tox_target(parent_selector)
  if (is.null(b)) b <- 1 - tox_target(parent_selector)
  stopifnot(a > 0, b > 0)

  beta_binom_tox_selector(
    parent_selector = parent_selector,
    dose = selector_factory$dose,
    tox_threshold = selector_factory$tox_threshold,
    confidence = selector_factory$confidence,
    a = a,
    b = b
  )
}

# -- Selector constructor ----------------------------------------------

beta_binom_tox_selector <- function(parent_selector, dose, tox_threshold,
                                    confidence, a, b) {
  l <- list(
    parent = parent_selector,
    dose = dose,
    tox_threshold = tox_threshold,
    confidence = confidence,
    a = a,
    b = b
  )
  class(l) <- c("beta_binom_tox_selector",
                 "derived_dose_selector",
                 "selector")
  l
}

# -- Internal helpers ---------------------------------------------------

.bb_exceedance_probs <- function(x) {
  nd <- num_doses(x)
  n_at <- n_at_dose(x)
  s_at <- tox_at_dose(x)

  exceedance <- rep(NA_real_, nd)
  for (j in seq_len(nd)) {
    if (n_at[j] == 0) next
    exceedance[j] <- 1 - pbeta(x$tox_threshold,
                                x$a + s_at[j],
                                x$b + n_at[j] - s_at[j])
  }
  exceedance
}

.bb_should_stop <- function(x) {
  exceedance <- .bb_exceedance_probs(x)
  dose <- x$dose

  if (is.character(dose)) {
    if (dose == "any") {
      # Elimination mode: individual doses are marked inadmissible via
      # dose_admissible(); the trial stops only when no admissible dose
      # can be recommended (i.e. all reachable doses are eliminated).
      admissible <- dose_admissible(x)
      rec <- recommended_dose(x$parent)
      if (is.na(rec)) return(TRUE)  # nocov
      return(!any(admissible[seq_len(rec)]))
    }
    if (dose == "recommended") {
      rec <- recommended_dose(x$parent)
      if (is.na(rec)) return(FALSE)  # nocov
      dose <- rec
    }
  }

  for (d in dose) {
    if (d >= 1 && d <= length(exceedance) &&
        !is.na(exceedance[d]) && exceedance[d] >= x$confidence) {
      return(TRUE)
    }
  }
  FALSE
}

# -- S3 methods ---------------------------------------------------------

#' @export
continue.beta_binom_tox_selector <- function(x, ...) {
  if (!continue(x$parent, ...)) return(FALSE)
  if (.bb_should_stop(x)) return(FALSE)
  TRUE
}

#' @export
recommended_dose.beta_binom_tox_selector <- function(x, ...) {
  if (.bb_should_stop(x)) return(NA)
  rec <- recommended_dose(x$parent, ...)
  if (is.na(rec)) return(NA)  # nocov

  admissible <- dose_admissible(x)
  if (admissible[rec]) return(rec)

  # Clamp to highest admissible dose at or below the recommendation
  candidates <- which(admissible[seq_len(rec)])
  if (length(candidates) == 0L) return(NA_integer_)  # nocov
  as.integer(max(candidates))
}

#' @export
dose_admissible.beta_binom_tox_selector <- function(x, ...) {
  exceedance <- .bb_exceedance_probs(x)
  parent_admissible <- dose_admissible(x$parent, ...)
  inadmissible <- !is.na(exceedance) & exceedance >= x$confidence

  # Cascade: if dose d is eliminated, all doses above d are also eliminated
  # (matches BOIN/TITEgBOIN elimination rule)
  first_elim <- which(inadmissible)[1]
  if (!is.na(first_elim) && first_elim < length(inadmissible)) {
    inadmissible[first_elim:length(inadmissible)] <- TRUE
  }

  parent_admissible & !inadmissible
}
