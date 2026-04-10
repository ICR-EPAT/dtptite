# TITE-BOIN selector — BOIN decision function and factory
#
# Provides the BOIN-specific lambda-based dose decision for the shared
# TITE-MAD engine. This is a barebone selector: no elimination, no
# suspension, no stopping rules — those are composed via decorators.

# -- Decision function -------------------------------------------------------

.boin_tite_decision <- function(s, n0, args) {
  if (n0 == 0) return("stay")
  phat <- s / n0
  if (phat <= args$lambda_e) return("escalate")
  if (phat >= args$lambda_d) return("deescalate")
  "stay"
}

# -- Factory -----------------------------------------------------------------

#' Create a TITE-BOIN dose-finding selector.
#'
#' Returns a selector factory for the BOIN design with TITE
#' (time-to-event) support via effective sample sizes. This is a barebone
#' selector: stopping, elimination, and overdose control are added via
#' piping (e.g. `stop_for_beta_binomial_toxicity`, `dont_skip_doses`).
#'
#' @param num_doses Number of doses under investigation.
#' @param target Target toxicity probability.
#' @param p.saf Safety boundary parameter (default `0.6 * target`).
#' @param p.tox Toxicity boundary parameter (default `1.4 * target`).
#' @param ... Extra arguments stored for future use.
#'
#' @return A `boin_tite_selector_factory` object.
#' @export
#'
#' @examples
#' # Barebone — pure dose decision
#' library(escalation)
#' get_boin_tite(5, 0.25) %>%
#'   fit(data.frame(dose = c(1,1,1), tox = c(0,0,0),
#'                  weight = c(1,1,1), cohort = c(1,1,1)))
get_boin_tite <- function(num_doses, target,
                          p.saf = 0.6 * target,
                          p.tox = 1.4 * target, ...) {
  stopifnot(is.numeric(num_doses), length(num_doses) == 1, num_doses >= 1)
  stopifnot(is.numeric(target), length(target) == 1, target > 0, target < 1)
  stopifnot(p.saf < target, p.tox > target)

  lambda_e <- log((1 - p.saf) / (1 - target)) /
    log(target * (1 - p.saf) / (p.saf * (1 - target)))
  lambda_d <- log((1 - target) / (1 - p.tox)) /
    log(p.tox * (1 - target) / (target * (1 - p.tox)))

  x <- list(
    num_doses = num_doses,
    target = target,
    lambda_e = lambda_e,
    lambda_d = lambda_d,
    extra_args = list(...)
  )
  class(x) <- c("boin_tite_selector_factory",
                 "tite_mad_selector_factory",
                 "tox_selector_factory",
                 "selector_factory")
  x
}

# -- Fit method --------------------------------------------------------------

#' @export
fit.boin_tite_selector_factory <- function(selector_factory, outcomes, ...) {
  tite_mad_selector(
    outcomes = outcomes,
    num_doses = selector_factory$num_doses,
    target = selector_factory$target,
    decision_fn = .boin_tite_decision,
    decision_args = list(
      lambda_e = selector_factory$lambda_e,
      lambda_d = selector_factory$lambda_d
    ),
    design_label = "boin_tite",
    ...
  )
}

# -- Simulation dispatch -----------------------------------------------------

#' @export
simulation_function.tite_mad_selector_factory <- function(selector_factory) {
  escalation:::phase1_tite_sim
}
