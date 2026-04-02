# Shared TITE engine for Model-Assisted Designs (MAD)
#
# Provides the common weight→ESS→decision logic for TITE variants of
# model-assisted dose-finding designs (BOIN, keyboard, mTPI). Each design
# supplies a pluggable decision function; all other machinery is shared.
#
# "MAD" clarifies this engine is for model-assisted designs only —
# model-based TITE (CRM, BLRM) uses weights in the likelihood, not ESS.

# -- Helper: per-dose summary from TITE data --------------------------------

.tite_dose_summary <- function(df, num_doses) {
  summary <- data.frame(
    dose = seq_len(num_doses),
    n    = integer(num_doses),
    s    = integer(num_doses),
    r    = integer(num_doses),
    c    = integer(num_doses),
    STFT = numeric(num_doses),
    n0   = numeric(num_doses)
  )

  if (nrow(df) == 0) return(summary)

  for (d in seq_len(num_doses)) {
    at_d <- df$dose == d
    if (!any(at_d)) next

    n_d <- sum(at_d)
    s_d <- sum(df$tox[at_d])
    r_d <- sum(df$weight[at_d] == 1)
    c_d <- n_d - r_d

    if (c_d > 0) {
      pending <- at_d & df$weight < 1
      STFT_d <- sum(df$weight[pending])
      n0_d <- r_d + STFT_d
    } else {
      STFT_d <- 0
      n0_d <- n_d
    }

    summary$n[d]    <- n_d
    summary$s[d]    <- s_d
    summary$r[d]    <- r_d
    summary$c[d]    <- c_d
    summary$STFT[d] <- STFT_d
    summary$n0[d]   <- n0_d
  }

  summary
}

# -- Constructor -------------------------------------------------------------

tite_mad_selector <- function(outcomes, num_doses, target,
                              decision_fn, decision_args,
                              design_label, ...) {

  # Parse outcomes
  if (is.character(outcomes)) {
    df <- escalation::parse_phase1_outcomes(outcomes, as_list = FALSE)
    if (nrow(df) > 0) df$weight <- 1
  } else if (is.data.frame(outcomes)) {
    df <- outcomes
    if (!is.integer(df$dose)) df$dose <- as.integer(df$dose)
    if (!is.integer(df$tox)) df$tox <- as.integer(df$tox)
    if ("cohort" %in% colnames(df) && !is.integer(df$cohort))
      df$cohort <- as.integer(df$cohort)
    if (!"weight" %in% colnames(df)) df$weight <- 1
  } else {
    stop("outcomes should be a character string or a data-frame.")
  }

  # Per-dose summary
  dose_summary <- .tite_dose_summary(df, num_doses)

  # Dose decision
  if (nrow(df) == 0) {
    rec_dose <- 1L
  } else {
    last_dose <- df$dose[nrow(df)]
    s_d <- dose_summary$s[last_dose]
    n0_d <- dose_summary$n0[last_dose]

    action <- decision_fn(s_d, n0_d, decision_args)

    if (action == "escalate") {
      rec_dose <- min(last_dose + 1L, num_doses)
    } else if (action == "deescalate") {
      rec_dose <- max(last_dose - 1L, 1L)
    } else if (action == "stay") {
      rec_dose <- last_dose
    } else {
      stop("decision_fn must return 'escalate', 'deescalate', or 'stay', got: '",
           action, "'")
    }
    rec_dose <- as.integer(rec_dose)
  }

  # Isotonic probability estimates via BOIN::select.mtd
  if (nrow(df) == 0) {
    isotonic_phat <- rep(NA_real_, num_doses)
  } else {
    mtd_result <- BOIN::select.mtd(
      target = target,
      npts = dose_summary$n,
      ntox = dose_summary$s
    )
    isotonic_phat <- suppressWarnings(as.numeric(mtd_result$p_est$phat))
  }

  l <- list(
    df = df,
    num_doses = num_doses,
    target = target,
    dose_summary = dose_summary,
    recommended_dose = rec_dose,
    isotonic_phat = isotonic_phat
  )
  class(l) <- c(paste0(design_label, "_selector"),
                 "tite_mad_selector",
                 "tox_selector",
                 "selector")
  l
}

# -- S3 methods on tite_mad_selector ----------------------------------------

#' @export
recommended_dose.tite_mad_selector <- function(x, ...) {
  x$recommended_dose
}

#' @export
continue.tite_mad_selector <- function(x, ...) {
  TRUE
}

#' @export
tox_target.tite_mad_selector <- function(x, ...) {
  x$target
}

#' @export
num_doses.tite_mad_selector <- function(x, ...) {
  x$num_doses
}

#' @export
num_patients.tite_mad_selector <- function(x, ...) {
  nrow(x$df)
}

#' @export
cohort.tite_mad_selector <- function(x, ...) {
  if (nrow(x$df) == 0) return(integer(0))
  x$df$cohort
}

#' @export
doses_given.tite_mad_selector <- function(x, ...) {
  if (nrow(x$df) == 0) return(integer(0))
  x$df$dose
}

#' @export
tox.tite_mad_selector <- function(x, ...) {
  if (nrow(x$df) == 0) return(integer(0))
  x$df$tox
}

#' @export
weight.tite_mad_selector <- function(x, ...) {
  if (nrow(x$df) == 0) return(numeric(0))
  x$df$weight
}

#' @export
tox_at_dose.tite_mad_selector <- function(x, ...) {
  x$dose_summary$s
}

#' @export
mean_prob_tox.tite_mad_selector <- function(x, ...) {
  x$isotonic_phat
}

#' @export
dose_admissible.tite_mad_selector <- function(x, ...) {
  rep(TRUE, x$num_doses)
}

#' @export
prob_tox_exceeds.tite_mad_selector <- function(x, threshold, ...) {
  n_at <- x$dose_summary$n
  s_at <- x$dose_summary$s
  result <- 1 - pbeta(threshold, 1 + s_at, 1 + n_at - s_at)
  result[n_at == 0] <- NA_real_
  result
}

#' @export
prob_tox_quantile.tite_mad_selector <- function(x, p, ...) {
  n_at <- x$dose_summary$n
  s_at <- x$dose_summary$s
  result <- qbeta(p, 1 + s_at, 1 + n_at - s_at)
  result[n_at == 0] <- NA_real_
  result
}

#' @export
median_prob_tox.tite_mad_selector <- function(x, ...) {
  prob_tox_quantile(x, 0.5)
}

#' @export
supports_sampling.tite_mad_selector <- function(x, ...) {
  FALSE
}
