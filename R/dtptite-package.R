#' @section Funding:
#' Xiaoran Lai, NIHR Development and Skills Enhancement Award (DSE),
#' NIHR500663, is funded by the NIHR for this research project. The views
#' expressed in this publication are those of the author(s) and not
#' necessarily those of the NIHR, NHS or the UK Department of Health and
#' Social Care.
#'
#' @keywords internal
#' @importFrom escalation recommended_dose fit doses_given continue
#' @importFrom escalation num_patients mean_prob_tox weight num_doses tox_target
#' @importFrom escalation dose_admissible supports_sampling
#' @importFrom escalation n_at_dose tox_at_dose
#' @importFrom escalation simulation_function
#' @importFrom BOIN get.boundary select.mtd
#' @importFrom escalation tox cohort prob_tox_quantile median_prob_tox
#' @importFrom escalation prob_tox_exceeds empiric_tox_rate
#' @importFrom escalation parse_phase1_outcomes trial_duration
#' @importFrom dplyr bind_rows group_by summarise left_join n
#' @importFrom tibble as_tibble tibble
#' @importFrom stats pbeta qbeta runif quantile median qnorm
"_PACKAGE"

# Local definition rather than base R's, which only exists from R 4.4.
`%||%` <- function(x, y) if (is.null(x)) y else x

