# Plan — `feature/time-to-first-mtd`

## Metric

Time until the **first patient is dosed at the true MTD**. Applies to any design —
DTP or not, TITE or not.

## Architecture

Computed post-hoc from the **final fit's patient-level data frame**, which is
cumulative and carries `dose` + `time`. No change to `simulate.R`, no
`return_all_fits = TRUE`, no memory cost. Verified across DTP+TITE-CRM,
TITE-CRM, plain CRM, DTP+TITE-BOIN, and plain BOIN.

**Time semantics:** the `time` column is *dosing* time. For non-DTP designs that
equals arrival; for DTP-queued patients it is the wait-end instant
(`simulate.R:412`, `:424` overwrite queued arrivals with `time_now`), so DTP's
queue delay is fully charged. Origin is trial start, `t = 0`.

## Surface — four exports, one internal

### `scenario_mtd(true_prob_tox, target, tol = 0.05)`

Tolerance band, `|p - target| <= tol` inclusive. Multiple qualifying doses →
**lowest**. Empty band → `NA`, with a warning distinguishing all-doses-toxic from
no-dose-within-tolerance. Target retrieved via `tox_target()`.

### `time_to_first_mtd(x, mtd = NULL, tol = 0.05)`

S3 on `simulations` / `simulations_collection`. Returns a classed tibble, one row
per replicate, with `num_sims` stamped as an attribute:

```
 [design]  replicate  mtd  time_to_mtd  reached  n_before_mtd
```

`mtd = NULL` derives; an integer overrides; `NA` asserts no MTD. First patient
(`min`), not full cohort.

### `summary.time_to_first_mtd(object, probs = NULL)`

```
 design num_sims mtd num_reached reached_fraction median_time_to_mtd median_n_before_mtd
    DTP      200   3         191            0.955              11.20                   6
 no-DTP      200   3          55            0.275                 NA                  12
```

`probs` appends `q10_time_to_mtd`, `q25_time_to_mtd`, … `median_time_to_mtd` is
permanent; `probs` adds to it rather than replacing it. Reads the denominator
from the `num_sims` attribute, warning if `nrow()` disagrees — this stops
`filter(reached) |> summary()` silently reporting `reached_fraction = 1`.

### `plot.time_to_first_mtd(x, axis = c("time", "patients"), ci = FALSE)`

Cumulative-incidence step curve, one per design, plateauing at the reach rate.
Calendar time by default. Optional pointwise Wilson intervals for Monte Carlo
error, off by default. ggplot2 guarded by `requireNamespace()`.

### `.patient_dose_times(fit)` (internal)

Recursive `$parent` walk following `.find_dtp_factory()`'s idiom
(`simulate.R:77`). Errors loudly rather than returning `NULL`, since
`model_frame()` drops `time` and this is the one point of contact with
escalation internals.

## Censoring

Empirical cumulative incidence over all replicates. Non-reachers stay in the
denominator; no-MTD scenarios are excluded entirely (`NA`, not `FALSE`).

**No Kaplan–Meier** — the trial horizon is a design feature, not a nuisance, and
early toxicity stops are competing events with no counterfactual continuation.

Quantiles reuse `stats::quantile()` by coding non-reachers as `+Inf` and mapping
`Inf` back to `NA`, so any `p` beyond the reach fraction is correctly undefined.
This produces the `NA` median for no-DTP at 27.5% reach, while `q10` and `q25`
remain estimable at 11.16 and 16.20 versus DTP's 9.59 and 10.32.

## Reuse

- `bind_rows(lapply(x, f), .id = "design")` verbatim from
  `dtp_wait_events.simulations_collection`
- The `x$fits[[i]][[length(...)]]` last-fit idiom
- `num_*` / `*_fraction` naming, with `reached_fraction` parallel to
  `wait_fraction`
- `.wait_events_empty()`'s zero-row pattern
- Existing NAMESPACE imports only — no new Imports
- Test fixtures `skeleton`, `target`, `arrivals3`, and the all-toxic scenario
  from `test-dtp-simulation-metrics.R`

## Dropped

`reason` factor, `stop_origin`, overshoot detection, the
`criterion`/`basis`/`phase` arguments, unequal-`n` warning, and a standalone
`time_to_first_mtd_summary()`.

## Docs

Roxygen states only what the function computes — event definition, band rule and
`tol`, time origin, `NA` meaning not reached by trial end, denominator
composition, column meanings. No interpretation guidance.

## Tests

| | |
|---|---|
| (a) | `.patient_dose_times()` across all five design classes + loud failure — the escalation-upgrade regression guard |
| (b) | `scenario_mtd()` exhaustively, no simulation: lowest-wins, all-toxic, all-safe, non-monotone, `== tol` boundary inclusive, invalid overrides |
| (c) | First-hit logic on synthetic frames, no simulation: single hit, repeated hits take `min`, never dosed, DTP shared-timestamp ties |
| (d) | Censoring arithmetic: `nrow == num_sims`, `reached_fraction` matches, plateau matches, median is `NA` below 50% reach |
| (e) | No-MTD → `mtd = NA`, `reached = NA`, excluded from denominator |
| (g) | Collection shape and `design` column |
| (h) | `axis = "patients"` consistency; `plot()` smoke test under `skip_if_not_installed("ggplot2")` |

Fixed seeds, `num_sims = 8` matching the existing suite, except (d) which needs
more replicates with a tolerance rather than exact figures.

## Files

- `R/time_to_first_mtd.R`
- `tests/testthat/test-time-to-first-mtd.R`
- roxygen → `NAMESPACE` + `man/`

Vignette coverage is out of scope here — it belongs with the dedicated
documentation issue.
