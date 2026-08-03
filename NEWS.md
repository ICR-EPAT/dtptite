# dtptite 0.0.0.9000

* `fit.dtp_selector_factory()` now skips the wait projection when the parent's
  recommendation is already the top dose. `.dtp_project()` only reports a wait
  when the projected dose rises strictly above the recommendation, which is
  unsatisfiable there, so this is behaviour-preserving and saves one parent fit
  per occurrence. A de-escalation off the top dose still projects, since
  recovering to the top dose is a genuine hit.

* Fixed wait-heavy DTP trials terminating before enrolling their planned sample
  (#32). The depth safety valve counted recorded trajectory entries, of which a
  DTP trial produces one per cohort *and* one per wait, so a wait-heavy
  cohort-size-1 trial hit the cap (`max_i = 30`) and stopped early with a valid
  recommendation. `phase1_dtp_tite_sim()` now defaults `i_like_big_trials = TRUE`;
  termination is guaranteed by the design's stopping rules.

* Fixed a DTP simulation bug (#30) where a trial stopped by a DLT *during* a
  wait recorded the stale pre-wait fit as its terminal state, so
  `recommended_dose()` returned a valid dose instead of `NA` and
  `trial_duration()` returned the pre-wait (too-short) time. The post-wait fit
  is now recorded, so the terminal fit and time reflect the DLT stop.
* `apply_dtp()` gains `projection_search`. The default `"binary"` screens the
  top of the candidate wait range and then bisects, replacing the day-by-day
  scan and cutting DTP simulation time by 3.1x (`t_max = 14`) to 5.9x
  (`t_max = 56`) with bit-identical results. Pass `projection_search =
  "linear"` for the previous exhaustive scan.
* `apply_dtp()` gains `verbose`. The "t_max covers all remaining follow-up"
  message fired on nearly every projection during simulation; it is now
  silenced automatically for the duration of a simulation run and can be
  turned off explicitly with `verbose = FALSE`.
* Initial package skeleton with `DESCRIPTION`, `NAMESPACE`, and package-level documentation.
