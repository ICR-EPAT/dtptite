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
* **Behaviour change.** Fixed the simulation clock (#32), which conflated the
  accrual process with the trial timeline. `time_now` was a single cursor
  serving as both "when did the last patient arrive" and "when is the model
  being fit", and arrivals were generated as offsets from it, so a minimum
  follow-up period or a DTP wait silently pushed recruitment. With cohort size
  1, 14-day accrual and `min_fup_time = 14`, patients were recruited at
  t = 14, 42, 70, 98 — a 28-day rhythm.

  The pool of available patients now runs on its own cursor, untouched by
  anything on the trial timeline. A cohort opens when its first patient
  arrives and the model is fit *then*, before anyone in that cohort is dosed,
  so the recommendation applies to the patients actually in front of it. DTP is
  an add-on to that update: a wait moves the earliest time a cohort may be
  dosed, never the time the model updates.

  Consequences for existing code:

  - `min_fup_time` is now an inter-cohort **stagger floor** rather than a delay
    added to accrual. Raising it no longer stretches the inter-arrival rhythm.
    The default of 0 leaves model updates purely accrual-driven, which is what
    most callers want; scripts passing `min_fup_time = accrual_gap` to
    compensate for the old behaviour should drop it.
  - Patients who present while the waiting room is full (capacity
    `queue_size`, default the cohort size) are no longer enrolled. This is
    reported in the new `missed` column of `dtp_wait_events()` and
    `num_missed` in `dtp_wait_summary()`. A patient is missed exactly when
    `wait_duration > queue_size * accrual_gap`.
  - Cohort members need not share a dosing time: only patients held through a
    wait are dosed at its end, the rest on arrival.
  - Non-DTP arms dispatch to `escalation::phase1_tite_sim`, which retains the
    old clock. Build baselines as `apply_dtp(t_max = 0, obswin)` so both arms
    run the same simulator; `misc/sensitivity/run_sweep.R` now does this.
  - **Stored operating characteristics are not comparable across this
    change.** Selection is broadly stable (PCS within ~2 points in testing),
    but trial durations fall substantially — roughly 1.7x at cohort size 1 and
    1.3x at cohort size 3 — and wait activity moves in opposite directions
    depending on cohort size, so it must be re-measured rather than
    extrapolated. Re-run `misc/sensitivity`.

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
