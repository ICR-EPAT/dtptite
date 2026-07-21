# dtptite 0.0.0.9000

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
