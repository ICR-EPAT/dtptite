## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Test environments

* local: Debian 12 (WSL2, conda environment), R 4.5.2
* GitHub Actions: ubuntu-latest (R devel, release, oldrel-1), macos-latest (R release),
  windows-latest (R release)

## Notes

* Three reference-comparison tests that run full simulations against 'TITEgBOIN'
  take about 6 minutes, so they are skipped on CRAN with `skip_on_cran()`. They run
  in CI.
