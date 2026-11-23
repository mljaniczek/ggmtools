## Submission

This is a new submission.

ggmtools collects generics, evaluation metrics and utilities shared by
packages that estimate multiple Gaussian graphical models ('spikeyglass'
and 'multiGGM', currently on GitHub). It depends only on base packages
('graphics', 'stats').

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Test environments

* local macOS (aarch64), R 4.5.2
<!-- TODO before submitting: run these and update this list.
* GitHub Actions: macOS, Windows and Ubuntu (R release), Ubuntu (R devel
  and R oldrel-1)
* win-builder (R devel)
-->

The 1 note is the new-submission note. (Locally there were also notes
about an unreachable time server and an old HTML Tidy, which come from
the check machine, not the package.)

## Notes for reviewers

* The `brca_pi3k` dataset is derived from open-access TCGA data via the
  Bioconductor package curatedTCGAData; the build script is in `data-raw/`
  of the GitHub repository and the source is documented in `?brca_pi3k`.
* `pcor_estimate.jgl()` is an S3 method for objects of class "jgl"
  produced by the JGL package, which was archived on CRAN in July 2025.
  The method only reads list elements of the fit, so ggmtools neither
  imports nor suggests JGL, and its examples use a small mock object.
