## Resubmission

This is a resubmission. In response to the CRAN reviewer comments:

* The `predict.zinb_gp_fit()` example is no longer wrapped in `\dontrun{}`.
  It now fits a tiny model and predicts at new locations, and runs in well
  under five seconds. No examples use `\dontrun{}` or `\donttest{}`.
* The package no longer writes to the console with `print()` or `cat()`.
  MCMC progress reports are off by default and are emitted with `message()`
  when `print_progress = TRUE`, so they can be suppressed with
  `suppressMessages()`. Requesting a model with no Gaussian process now
  signals an error with `stop()` instead of printing text and returning
  `NULL`.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Test environments

* Local Windows 11, R 4.4.1: 0 errors, 0 warnings, 2 notes. The additional
  note stated that the check could not verify the current time.
* macOS 26.6, Apple Silicon, R 4.6.1 Patched: 0 errors, 0 warnings, 0 additional notes.
  The PDF manual was built successfully.
* win-builder R-release: 0 errors, 0 warnings, 0 additional notes.
* win-builder R-devel: 0 errors, 0 warnings, 0 additional notes.

## Reverse dependencies

There are currently no downstream dependencies for this package.
