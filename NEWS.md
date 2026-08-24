# nlmixr2llp 0.2

* `runLLPControl()` gains an `rxThreads` option, forwarded to
  `nlmixr2utils::.withWorkerPlan()`, so parallel profiling (`workers > 1` or
  `workers = "auto"`) can avoid nlmixr2utils's rxode2-thread oversubscription
  guard instead of always aborting.

* Accommodated `nlmixr2est`'s full theta+sigma+Omega covariance: residual-error
  thetas (e.g. `add.sd`) are no longer misclassified as structural just
  because they're now present in `fit$cov` (asymmetry-warning threshold fix),
  and `llpOmegaSE()` now prefers nlmixr2est's reported `om.<eta>` SE over the
  Wishart chi-squared approximation when available.

* A re-fit that fails during profiling (e.g. from an upstream `nlmixr2est`/
  `rxode2` API break, not just ordinary non-convergence) now surfaces its
  original error message via a warning instead of silently producing an
  opaque all-NA profile row.

* The omega-diagonal profiling path (`buildFixedOmegaModel()`,
  `llpRunOneOmega()`, `runLLP()` with an omega name in `which`) now has real
  test coverage against an actual fit instead of being entirely skipped.

* `nlmixr2utils` is now required to be at least version 0.3, the version that
  introduced the `.withWorkerPlan()` `rxThreads` argument and
  `resolveRxThreads()` this package now relies on.

* `runLLP()`'s parallel dispatch now uses `nlmixr2utils::.plap()` instead of
  a hand-rolled `future.apply::future_lapply()` loop, picking up `progressr`
  progress reporting for free when that package is installed.

* `llpProfileFixed(method = "fixed")` now gives a clear error naming the
  unrecognized option(s) if `control` is non-empty, instead of R's raw
  "unused argument" error.

* `llpFixedControl()`, `llpProfileFixed()`, `llpProfileFixedSingle()`, and `runLLPControl()` now replace the old conflicting helper names, and the deprecated `profileLlp()` path has been removed.

* Initial package split from `nlmixr2extra`, providing `runLLP()`,
  `runLLPControl()`, `llpProfileFixed()`, control helpers, S3 methods, tests,
  and the LLP vignette as a standalone package depending on `nlmixr2utils`.
