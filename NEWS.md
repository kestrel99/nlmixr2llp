# nlmixr2llp 0.1

* `llpFixedControl()`, `llpProfileFixed()`, `llpProfileFixedSingle()`, and `runLLPControl()` now replace the old conflicting helper names, and the deprecated `profileLlp()` path has been removed.

* Initial package split from `nlmixr2extra`, providing `runLLP()`,
  `runLLPControl()`, `llpProfileFixed()`, control helpers, S3 methods, tests,
  and the LLP vignette as a standalone package depending on `nlmixr2utils`.
