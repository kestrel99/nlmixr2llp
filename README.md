# nlmixr2llp

`nlmixr2llp` provides log-likelihood profiling tools for `nlmixr2` population
PK/PD models.

The package focuses on PsN-style likelihood profiling workflows for `nlmixr2` models, 
including:

* `runLLP()` for parameter-wise log-likelihood profiling
* `runLLPControl()` for tuning convergence, initial guesses, and parallel workers
* `llpProfileFixed()` for estimating OFV at user-specified fixed parameter values
* S3 `print()`, `confint()`, and `plot()` methods for profile results

`nlmixr2llp` is designed to work alongside `nlmixr2utils`, which provides the
shared `nlmixr2()`/`ini()`/`model()` helpers and worker-plan infrastructure
used across the split `nlmixr2` extension packages.

## Installation

The package is not on CRAN. Install it together with `nlmixr2utils`.

Using `pak`:

```r
pak::pkg_install(c(
  "kestrel99/nlmixr2utils",
  "kestrel99/nlmixr2llp"
))
```

Using `remotes`:

```r
remotes::install_github("kestrel99/nlmixr2utils")
remotes::install_github("kestrel99/nlmixr2llp")
```

## Basic Use

```r
library(nlmixr2data)
library(nlmixr2utils)
library(nlmixr2llp)

one_cmt <- function() {
  ini({
    tka <- log(1.57)
    tcl <- log(2.72)
    tv <- fixed(log(31.5))
    eta.ka ~ 0.6
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl)
    v <- exp(tv)
    cp <- linCmt()
    cp ~ add(add.sd)
  })
}

fit <- nlmixr2(
  one_cmt,
  data = nlmixr2data::theo_sd,
  est = "focei",
  control = list(print = 0L)
)

prof <- runLLP(fit, which = c("tka", "tcl", "add.sd", "eta.ka"))
print(prof)
confint(prof)
plot(prof)
```

## How this differs from `nlmixr2extra`

The core statistical goal is the same as in `nlmixr2extra`: move one
parameter at a time away from the MLE, re-estimate the remaining parameters,
and locate the lower and upper parameter values that correspond to a target
increase in OFV. The main differences in `nlmixr2llp` are in the public
workflow, the supported parameter types, and the diagnostics returned with the
profile.

Compared with `profileLlp()` and `llpControl()` in `nlmixr2extra`,
`nlmixr2llp`:

* uses `runLLP()` and `runLLPControl()` as the primary public interface, with
  the legacy `profileLlp()` path removed and the fixed-evaluation helpers
  renamed to `llpProfileFixed()` and `llpFixedControl()`
* can profile supported OMEGA diagonal variance parameters in addition to the
  usual fixed-effect parameters; the old `nlmixr2extra` workflow was limited
  to parameters returned by `fixef(fit)`
* chooses initial profile steps using a more PsN-like rule:
  `normq * SE` when a standard error is available, with
  `normq * rseTheta/100 * |MLE|` as fallback, instead of relying only on the
  older RSE-based width heuristic
* handles parameter bounds more carefully, using soft clamping of initial
  guesses and explicit boundary tracking, which is especially important for
  asymmetric profiles and variance parameters near zero
* supports parameter-wise parallel execution through `workers`, with isolated
  temporary fit directories to avoid file collisions between simultaneous
  profiling jobs
* returns an `nlmixr2Profile` object rather than only a plain `data.frame`;
  the row data remain profile-table compatible, but the object now also carries
  status metadata used by `print()`, `confint()`, and `plot()`
* records extra profile diagnostics such as whether a bound was approached,
  whether a direction hit `itermax`, and the lower/upper interval ratio, and
  warns on strongly asymmetric profiles

One practical consequence is that the current package is better suited for
production profiling runs where you need to inspect convergence behavior, work
with omega variance profiles, or spread multiple parameter profiles across
workers without manually orchestrating them.

## References

The implementation follows the [PsN LLP workflow](https://github.com/UUPharmacometrics/PsN/releases/download/v5.7.0/llp_userguide.pdf) and related likelihood
profiling literature, adapted for `nlmixr2` model fits.

For a worked example, see:
`vignette("runLLP", package = "nlmixr2llp")`.
