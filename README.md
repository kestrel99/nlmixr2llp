# nlmixr2llp

`nlmixr2llp` provides log-likelihood profiling tools for `nlmixr2` population
PK/PD models.

The package focuses on PsN-style likelihood profiling workflows for FOCE-I
fits, including:

* `runLLP()` for parameter-wise log-likelihood profiling
* `llpControl()` for tuning convergence, initial guesses, and parallel workers
* `profileFixed()` for estimating OFV at user-specified fixed parameter values
* S3 `print()`, `confint()`, and `plot()` methods for profile results

`nlmixr2llp` is designed to work alongside `nlmixr2utils`, which provides the
shared `nlmixr2()`/`ini()`/`model()` helpers and worker-plan infrastructure
used across the split `nlmixr2` extension packages.

## Installation

The package is not on CRAN. Install it together with `nlmixr2utils`.

Using `pak`:

```r
pak::pkg_install(c(
  "nlmixr2/nlmixr2utils",
  "nlmixr2/nlmixr2llp"
))
```

Using `remotes`:

```r
remotes::install_github("nlmixr2/nlmixr2utils")
remotes::install_github("nlmixr2/nlmixr2llp")
```

If you are working locally with the split repositories:

```r
devtools::install_local("../nlmixr2utils")
devtools::install_local("../nlmixr2llp")
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

## References

The implementation follows the PsN LLP workflow and related likelihood
profiling literature, adapted for `nlmixr2` FOCE-I fits.

For a worked example, see:
`vignette("runLLP", package = "nlmixr2llp")`.
