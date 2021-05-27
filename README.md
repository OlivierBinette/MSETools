
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MSETools

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/MSETools)](https://CRAN.R-project.org/package=MSETools)
<!-- badges: end -->

**MSETools** provides a unified interface to multiple systems estimation
(MSE) software. It implements best usage practices and computational
speedups. Data from multiple system estimation studies of human
trafficking has been reproduced for illustrations and analyses.

The following sections of the README file provide a general overview of
**MSETools**, including some of the available functions, installation
instructions, and simple examples. For a more in-depth introduction,
consult the following vignettes and user manual:

-   [Computing prevalence estimates using **MSETools**]()
-   [Working with MSETools on a computer cluster]()
-   [**MSETools** user manual)]()

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("OlivierBinette/MSETools")
```

## Summary

### Implemented models

-   **dga:** The `dga()` function provides an interface to the `dga`
    package of Lum, Johndrow and Ball (2015), here re-implemented using
    Rcpp and extended to allow more flexible priors as well as Bayesian
    fractional posterior distributions. The package implements
    decomposable graphical models with hyper-Dirichlet priors and
    Bayesian model averaging.
-   **LCMCR:** The `lcmcr()` function provides an interface to the
    `LCMCR` package of Manrique-Vallier (2020), which implements the
    latent class model of Manrique-Vallier (2016). By default,
    **MSETools** initializes 500 parallel MCMC chains to provide
    cross-replication stability – this is necessary since the LCMCR
    Gibbs sampler fails to converge in some cases. Convergence
    diagnostics are available through `MSETools::diagnostics()`.
-   **SparseMSE:** The `sparsemse()` function provides an interface to
    the `SparseMSE` package of Chan, Silverman and Vincent (2019), which
    implements a Poisson log-linear approach with stepwise model
    selection and bootstrap confidence intervals. An option has been
    added to parallelize bootstrap replications across available cores.

### Helper functions

-   `estimates()` computes point estimates and confidence intervals for
    a list of models.
-   `batch.estimates()` compute estimates as a slurm job array (for use
    in a cluster).
-   `summary()` provides a result summary of fitted models.
-   `diagnostics()` provides convergence diagnostics for `lcmcr`
    objects.

### Datasets

-   **UK:**
-   **New Orleans:**
-   **Western U.S.:**
-   **Netherlands:**
-   **Australia:**

## Examples

``` r
library(MSETools)
```

Define a list of models fitted on the UK dataset:

``` r
models = list(lcmcr(UK), sparsemse(UK), dga(UK), independence(UK))
```

Compute estimates:

``` r
estimates(models)
```

Parallelize on a cluster:

``` r
batch.estimates(models, njobs=4)
```

## References
