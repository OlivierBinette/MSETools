
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MSETools

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**MSETools** provides a unified interface to multiple systems estimation
(MSE) software. It implements best usage practices and computational
speedups. Data from multiple system estimation studies of human
trafficking has been reproduced for illustrations and analyses. Analyses
from Binette and Steorts (2021) are contained in the “analyses” folder.

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

See Binette and Steorts (2021) for a description of the datasets
reproduced herein.

## Examples

``` r
library(MSETools)
```

Define a list of models fitted to the UK dataset:

``` r
models = list(lcmcr(UK), sparsemse(UK), dga(UK), independence(UK))
```

Compute estimates:

``` r
estimates(models)
#> [[1]]
#>    N.hat    N.lwr    N.upr 
#> 20232.00 10934.95 32842.27 
#> 
#> [[2]]
#>     N.hat     N.lwr     N.upr 
#> 11312.990  9185.358 15610.886 
#> 
#> [[3]]
#> N.hat N.lwr N.upr 
#> 23016 10408 33481 
#> 
#> [[4]]
#>    N.hat    N.lwr    N.upr 
#> 13444.13 12004.03 15160.33
```

Parallelize on a computing cluster:

``` r
batch.estimates(models, njobs=4)
```

## References

-   Binette, O. and Steorts, Rebecca C. (2021) On the Reliability of
    Multiple Systems Estimation for the Quantification of Modern
    Slavery.
-   Lax Chan, Bernard Silverman and Kyle Vincent (2019). SparseMSE:
    ‘Multiple Systems Estimation for Sparse Capture Data’. R package
    version 2.0.1. <https://CRAN.R-project.org/package=SparseMSE>
-   James Johndrow, Kristian Lum and Patrick Ball (2021). dga:
    Capture-Recapture Estimation using Bayesian Model Averaging. R
    package version 2.0.1. <https://CRAN.R-project.org/package=dga>
-   Daniel Manrique-Vallier (2020). LCMCR: Bayesian Non-Parametric
    Latent-Class Capture-Recapture. R package version 0.4.11.
    <https://CRAN.R-project.org/package=LCMCR>
