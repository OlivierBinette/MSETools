
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MSETools

<!-- badges: start -->
<!-- badges: end -->

This repository contains all of the code and analyses for the paper
titled [“On the Reliability of Multiple Systems Estimation for the
Quantification of Modern Slavery” (Binette and Steorts,
2021)](https://arxiv.org/abs/2112.01594). It is structured as follows:

-   The R package **MSETools** provides a unified interface to multiple
    systems estimation (MSE) software. It implements best usage
    practices and computational speedups. Data from multiple system
    estimation studies of human trafficking has been reproduced for
    illustrations and analyses.

-   The **analyses** folder contains the analyses and figures for
    Binette and Steorts (2021). Each analysis is provided as an Rmd
    document which can be knitted on any platform using the included
    cache. Figures are saved to png and pdf format into subfolders.
    Cache can be regenerated by knitting the Rmd documents on computing
    clusters using SLURM. To run the entire analysis from scratch, make
    sure that the **MSETools** package is installed (see instructions
    below) and use the command
    `cd analyses && make clear_cache && srun make`. Long-running
    programs on a cluster can be run within a detachable terminal
    (e.g. `tmux`) to avoid connection issues.

## Installation

You can install the development version of **MSETools** from
[GitHub](https://github.com/):

``` r
# install.packages("devtools")
devtools::install_github("OlivierBinette/MSETools")
```

## Summary of MSETools

### Implemented models

-   **dga:** The `dga()` function provides an interface to the `dga`
    package of Lum, Johndrow and Ball (2015). The package implements
    decomposable graphical models with hyper-Dirichlet priors and
    Bayesian model averaging. Here it has been re-implemented in Rcpp
    and extended to allow more flexible prior distributions.
-   **LCMCR:** The `lcmcr()` function provides an interface to the
    `LCMCR` package of Manrique-Vallier (2020), which implements the
    latent class model of Manrique-Vallier (2016). By default,
    **MSETools** initializes 200 parallel MCMC chains to provide
    cross-replication stability – this is necessary since the LCMCR
    Gibbs sampler fails to converge in some cases. Convergence
    diagnostics are available through `MSETools::diagnostics()`.
-   **SparseMSE:** The `sparsemse()` function provides an interface to
    the `SparseMSE` package of Chan, Silverman and Vincent (2019), which
    implements a Poisson log-linear approach with stepwise model
    selection and bootstrap confidence intervals.

### Helper functions

-   `estimates()` computes point estimates and confidence intervals for
    a list of models.
-   `batch.estimates()` compute estimates as a SLURM job array for use
    in a cluster.
-   `diagnostics()` computes convergence diagnostics for `lcmcr`
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
```

Parallelize the computation of estimates on a computing cluster:

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
