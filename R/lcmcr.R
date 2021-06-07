#' Fit a latent class model
#'
#' Fit the latent class model of Manrique-Vallier (2016). This function initializes the MCMC chains and sampling is postponed until a call to the \code{estimates()} function. By default, 200 independent chains are initialized.
#'
#' @usage lcmcr(data, K=10, seeds=1:200, lazy=TRUE, ...)
#'
#' @param data object of class "MSEdata" representing list inclusion pattern counts.
#' @param K maximum number of latent classes.
#' @param seeds numeric vector of RNG seeds for the MCMC chains. A chain is initialized for vector element.
#' @param lazy whether or not to wait before initializing the chains with lcmcr objects. Default is TRUE.
#' @param ... other parameters passed to the \code{LCMCR::lcmCR()} function for each chain.
#'
#' @value Object of class "lcmcr" containing the list of chains "lcmcr_chains", and the arguments "args" passed to this function.
#'
#' @examples
#' lcmcr_fit <- lcmcr(UK)
#' estimates(lcmcr_fit, mc.cores = 1)
#' @import assert LCMCR
#' @export
lcmcr <- function(data, K = 10, seeds = 1:200, lazy=TRUE, ...) {
  assert(is.MSEdata(data))
  args <- list(data=data, K=K, seeds=seeds, lazy=lazy, ...)

  patterns <- lapply(data[colnames(data) != "count"], factor)
  df <- cbind(as.data.frame(patterns), Freq = data[, "count"])

  if (lazy) {
    lcmcr_chains <- NULL
  } else {
    lcmcr_chains <- lapply(seeds, function(seed) {
      lcmcr <- lcmCR(df, tabular = TRUE, K = K, seed = seed, buffer_size=1, ...)

      # Trace initialization
      lcmcr$Set_Trace("prob_zero")
      lcmcr$Set_Trace("n0")
      lcmcr$Set_Trace("k_star")
      #lcmcr$Set_Trace("alpha")
      #lcmcr$Set_Trace("log_nuK")

      return(lcmcr)
    })
  }

  lcmcr_fit <- list(
    lcmcr_chains = lcmcr_chains,
    args = args
  )
  structure(lcmcr_fit, class=c("lcmcr", "MSEfit"))
}

is.lcmcr <- function(obj) {
  inherits(obj, "lcmcr")
}


#' Get MCMC traces for LCMCR
#'
#' Run MCMC and return traces for the "n0", "prob_zero" and "k_star" parameters of LCMCR objects. For each of these parameters, a matrix is returned where each row represents an iteration and each column a chain.
#'
#' @usage MCMCtraces(fit, burnin = 0, nSamples = 100, thinning = 1000, mc.cores=detectCores())
#'
#' @param fit object of class "lcmcr".
#' @param burnin number of burn-in iterations.
#' @param nSamples number of samples to be returned.
#' @param thinning number of iterations between each sample.
#' @param mc.cores number of cores to use for computation.
#'
#' @import purrr dplyr parallel
#' @export
MCMCtrace <- function(fit, burnin = 0, nSamples = 100, thinning = 1000, mc.cores=detectCores()) {
  assert(is.lcmcr(fit))

  # Always re-run.
  fit <- lcmcr(fit$args$data, K = fit$args$K, seeds = fit$args$seeds, lazy=FALSE)
  traces <- mclapply(fit$lcmcr_chains, function(lcmcr) {
    lcmcr$Update(burnin, output = FALSE)
    lcmcr$Change_Trace_Length(nSamples)
    lcmcr$Change_SubSamp(thinning)
    lcmcr$Activate_Tracing()

    # Sample
    lcmcr$Update(nSamples * thinning, output = FALSE)

    trace = list(n0 = lcmcr$Get_Trace("n0"),
                 prob_zero = lcmcr$Get_Trace("prob_zero"),
                 k_star = lcmcr$Get_Trace("k_star"))
    return(trace)
  }, mc.cores = mc.cores, mc.preschedule = FALSE)

  n0 = sapply(map(traces, "n0"), function(x) x)
  prob_zero = sapply(map(traces, "prob_zero"), function(x) x)
  k_star = sapply(map(traces, "k_star"), function(x) x)

  return(list(n0 = n0,
              prob_zero = prob_zero,
              k_star = k_star))
}

#' Compute MCMC convergence diagnostics for MCMC trace
#'
#' Computes MCMC convergence diagnostics, R hat, effective sample size, and number of iterations, for MCMC trace obtained from the `MCMCtrace()` function.
#'
#' @usage diagnostics(traces)
#'
#' @param traces output from the `MCMCtrace()` function.
#'
#' @importFrom rstan Rhat ess_bulk
#' @export
diagnostics <- function(traces) {

  functs = list(Rhat = rstan::Rhat,
                n_ess = rstan::ess_bulk,
                n_iter = function(x) length(x))

  lapply(functs, function(FUN) {
    sapply(traces, function(trace) FUN(trace))
  }) %>%
    data.frame()
}

