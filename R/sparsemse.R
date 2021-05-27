#' Fit the SparseMSE model
#'
#' Fit the SparseMSE approach of Chan et al. (2020). Note: computationally-intensive bootstrap confidence intervals are only computed as part of a call to the \code{estimates()} function.
#'
#' @usage sparsemse(data, pthresh=0.02)
#'
#' @param data object of class "MSEdata" representing list inclusion pattern counts.
#' @param pthresh p-value threshold for forward stepwise p-value thresholding. Default is 0.02.
#'
#' @return Object of class "sparsemse" containing the fitted model "fit" as well as the arguments "args" passed to this function. Use the function \code{estimates()} to recover point and interval estimates.
#'
#' @examples
#' sparsemse_fit <- sparsemse(UK)
#' estimates(sparsemse_fit)
#' @import assert
#' @export
sparsemse <- function(data,
                      pthresh = 0.02) {
  assert(is.MSEdata(data))
  assert(is.numeric(pthresh), length(pthresh) == 1)

  args <- list(data=data, pthresh=pthresh)

  fit <- SparseMSE::estimatepopulation.0(data, pthresh = pthresh, quantiles = c())

  sparsemse_fit <- list(
    fit = fit,
    args = args
  )
  structure(sparsemse_fit, class=c("MSEfit", "sparsemse"))
}

is.sparsemse <- function(x) {
  inherits(x, "sparsemse")
}

getGLM <- function(sparsemse_fit) {
  assert(is.sparsemse(sparsemse_fit))

  return(sparsemse_fit$fit$MSEfit$fit)
}
