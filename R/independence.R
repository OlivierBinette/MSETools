#' Fit an independence model
#'
#' Fit a (sparse) independence model to the data using the \code{SparseMSE::modelfit()} function.
#'
#' @usage independence(data)
#'
#' @param data object of class "MSEdata" representing list inclusion pattern counts.
#'
#' @return Object of class "independence" containing the fitted model "fit" and the arguments "args" passed to this function.
#'
#' @examples
#' indep_fit <- independence(UK)
#' estimates(indep_fit)
#' @import SparseMSE
#' @export
independence <- function(data) {
  assert(is.MSEdata(data))
  args <- as.list(environment())

  fit <- SparseMSE::modelfit(data, mX = NULL)

  independence_fit <- list(fit = fit, args = args)
  structure(independence_fit, class=c("MSEfit", "independence"))
}

is.independence <- function(fit) {
  inherits(fit, "independence")
}
