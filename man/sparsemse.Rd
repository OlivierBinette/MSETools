% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparsemse.R
\name{sparsemse}
\alias{sparsemse}
\title{Fit the SparseMSE model}
\usage{
sparsemse(data, pthresh=0.02)
}
\arguments{
\item{data}{object of class "MSEdata" representing list inclusion pattern counts.}

\item{pthresh}{p-value threshold for forward stepwise p-value thresholding. Default is 0.02.}
}
\value{
Object of class "sparsemse" containing the fitted model "fit" as well as the arguments "args" passed to this function. Use the function \code{estimates()} to recover point and interval estimates.
}
\description{
Fit the SparseMSE approach of Chan et al. (2020). Note: computationally-intensive bootstrap confidence intervals are only computed as part of a call to the \code{estimates()} function.
}
\examples{
sparsemse_fit <- sparsemse(UK)
estimates(sparsemse_fit, mc.cores = 1)

}
