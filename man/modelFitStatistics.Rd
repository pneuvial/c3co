% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelFitStats.R
\name{modelFitStatistics}
\alias{modelFitStatistics}
\title{Print statistics of a c3co model}
\usage{
modelFitStatistics(Y, Yhat, What, Zhatt)
}
\arguments{
\item{Y}{Original n-by-J matrix of data for n samples and J segments.}

\item{Yhat}{Estimated n-by-J matrix for n samples and J segments.}

\item{What}{Inferred n-by-K weight matrix for n samples and K subclones.}

\item{Zhatt}{The transposed version of the inferred n-by-J subclones
matrix 'Zhat' for n samples and J segments.}
}
\value{
A named numeric vector with the main statistics of the inferred model:
\itemize{
\item \code{BIC}   : Bayesian Information Criterion
\item \code{PVE}   : Percentage of Variation Explained (Nowak et al., 2011)
\item \code{logLik}: Log Likelihood
\item \code{loss}  : Loss, which equals \eqn{RSS / (n * J)},
where RSS is the residual sum of squares
}
}
\description{
Print statistics of a c3co model
}
\references{
Nowak, G., Hastie, T., Pollack, J. R., & Tibshirani, R. (2011).
A fused lasso latent feature model for analyzing multi-sample
aCGH data. Biostatistics, 12(4), 776-791
}
