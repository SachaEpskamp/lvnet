\name{lvglasso}
\alias{lvglasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Latent variable graphical LASSO
}
\description{
The lvglasso algorithm to estimate network structures containing latent variables, as proposed by Yuan (2012). Uses the glasso package (Friedman, Hastie and Tibshirani, 2014) and mimics input and output of the \code{\link[glasso]{glasso}} function.
}
\usage{
lvglasso(S, nLatents, rho = 0, thr = 1e-04, maxit = 10000, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
Sample variance-covariance matrix
}
  \item{nLatents}{
Number of latent variables.
}
  \item{rho}{
The LASSO tuning parameter
}
  \item{thr}{
The threshold to use for convergence
}
  \item{maxit}{
Maximum number of iterations
}
  \item{lambda}{
The lambda argument containing factor loadings, only used for starting values!}
}

\value{
A list of class \code{lvglasso} containing the following elements:
\item{w}{The estimated variance-covariance matrix of both observed and latent variables}
\item{wi}{The estimated inverse variance-covariance matrix of both observed and latent variables}
\item{pcor}{Estimated partial correlation matrix of both observed and latent variables}
\item{observed}{Logical vector indicating which elements of w, wi and pcor are observed}
\item{niter}{The number of iterations used}
\item{lambda}{The estimated lambda matrix, when result is transformed to EFA model}
\item{theta}{The estimated theta matrix}
\item{omega_theta}{The estimated omega_theta matrix}
\item{psi}{The estimated psi matrix}
}
\references{
Yuan, M.  (2012).  Discussion:  Latent variable graphical model selection via convex optimization.The Annals of Statistics,40, 1968-1972.

Jerome Friedman, Trevor Hastie and Rob Tibshirani (2014). glasso: Graphical lasso-estimation of Gaussian graphical models. R package version 1.8. http://CRAN.R-project.org/package=glasso
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
