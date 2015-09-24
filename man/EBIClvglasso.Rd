\name{EBIClvglasso}
\alias{EBIClvglasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Latent variable graphical LASSO using EBIC to select optimal tuning parameter
}
\description{
This function minimizes the Extended Bayesian Information Criterion (EBIC; Chen and Chen, 2008) to choose the lvglasso tuning parameter. See \code{\link{lvglasso}}
}
\usage{
EBIClvglasso(S, n, nLatents, gamma = 0.5, nRho = 100, lambda, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
Sample variance-covariance matrix
}
  \item{n}{
Sample Size
}
  \item{nLatents}{
Number of latent variables
}
  \item{gamma}{
EBIC hyper-parameter
}
  \item{nRho}{
Number of tuning parameters to test
}
  \item{lambda}{
The lambda argument containing factor loadings, only used for starting values!
}
  \item{\dots}{
Arguments sent to \code{\link{lvglasso}}
}
}

\value{
The optimal result of \code{\link{lvglasso}}, with two more elements:
\item{rho}{The selected tuning parameter}
\item{ebic}{The optimal EBIC}
}
\references{
Chen, J., & Chen, Z. (2008). Extended Bayesian information criteria for model selection with large model spaces. Biometrika, 95(3), 759-771.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}


\seealso{
\code{\link{lvglasso}}
}
