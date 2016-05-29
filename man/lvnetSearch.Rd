\name{lvnetSearch}
\alias{lvnetSearch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Step-wise exploratory search for optimal fitting model
}
\description{
Performs stepwise search to optimize the structure of omega_theta, omega_psi, theta or psi. Starts at empty or full structure and iteratively adds or removes edges to optimize the criterion.
}
\usage{
lvnetSearch(data, matrix = c("omega_theta", "omega_psi", "theta", "psi"), 
          criterion = c("bic", "ebic","chisq","aic"), 
          start =  c("default","empty","full"),  alpha = 0.05, lambda, sampleSize, 
          maxIter,  nCores = 1, maxChange = 1, ...,  verbose = TRUE, file, 
          startValues = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
The data argument as used in \code{\link{lvnet}}
}
  \item{matrix}{
Character string indicating the matrix to be optimized. Can be \code{"omega_theta"}, \code{"omega_psi"}, \code{"theta"} and \code{"psi"}.
}
  \item{criterion}{
Character string indicating the criterion to be used. \code{"AIC"} and \code{"BIC"} optimize the AIC or BIC respectively, and \code{"chisq"} performs chi-square tests to see if adding an edge significantly improves model fit or removing an edges does not significantly reduce model fit.
}
  \item{start}{
A character string indicating the structure of the matrix at the start of the algorithm. \code{"empty"} starts with a matrix with only zeroes and \code{"full"} starts with a matrix in which all elements are free to estimate. \code{"lvglasso"} employs the lvglasso algorithm (\code{\link{EBIClvglasso}} to find a starting structure for omega_theta and \code{"glasso"} employs the glasso algorithm to find a starting point for omega_psi (\code{\link[qgraph]{EBICglasso}}). \code{"default"} will lead to a full matrix if omega_psi or psi is optimized, and an empty matrix if omega_theta or theta is optimized.
}
  \item{alpha}{
The alpha level for chi-square significance testing.
}
  \item{lambda}{
The lambda argument as used in \code{\link{lvnet}}
}
  \item{sampleSize}{
The sample size, only used if \code{data} is a covariance matrix.
}
  \item{maxIter}{
The maximum number of edges to test. Defaults to M(M-1)/2
}
\item{nCores}{
Number of cores to use in parallel estimation.
}
\item{maxChange}{
Set to higher than one to change multiple edges in each run. Each iteration, maxChange is reset to max(number of changed edges - 1, 1). Can result in instable results when searching \code{"omega_theta"}.
}
  \item{\dots}{
Arguments sent to \code{\link{lvnet}}
}
  \item{verbose}{
Logical if progress should be printed to the consile.
}
  \item{file}{
An optional character string containing a file name to store temporary results in.
}
  \item{startValues}{
A list containing start values as used in \code{\link{lvnet}}
}
}
\value{
An object of class \code{lvnetSearch}, which is a list containing:
\item{best}{The \code{lvnet} object of the best fitting model}
\item{modList}{A list containing the chain of fitted models}
\item{niter}{The number of iterations used}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

\seealso{
\code{\link{lvnet}}
}
\examples{
# Load dataset:
library("lavaan")
data(HolzingerSwineford1939)
Data <- HolzingerSwineford1939[,7:15]

# Measurement model:
Lambda <- matrix(0, 9, 3)
Lambda[1:3,1] <- NA
Lambda[4:6,2] <- NA
Lambda[7:9,3] <- NA

# Search best fitting omega_psi:
\dontrun{
res <- lvnetSearch(Data, "omega_psi", lambda = Lambda)
res$best
}
}

