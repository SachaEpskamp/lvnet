\name{lvnetLasso}
\alias{lvnetLasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
LASSO model selection
}
\description{
This function runs \code{lvnet} for a number of different tuning parameters, selects the best model based on some criterion and refits that model to obtain accurate parameter estimates. The \code{\link{lassoSelect}} function can afterwards be used to select a different model.
}
\usage{
lvnetLasso(data, lassoMatrix, lassoTol = 1e-04, nTuning = 20, tuning.min = 0.01, 
    tuning.max = 0.5, criterion = c("bic", "aic", "ebic"), verbose = TRUE, refit = TRUE, 
    nCores = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The data argument as used in \code{\link{lvnet}}
}
  \item{lassoMatrix}{
Vector indicating the matrix or matrices to use in LASSO optmimization
}
  \item{lassoTol}{
Tolerance for absolute values to be treated as zero in counting parameters.
}
  \item{nTuning}{
Number of tuning parameters to estimate. 
}
  \item{tuning.min}{
Minimal tuning parameter
}
  \item{tuning.max}{
Maximal tuning parameter
}
  \item{criterion}{
Criterion to use in model selection
}
  \item{verbose}{
Should progress be printed to the console?
}
  \item{refit}{
Logical, should the best fitting model be refitted without LASSO regularization?
}
  \item{nCores}{
Number of cores to use in parallel computing.
}
  \item{\dots}{
Arguments sent to \code{\link{lvnet}}
}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
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

# Search best fitting omega_theta:
\dontrun{
res <- lvnetLasso(Data, "omega_theta", lambda = Lambda)
res$best
summary(res)
}
}
