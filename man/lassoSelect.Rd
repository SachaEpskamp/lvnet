\name{lassoSelect}
\alias{lassoSelect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Update lvnatLasso results to select a different model
}
\description{
This function can be used to select a model using any fit index
}
\usage{
lassoSelect(object, select, minimize = TRUE, refit = TRUE, lassoTol = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An \code{lvnetLasso} object
}
  \item{select}{
A raw R expression using names used in the \code{object$fitMeasures} part of the output of \code{lvnet}
}
  \item{minimize}{
Logical. Minimize or maximize?
}
  \item{refit}{
Logical. Should the new best model be refitted.
}
  \item{lassoTol}{
Tolerance for absolute values to be treated as zero in counting parameters.
}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
\examples{
\dontrun{
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
res <- lvnetLasso(Data, "omega_theta", lambda = Lambda)
res$best
summary(res)

# Update to use EBIC:
resEBIC <- lassoSelect(res, ebic)
summary(resEBIC)

# Update to use minimal fitting model with RMSEA < 0.05:
resMinimal <- lassoSelect(res, df * (rmsea < 0.05), minimize = FALSE)
summary(resMinimal)
}
}
