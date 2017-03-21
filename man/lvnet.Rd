\name{lvnet}
\alias{lvnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confirmatory Latent Variable Network Models
}
\description{
This function utilizes OpenMx (Boker et al., 2011, 2014) to confirmatory test latent variable network models between P manifests and M latents. See the details section for information about the modeling framework used. All the input matrices can be assigned R matrices with numbers indicating fixed values and NA indicating a value is free to estimate.
}
\usage{
lvnet(data, lambda, beta, omega_theta, delta_theta, omega_psi, delta_psi, psi, theta, 
        sampleSize, fitInd, fitSat, startValues = list(), scale = FALSE, nLatents,
        lasso = 0,  lassoMatrix, lassoTol = 1e-4, ebicTuning = 0.5, 
        mimic = c("lavaan","lvnet")
        )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
An N (sample size) x P matrix or data frame containing the raw data, or a P x P variance-covariance matrix.
}
  \item{lambda}{
A P x M matrix indicating factor loadings. Defaults to a full NA P x M matrix if psi or omega_psi is not missing, or a P x 0 dummy matrix.
}
  \item{beta}{
An M x M  matrix indicating linear effects between latent variables. Defaults to an M x M  matrix containing only zeroes.
}
  \item{omega_theta}{
A P x P matrix encoding the residual network structure. By default, theta is modeled instead.
}
  \item{delta_theta}{
A P x P diagonal scaling matrix. Defaults to NA on all diagonal elements. Only used if omega_theta is modeled.
}
  \item{omega_psi}{
An M x M  matrix containing the latent network structure. Dy default, psi is modeled instead.
}
  \item{delta_psi}{
A diagonal M x M scaling matrix. Defaults to an identity matrix. Only used if omega_psi is modeled.
}
  \item{psi}{
An M x M  variance-covariance matrix between latents and latent residuals. Defaults to a full NA matrix.
}
  \item{theta}{
A P x P variance-covariance matrix of residuals of the observed variables. Defaults to a diagonal matrix containing NAs
}
  \item{sampleSize}{
The sample size, only used if \code{data} is assigned a variance-covariance matrix.
}
  \item{fitInd}{
The fit of the independence model. Used to speed up estimation fitting multiple models.
}
  \item{fitSat}{
The fit of the saturated model. Used to speed up estimation fitting multiple models.
}
  \item{startValues}{
An optional named list containing starting values of each model. e.g., \code{list(lambda = matrix(1,9,3))} would set the starting values of a 10 x 3 lambda matrix to ones.
}
\item{scale}{
Logical, should data be standardized before running lvnet?
}
\item{nLatents}{
The number of latents. Allows for quick specification when \code{lambda} is missing. Not needed is \code{lambda} is assigned.
}
\item{lasso}{
The LASSO tuning parameter.
}
\item{lassoMatrix}{
Character vector indicating the names of matrices to apply LASSO regularization on. e.g., \code{"omega_psi"} or  \code{"omega_theta"}.
}
\item{lassoTol}{
Tolerance for absolute values to be treated as zero in counting parameters.
}
\item{ebicTuning}{Tuning parameter used in extended Bayesian Information Criterion.}
\item{mimic}{If set to \code{"lavaan"} (default), covariance matrix is rescaled and N is used rather than N - 1 in likelihood computation.}
}
\details{
The modeling framework follows the all-y LISREL framework for Structural Equation Models (SEM; Hayduk, 1987) to model relationships between P observed variables and M latent variables:

sigma = lambda * (I - beta)^(-1) psi (I - beta)^(-1 T) * lambda^T + theta

Where Sigma is the P x P model-implied covariance matrix, lambda a P x M matrix of factor loadings, B an M x M matrix containing regression effects between latent variables, Psi a M x M covariance matrix of the latent variables/residuals and Theta a P x P covariance matrix of residuals of the observed indicators.

The lvnet function allows for two extensions of this modeling framework. First, psi can be chosen to be modeled as follows:

psi = delta_psi (I - omega_psi)^(-1) delta_psi

In which delta_psi is a M x M diagonal scaling matrix and omega_psi a M x M matrix containing zeroes on the diagonal and partial correlation coefficients on the offdiagonal values of two latent variables conditioned on all other latent variables. omega_psi therefore corresponds to a Gaussian Graphical Model, or a network structure.

Similarly, theta can be chosen to be modeled as follows:

theta = delta_theta (I - omega_theta)^(-1) delta_theta

In which delta_theta is a P x P diagonal scaling matrix and omega_theta a P x P matrix containing zeroes on the diagonal and partial correlation coefficients on the offdiagonal values of two residuals conditioned on all other residuals.

Modeling  omega_psi is termed Latent Network Modeling (LNM) and modeling omega_theta is termed Residual Network Modeling (RNM). lvnet automatically chooses the appropriate modeling framework based on the input.
}
\value{
An \code{lvnet} object, which is a list containing the following elements:
\item{matrices}{A list containing thee estimated model matrices}
\item{sampleStats}{A list containing the covariance matrix (\code{covMat}) and sample size {\code{sampleSize}}}
\item{mxResults}{The OpenMx object of the fitted model}
\item{fitMeasures}{A named list containing the fit measures of the fitted model}
}
\references{
Boker, S. M., Neale, M., Maes, H., Wilde, M., Spiegel, M., Brick, T., ... Fox, J. (2011).  OpenMx:  an open source extended structural equation modelingframework. Psychometrika, 76(2), 306-317

Boker, S. M., Neale, M. C., Maes, H. H., Wilde, M. J., Spiegel, M., Brick, T. R., ..., Team OpenMx. (2014). Openmx 2.0 user guide [Computer software manual].

Hayduk, L. A.  (1987).Structural equation modeling with LISREL: Essentials advances.  Baltimore, MD, USA: Johns Hopkins University Press.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
\seealso{
\code{\link{lvnetSearch}}
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

# Fit CFA model:
CFA <- lvnet(Data, lambda = Lambda)

# Latent network:
Omega_psi <- matrix(c(
  0,NA,NA,
  NA,0,0,
  NA,0,0
),3,3,byrow=TRUE)

# Fit model:
LNM <- lvnet(Data, lambda = Lambda, omega_psi=Omega_psi)

# Compare fit:
lvnetCompare(cfa=CFA,lnm=LNM)

# Summary:
summary(LNM)

# Plot latents:
plot(LNM, "factorStructure")
}


