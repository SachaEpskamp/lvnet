\name{lav2lvnet}
\alias{lav2lvnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert lavaan model to lvnet model matrices
}
\description{
This function can be used to easily generate input matrices for lvnet based on a lavaan model.
}
\usage{
lav2lvnet(model, data, std.lv = TRUE, lavaanifyOps = list(auto = TRUE, std.lv = std.lv))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
Lavaan model syntax
}
  \item{data}{
The dataset. Only used to extract order of variables names from the columnnames.
}
  \item{std.lv}{
Should the model be identified by constraining latent variable variance to 1. Defaults to \code{TRUE} unlike lavaan! This is because the starting values work better for this identification.
}
  \item{lavaanifyOps}{
A list with other options sent to \code{\link[lavaan]{lavaanify}}
}
}

\value{
A list with the model matrices for \code{lambda}, \code{psi}, \code{theta} and  \code{beta}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

\examples{
\dontrun{
library("lavaan")

# Load dataset:
data(HolzingerSwineford1939)
Data <- HolzingerSwineford1939[,7:15]

# lavaan model
HS.model <- '
visual  =~ x1 + x2 + x3
textual =~ x4 + x5 + x6
speed   =~ x7 + x8 + x9 '

# fit via lavaan:
lavFit <- cfa(HS.model, HolzingerSwineford1939[7:15],std.lv=TRUE)

# Fit via lvnet:
mod <- lav2lvnet(HS.model, HolzingerSwineford1939[7:15])
lvnetFit <- lvnet(Data, lambda = mod$lambda, psi = mod$psi)

# Compare:
Compare <- data.frame(
  lvnet = round(unlist(lvnetFit$fitMeasures)[c("npar","df","chisq","fmin","aic","bic",
                                              "rmsea","cfi","tli","nfi","logl")],3),
  lavaan = round(fitMeasures(lavFit)[c("npar","df","chisq","fmin","aic","bic","rmsea",
                                      "cfi","tli","nfi","logl")],3))

Compare
}
}
