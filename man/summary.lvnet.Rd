\name{summary.lvnet}
\alias{summary.lvnet}
\alias{summary.lvnetSearch}
\alias{summary.lvnetLasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary method for lvnet
}
\description{
Plot method for \code{lvnet}. For \code{lvnetSearch} and \code{lvnetLasso} objects this is simply defined as \code{summary(object$best, ...)}
}
\usage{
 \method{summary}{lvnet}(object, include = c("input", "chisq", "infcrit", "fitindices", 
        "rmsea", "parests"), digits = 3, ...)
  \method{summary}{lvnetLasso}(object, ...)
 \method{summary}{lvnetSearch}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An \code{lvnet} object
}
  \item{include}{
Vector indicating what to include? \code{"input"} for the input used, \code{"chisq"} for the chi-square fit, \code{"infcrit"} for information criteria, \code{"fitindices"} for fit indices,  \code{"rmsea"} for the RMSEA, ans \code{"parests"} for parameter estimates.
}
  \item{digits}{
Number of digits to round to.
}
  \item{\dots}{
Not used.
}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

