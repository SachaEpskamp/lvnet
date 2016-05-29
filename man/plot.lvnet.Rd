\name{plot.lvnet}
\alias{plot.lvnet}
\alias{plot.lvnetSearch}
\alias{plot.lvnetLasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot model matrices
}
\description{
Plot method for \code{lvnet}. For \code{lvnetSearch} and \code{lvnetLasso} objects this is simply defined as \code{plot(object$best, ...)}
}
\usage{
 \method{plot}{lvnet}(x, what = c("factorStructure", "residual", "latent"), partial, 
        layout = "circle", ...)
 \method{plot}{lvnetLasso}(x, ...)
 \method{plot}{lvnetSearch}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An \code{lvnet} object.
}
  \item{what}{
What to plot? \code{"factorStructure"} plots the factor loadings and latent correlations or network. \code{"residual"} the residual correlations or network and \code{"latent"} the latent correlations or network.
}
  \item{partial}{
Plot partial correlations instead of correlations? Defaults to \code{TRUE} if \code{omega_psi} or \code{omega_theta} is estimated.
}
  \item{layout}{
The layour argument as used in \code{\link[qgraph]{qgraph}}
}
  \item{\dots}{
Arguments sent to \code{\link[qgraph]{qgraph}}
}
}

\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
