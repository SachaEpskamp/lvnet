\name{lvnetRefit}
\alias{lvnetRefit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Refit lvnet model to new data
}
\description{
Obtain fit indices from the estimated model parameters on a new dataset.
}
\usage{
lvnetRefit(lvnetObject, data, sampleSize)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lvnetObject}{
Output of \code{\link{lvnet}}.
}
  \item{data}{
New dataset or variance-covariance matrix.
}
  \item{sampleSize}{
Sample size (if \code{data} is a variance-covariance matrix).
}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}