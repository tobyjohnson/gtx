\name{moments2}
\alias{is.moments2}
\alias{moments2}
\title{Class for summary statistic matrix of second moments.}
\description{
  Objects of this class are a matrix of second moments between a set of
  variables to be used in regression analyses.
}
\usage{
is.moments2(object)
}
\arguments{
  \item{object}{an object to test.}
}
\details{
  An object of class \code{moments2} is square matrix, with an attribute
  \dQuote{\code{vscale}} set to NULL or 1 depending on whether the
  moments are unweighted (for normal linear model analyses) or weighted
  (for GLS approximation of generalised linear model analysis).

  Typically created by calling \code{\link{make.moments2}()} and used in
  calls to \code{\link{lm.moments2}()} and \code{\link{est.moments2}()}.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
