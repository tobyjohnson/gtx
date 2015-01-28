\name{coeff.extract}
\alias{coeff.extract}
\alias{coeff.extract.default}
\title{Coefficient extractor for fitted models.}
\description{
  Tool to extract coefficients and standard errors from different fitted
  model objects, and return them in a standard format.
}
\usage{
\method{coeff.extract}{default}(object)
}
\arguments{
  \item{object}{a fitted model object.}
}
\details{
  coeff.extract is a generic function to extract a two column data frame
  with columns \dQuote{Estimate} and \dQuote{Std Error}, and one row for
  each coefficient in the fitted model \code{object}.

  The default method simply packages the values returned by using the generic
  functions \code{coef()} and \code{vcov()}.  It therefore works for all classes that
  these methods work for, including lm, glm, and coxph.

  Written as a generic function so that other fitted model objects might
  be handled.
}
\value{
  A data frame with two variables, always named \dQuote{Estimate}
  and \dQuote{Std. Error}.  Row names indicate model terms.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}