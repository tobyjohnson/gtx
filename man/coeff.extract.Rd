\name{coeff.extract}
\alias{coeff.extract}
\title{Coefficient extractor for fitted models.}
\description{
  Tool to extract coefficients and standard errors from different fitted
  model objects, and return them in a standard format.
}
\usage{
coeff.extract(object)
}
\arguments{
  \item{object}{a fitted model object of class lm, glm, or coxph.}
}
\details{
  coeff.extract detects the class of object and handles it appropriately
  to extract a two column data frame with columns \dQuote{Estimate} and
  \dQuote{Std Error}.

  This function is useful because there does not seem to otherwise be a
  uniform interface to all these classes.
}
\author{
  Toby Johnson \email{T.Johnson@qmul.ac.uk}
}