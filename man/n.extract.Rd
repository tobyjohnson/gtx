\name{n.extract}
\alias{n.extract}
\alias{n.extract.default}
\alias{n.extract.lm}
\title{Sample size extractor for fitted models.}
\description{
  Tool to extract sample size from different fitted
  model objects, and return it in a standard format.
}
\usage{
\method{n.extract}{default}(object)
}
\arguments{
  \item{object}{a fitted model object.}
}
\details{
  n.extract is a generic function to extract a named scalar \dQuote{n} from the fitted model \code{object}.
  
  Currently the default method assumes \code{object} is a list and tries
  to extract an element called \dQuote{n}.  A method is also provided
  for class \code{lm} (which works also for \code{glm}) which calculates
  the sample size from the number of non-missing residuals.
}
\value{
  The sample size that was analysed.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}