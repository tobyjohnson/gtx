\name{test.extract}
\alias{test.extract}
\alias{test.extract.htest}
\alias{test.extract.survdiff}
\alias{test.extract.lm}
\alias{test.extract.glm}
\alias{test.extract.coxph}
\title{Test statistic extractor.}
\description{
  A function to extract a test statistic or summary from different test and fitted
  model objects, and return in a standard format.
}
\usage{
\method{test.extract}{htest}(object, partial = NULL)
\method{test.extract}{survdiff}(object, partial = NULL)
\method{test.extract}{lm}(object, partial = NULL)
\method{test.extract}{glm}(object, partial = NULL)
\method{test.extract}{coxph}(object, partial = NULL)
}
\arguments{
  \item{object}{a fitted model object.}
  \item{partial}{character name of coefficient for partial test.}
}
\details{
  The \code{partial} argument only works for numeric terms, not
  factors.  Two level factors should be recoded as 0/1.  Factors with
  more than two levels could in principle be detected and handled
  appropriately, but this is not currently implemented.
    
  test.extract is implemented as a generic function, but handling of
  different tests of class htest requires if-else-if

  The returned effect might be an effect size, but its main purpose is
  to indicate effect direction.

  The df returned is always the df in the LRT sense, i.e. the number of
  additional free parameters.  Hence not the same as the df reported by
  some htest, e.g. not the same as the df for a t-test.
}
\value{
  A named real vector with three elements, named \dQuote{pval},
  \dQuote{effect} and \dQuote{df}.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}