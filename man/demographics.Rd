\name{demographics}
\alias{demographics}
\alias{demographics.numeric}
\alias{demographics.logical}
\alias{demographics.factor}
\alias{demographics.Surv}
\alias{demographics.default}
\alias{demographics.data.frame}
\title{Easy demographics tables.}
\description{
  Functions for easy construction of demographics tables.
}
\usage{
\method{demographics}{numeric}(object, by, style, digits)
\method{demographics}{logical}(object, by, style, digits)
\method{demographics}{factor}(object, by, style, digits)
\method{demographics}{Surv}(object, by, style, digits)
\method{demographics}{default}(object, by, style, digits)
\method{demographics}{data.frame}(object, by, style, digits)
}
\arguments{
  \item{object}{an object to summarise, either a vector or data frame.}
  \item{by}{a list of logical vectors, or a character vector, specifying
    the groups to summarize by.}
  \item{style}{details}
  \item{digits}{numeric precision}
}
\value{
  A character matrix.  See the vignette.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}