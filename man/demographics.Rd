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
\method{demographics}{numeric}(object, by, style)
\method{demographics}{logical}(object, by, style)
\method{demographics}{factor}(object, by, style)
\method{demographics}{Surv}(object, by, style)
\method{demographics}{default}(object, by, style)
\method{demographics}{data.frame}(object, by, style)
}
\arguments{
  \item{object}{an object to summarise, either a vector or data frame.}
  \item{by}{a list of logical vectors, or a character vector, specifying
    the groups to summarize by.}
  \item{style}{details}
}
\value{
  A character matrix.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}