\name{positivise}
\alias{positivise}
\title{Ensure numeric values are positive}
\description{
  Convenience function to convert non-positive values to missing.
}
\usage{
positivise(x, strict = FALSE)
}
\arguments{
  \item{x}{a vector of numeric values}
  \item{strict}{Output must be strictly positive}
}
\value{
  The vector with negative values replaced by NA.  If
  \code{strict=TRUE}, zero values are also replaced by NA.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
