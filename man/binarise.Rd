\name{binarise}
\alias{binarise}
\title{Derive zero-or-one binary variable}
\description{
  Convenience function to derive a zero-or-one binary vector from two
  logical vectors.
}
\usage{
binarise(is0, is1)
}
\arguments{
  \item{is0}{a logical vector}
  \item{is1}{a logical vector of same length as is0}
}
\value{
  Returns a vector of 0/1/NA of the same length as the input arguments,
  with 0 when \code{is0} is TRUE and \code{is1} is FALSE, 1 when
  \code{is1} is TRUE and \code{is0} is FALSE, and NA otherwise.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
