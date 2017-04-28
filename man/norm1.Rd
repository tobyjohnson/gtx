\name{norm1}
\alias{norm1}
\title{Normalise to sum to 1.}
\description{
  Normalises a vector so that the elements sum to 1, in a relatively
  numerically stable way.
}
\usage{
norm1(x, log = FALSE)
}
\arguments{
  \item{x}{Values to be normalised.}
  \item{log}{Whether values of \code{x} are logged and should be
    exponentiated when normalising.}
}
\value{
  A vector of normalised values that sums to 1.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}

