\name{mincover}
\alias{mincover}
\title{Compute minimum size of cover of overlapping intervals.}
\description{
  Compute minimum size of cover of overlapping intervals.
}
\usage{
mincover(x.begin, x.end, trunc.begin = NULL, trunc.end = NULL)
}
\arguments{
  \item{x.begin}{a vector of interval begin positions.}
  \item{x.end}{a vector of interval end positions.}
  \item{trunc.begin}{truncaction value.}
  \item{trunc.end}{truncation value.}
}
\details{
Note all arguments are coerced to integers.
}
\value{
  An integer, the total size covered by the union of all intervals.
}
\examples{
mincover(c(1, 5, 10, 11, 22), c(8, 17, 13, 19, 25))
## first to fourth intervals all overlap
## third interval 10:13 entirely inside second interval 5:17
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}