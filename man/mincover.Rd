\name{mincover}
\alias{mincover}
\title{Compute minimum size of cover of overlapping intervals.}
\description{
  Compute minimum size of cover of overlapping intervals.
}
\usage{
mincover(x.begin, x.end)
}
\arguments{
  \item{x.begin}{a vector of interval begin positions.}
  \item{x.end}{a vector of interval end positions.}
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