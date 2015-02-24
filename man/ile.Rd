\name{ile}
\alias{ile}
\title{Split data at tertiles, quartiles etc}
\description{
  Transforms numeric data to factor with equally sized groups, split at
  quantiles (tertiles for three groups, quantiles for four groups, etc.)
}
\usage{
ile(x, levels.out)
}
\arguments{
  \item{x}{a vector of values to be transformed}
  \item{levels.out}{number or names of the levels to be output}
}
\value{
  A factor.
}
\examples{
  x <- rnorm(15)
  ile(x, levels.out = c("Low", "Medium", "High"))
  ile(x, 4)
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
