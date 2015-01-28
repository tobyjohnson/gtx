\name{gclambda}
\alias{gclambda}
\title{Compute genomic control coefficient}
\description{
  Computes a measure from departure from uniform [0,1] distribution
  called the genomic control coefficient, lambda.
}
\usage{
gclambda(p, df = 1, prob = 0.5)
}
\arguments{
  \item{p}{A vector of observed p-values}
  \item{df}{Degrees of freedom of chi-squared distribution}
  \item{prob}{The quantile of the chi-squared distribution}
}
\value{
  The genomic control coefficient. 
}
\examples{
gclambda(runif(1000))
gclambda(rbeta(1000, shape1 = 0.6, shape2 = 1))
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
