\name{twopq}
\alias{twopq}
\title{Calculate binomial variance 2p(1-p)}
\description{
  Convenience function to calculate variance of a binomial variable with
  n=2 and p=p, equal to the variance of coded genotype under an additive
  genetic model and Hardy-Weinberg proportions.
}
\usage{
twopq(p)
}
\arguments{
  \item{p}{a vector of allele frequencies}
}
\value{
  The vector of binomial variances, 2p(1-p).
}
\examples{
  MAF <- rbeta(100, shape1 = 0.6, shape2 = 1)
  plot(MAF, twopq(MAF), xlim = c(0, 1), ylim = c(0, 0.5))
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
