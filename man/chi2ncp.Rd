\name{chi2ncp}
\alias{chi2ncp}
\alias{chi2power}
\title{Compute non-centrality parameter of chi squared distribution.}
\description{
  Computes the non-centrality parameter of chi squared distribution for
  specified alpha and beta, such that there is probability (power) beta
  of exceeding the critical value for probability (size) alpha for a
  central chi squared distribution
}
\usage{
chi2ncp(alpha, beta, df = 1)
chi2power(alpha, ncp, df = 1)
}
\arguments{
  \item{alpha}{Test size; required tail area for central chi squared distribution.}
  \item{beta}{Test power; required tail area for non-central chi squared distribution.}
  \item{df}{Degrees of freedom for both chi squared distributions.}
  \item{ncp}{Non-centrality parameter}
}
\details{
  See the examples.
} 
\value{
  \code{chi2ncp} returns the non-centrality parameter (\code{ncp}) for
  given \code{alpha} and \code{beta}.  \code{chi2power} returns the
  power (\code{beta}) for given \code{alpha} and \code{ncp}.
}
\examples{
## NCP needed to have 0.80 power for a 0.05 size test
chi2ncp(.05, .8)
## NCP needed to have 0.80 power for genome-wide significance
chi2ncp(5e-08, .8)
## test
critval <- qchisq(5e-08, lower.tail = FALSE, df = 1)
pchisq(critval, ncp = chi2ncp(5e-08, .8), lower.tail = FALSE, df = 1)
chi2power(5e-08, ncp = chi2ncp(5e-08, .8))
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
