\name{chi2ncp}
\alias{chi2ncp}
\title{Compute non-centrality parameter of chi squared distribution.}
\description{
  Computes the non-centrality parameter of chi squared distribution for
  specified alpha and beta, such that there is probability (power) beta
  of exceeding the critical value for probability (size) alpha for a
  central chi squared distribution
}
\usage{
chi2ncp(alpha, beta, df = 1)
}
\arguments{
  \item{alpha}{Required tail area for central chi squared distribution.}
  \item{beta}{Required tail area for non-central chi squared distribution.}
  \item{df}{Degrees of freedom for both chi squared distributions.}
}
\details{
  See the examples.
} 
\value{
  The non-centrality parameter.
}
\examples{
## 0.80 power for 0.05 size test
chi2ncp(.05, .8)
## 0.80 power for genome-wide significance
chi2ncp(5e-08, .8)
## test
critval <- qchisq(5e-08, lower.tail = FALSE, df = 1)
pchisq(critval, ncp = chi2ncp(5e-08, .8), lower.tail = FALSE, df = 1)
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
