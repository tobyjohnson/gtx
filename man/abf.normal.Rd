\name{abf.normal}
\alias{abf.normal}
\title{Calculate approximate Bayes factor (ABF) for normal prior.}
\description{
  Calculates an approximation to the Bayes Factor for an alternative
  model where the parameter beta is a priori normal, against a smaller
  model where beta is zero, by approximating
  the likelihood function with a normal distribution.
}
\usage{
abf.normal(beta, se, priorscale, gridrange = 3, griddensity = 20)
}
\arguments{
  \item{beta}{Vector of effect size estimates.}
  \item{se}{Vector of associated standard errors.}
  \item{priorscale}{Scalar specifying the scale (standard deviation) of the prior
    on true effect sizes.}
  \item{gridrange}{Parameter controlling range of grid for numerical integration.}
  \item{griddensity}{Parameter controlling density of points in grid for numerical integration.}
}
\value{
  A vector of approximate Bayes factors.
}
\details{
  This uses the same normal approximation for the likelihood function as
    \dQuote{Bayes factors for genome-wide association studies:
    comparison with P-values} by John Wakefield, 2009, Genetic Epidemiology
    33(1):79-86 at \url{http://dx.doi.org/10.1002/gepi.20359}.  In that
    work, an analytical expression for the approximate Bayes factor was
    derived, which is implemented in \code{\link{abf.Wakefield}}.  This
    function uses a numerical algorithm to calculate
    the (approximate) Bayes factor, which may be a useful starting point
    if one wishes to
    change the assumptions so that the analytical expression of
    Wakefield (2009) no longer applies (as in \code{\link{abf.t}}).
}
\examples{
data(agtstats)
agtstats$pval <- with(agtstats, pchisq((beta/se.GC)^2, df = 1, lower.tail = FALSE))
max1 <- function(bf) return(bf/max(bf, na.rm = TRUE))
agtstats$BF.normal <- with(agtstats, max1(abf.Wakefield(beta, se.GC, 0.05)))
agtstats$BF.numeric <- with(agtstats, max1(abf.normal(beta, se.GC, 0.05)))
with(agtstats, plot(BF.normal, BF.numeric)) # excellent agreement
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
