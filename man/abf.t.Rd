\name{abf.t}
\alias{abf.t}
\title{Calculate approximate Bayes factor (ABF) for t distribution prior.}
\description{
  Calculates an approximation to the Bayes Factor for an alternative
  model where the parameter beta is a priori t distributed, against a smaller
  model where beta is zero, by approximating
  the likelihood function with a normal distribution.
}
\usage{
abf.t(beta, se, priorscale, df = 1, gridrange = 3, griddensity = 20)
}
\arguments{
  \item{beta}{Vector of effect size estimates.}
  \item{se}{Vector of associated standard errors.}
  \item{priorscale}{Scalar specifying the scale (standard deviation) of the prior
    on true effect sizes.}
  \item{df}{Degrees of freedom for t distribution prior.}
  \item{gridrange}{Parameter controlling range of grid for numerical integration.}
  \item{griddensity}{Parameter controlling density of points in grid for numerical integration.}
}
\value{
  A vector of approximate Bayes factors.
}
\details{
  This uses the same normal approximation for the likelihood function as
    \dQuote{Bayes factors for genome-wide association studies:
    comparison with P-values} by John Wakeley, 2009, Genetic Epidemiology
    33(1):79-86 at \url{http://dx.doi.org/10.1002/gepi.20359}.  However,
    in contrast to that work, a t distribution is used for the prior,
    which means it is necessary to use a numerical algorithm to calculate
    the (approximate) Bayes factor.
    }
\examples{
data(agtstats)
agtstats$pval <- with(agtstats, pchisq((beta/se.GC)^2, df = 1, lower.tail = FALSE))
max1 <- function(bf) return(bf/max(bf, na.rm = TRUE))
agtstats$BF.normal <- with(agtstats, max1(abf.Wakefield(beta, se.GC, 0.05)))
agtstats$BF.t <- with(agtstats, max1(abf.t(beta, se.GC, 0.0208)))
with(agtstats, plot(-log10(pval), log(BF.normal)))
with(agtstats, plot(-log10(pval), log(BF.t)))
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
