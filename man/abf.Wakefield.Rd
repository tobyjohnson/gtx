\name{abf.Wakefield}
\alias{abf.Wakefield}
\title{Calculate approximate Bayes factor (ABF) using method of
  Wakefield (2009).}
\description{
  Calculates an approximation to the Bayes Factor for an alternative
  model where the parameter beta is a priori normal, against a smaller
  model where beta is zero, by approximating
  the likelihood function with a normal distribution.
}
\usage{
abf.Wakefield(beta, se, priorsd)
}
\arguments{
  \item{beta}{Vector of effect size estimates.}
  \item{se}{Vector of associated standard errors.}
  \item{priorsd}{Scalar specifying the standard deviation of the prior
    on true effect sizes.}
}
\value{
  A vector of approximate Bayes factors.
}
\details{
  See \dQuote{Bayes factors for genome-wide association studies:
    comparison with P-values} by John Wakefield, 2009, Genetic
  Epidemiology 33(1):79-86 at \url{http://dx.doi.org/10.1002/gepi.20359}.
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
