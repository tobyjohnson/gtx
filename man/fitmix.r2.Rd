\name{fitmix.r2}
\alias{fitmix.r2}
\title{For finite mixture of univariate Gaussian densities, computes
  proportion of variance explained by the mixture labels.}
\description{
  Computes the true R-squared if the mixture labels were known without
  error.  Currently only works if all components have the same variance.
}
\usage{
fitmix.r2(p, mu, sigma)
}
\arguments{
  \item{p}{Real vector of mixture proportions.}
  \item{mu}{Real vector of mixture component means.}
  \item{sigma}{Mixture component standard deviation, same for all components.}
}
\value{
  The R-squared value for a (hypothetical) regression of data values
  onto mixture label, as a \code{k} level factor.
}
\examples{
xx <- fitmix.simulate(100, c(0.49, 0.42, 0.09), c(0, 1, 2), c(.3, .3, .3))

## additive model, common variance, Hardy--Weinberg
fit.a <- fitmix(xx, 3, maxit = 10, restarts = 3,
                sigma.common = TRUE, p.binomial = TRUE, mu.additive = TRUE)
fitmix.plot(xx, fit.a$p, fit.a$mu, fit.a$sigma)
fitmix.r2(fit.a$p, fit.a$mu, unique(fit.a$sigma))
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
