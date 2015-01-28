\name{fitmix.plot}
\alias{fitmix.plot}
\title{Plot empirical density and components and total density for finite mixture of univariate Gaussian densities.}
\description{
  Draws a nice plot.
}
\usage{
fitmix.plot(x, p, mu, sigma)
}
\arguments{
  \item{x}{Real vector of data.}
  \item{p}{Real vector of mixture proportions.}
  \item{mu}{Real vector of mixture component means.}
  \item{sigma}{Real vector of mixture component standard deviations.}
}
\value{
  Returns an invisible null.
}
\examples{
xx <- fitmix.simulate(100, c(0.49, 0.42, 0.09), c(0, 1, 2), c(.3, .3, .3))

## additive model, common variance, Hardy--Weinberg
fit.a <- fitmix(xx, 3, maxit = 10, restarts = 3,
                sigma.common = TRUE, p.binomial = TRUE, mu.additive = TRUE)
## general (unrestricted) fit
fit.g <- fitmix(xx, 3, maxit = 10, restarts = 3)

oldpar <- par(mfrow = c(1, 2))
fitmix.plot(xx, fit.a$p, fit.a$mu, fit.a$sigma)
fitmix.plot(xx, fit.g$p, fit.g$mu, fit.g$sigma)
par(oldpar)
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
