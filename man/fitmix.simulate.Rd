\name{fitmix.simulate}
\alias{fitmix.simulate}
\title{Simulate from finite mixture of univariate Gaussian densities.}
\description{
  Simulates random variables from a user specified finite mixture of
  univariate Gaussian densities.
}
\usage{
fitmix.simulate(n, p, mu, sigma)
}
\arguments{
  \item{n}{Number of draws from target distribution.}
  \item{p}{Real vector of mixture proportions.}
  \item{mu}{Real vector of mixture component means.}
  \item{sigma}{Real vector of mixture component standard deviations.}
}
\value{
  A vector of simulated values.
}
\examples{
xx <- fitmix.simulate(100, c(0.49, 0.42, 0.09), c(0, 1, 2), c(.3, .3, .3))
plot(density(xx))
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
