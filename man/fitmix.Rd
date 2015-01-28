\name{fitmix}
\alias{fitmix}
\title{Fit finite mixture of univariate Gaussian densities to data.}
\description{
  Implementation of EM algorithm to fit k component univariate Gaussian
  mixture to data, for user specified value of k.  In constrast to
  general purpose (unconstrained) Gaussian mixture models, this function
  allows certain restrictions or parameter space reductions that make
  the model correspond more closely to a classical quantitative genetics model.
}
\usage{
fitmix(x, k, tol = 1e-6, maxit = 100, restarts = 20,
       p.binomial = FALSE, mu.additive = FALSE, sigma.common = FALSE)
}
\arguments{
  \item{x}{Real vector of data to be fitted by the model.}
  \item{k}{Number of components in the mixture.}
  \item{tol}{Threshold for log-likelihood increase
    for convergence.}
  \item{maxit}{Maximum number of iterations for each run of the EM algorithm.}
  \item{restarts}{Number of times to restart EM algorithm at random
    initial points.}
  \item{p.binomial}{Assume mixture proportions correspond to binomial
    expansion (corresponds to Hardy-Weinberg for k=3).}
  \item{mu.additive}{Assume means follow additive model.}
  \item{sigma.common}{Assume variances are same for all components.}
}
\details{
  Most likely applications are \code{k=2} with \code{sigma.common=TRUE}
  for a completely dominant or recessive genetic model, and \code{k=3}
  with \code{p.binomial=TRUE}, \code{mu.additive=TRUE} and
  \code{sigma.common=TRUE} for a codominant biallelic genetic model.

  Note that the unconstrained model has many global optima that allow
  unbounded increase in the log likelihood as one (or more) mixture
  components have sigma tend to zero as mu tends to one of the data
  observed values.  The EM algorithm rarely converges to these optima
  since they typically have very small domains of attraction.
  }
\value{
  A list.  The elements p, mu, and sigma, are each vectors of length
  \code{k}, and describe the mixture found that maximised the
  likelihood.  An element loglik is a vector of length \code{restarts+1}
  and gives the log likelihood at the end of each run of the EM
  algorithm.  Inspecting \code{loglik} may give some indication of
  whether a \dQuote{good} local optimum has been found.
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
