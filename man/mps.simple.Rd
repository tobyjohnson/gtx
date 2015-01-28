\name{ssm}
\alias{ssm.null}
\alias{ssm.QT}
\alias{ssm.LM}
\alias{ssm.CC}
\alias{ssm.CCpopulation}
\alias{ssm.CCsample}
\alias{ssm.GLM}
\alias{ssm.Surv}
\alias{ssm.CoxPH}
\title{Simple Simulation Modules (ssm) for Illustrating Modular Power Simulations}
\description{
  Some Simple Simulation Modules, which are functions that can be used
  as-is, or used as templates for more
  complicated functions, to use in Modular Power Simulations.
}
\usage{
ssm.null()
ssm.QT(x)
ssm.LM(x, y)
ssm.CC()
ssm.CCpopulation(x)
ssm.CCsample(x)
ssm.GLM(x, y)
ssm.Surv(x)
ssm.CoxPH(x, y)
}
\arguments{
  \item{x}{The return value of another function}
  \item{y}{The return value of another function}
}
\details{
  The idea of Modular Power Simulations is to provide a generic driver
  function (\code{\link{mps.driver}} documented elsewhere), to which
  functions such as those documented here are passed, to specify the
  simulation model and test statistic(s) to be evaluated.  The driver
  function handles rote looping over design points, and summarising of
  replicate simulation test statistics to calculate power at different
  alpha levels and test sidednesses.  As described in the documentation
  for \code{\link{mps.driver}}, a modular framework is intended to
  facilitate code re-use, because individual modules can be re-used or
  modified for different applications.
  
  The calling environment must resolve sampleSize, alleleFrequency,
  effectSize and dominanceCoeff, etc.
}
\value{
  Various return values. 
}
\examples{
design <- expand.grid(sampleSize = c(1000, 2000),
                      alleleFrequency = c(0.1, 0.3, 0.5),
                      effectSize = c(0.1, 0.2),
                      dominanceCoeff = 0)
mps.driver(design, ssm.null, ssm.QT, list(lm = ssm.LM),
           nrep = 100)
# run with larger nrep for better results
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
