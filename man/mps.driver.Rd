\name{mps.driver}
\alias{mps.driver}
\alias{mps.driver1}
\alias{mps.summary}
\title{Driver function for Modular Power Simulations}
\description{
  Driver function that calls modules (functions) for separate parts of a power
  simulation.
}
\usage{
mps.driver(design, xFun, yFun, tFun, nrep = 1000,
           alpha = 0.05, sided = "two-sided", verbose = FALSE)
mps.driver1(design1, xFun, yFun, tFun, nrep = 1000)
mps.summary(zz, alpha = 0.05, sided = "two-sided")
}
\arguments{
  \item{design}{A data frame of design points}
  \item{xFun}{A function to generate independent data, evaluated once
    for each design point} 
  \item{yFun}{A function to generate dependent data, evaluated
    \code{nrep} times for each design point} 
  \item{tFun}{A list of functions to test each realisation of the
    dependent data}
  \item{nrep}{The number of replicates for simulating dependent data
    for each design point}
  \item{alpha}{A vector of test significance levels}
  \item{sided}{A vector of test sidedness}
  \item{verbose}{Whether to print progress messages}
  \item{design1}{One row of a design matrix}
  \item{zz}{Output of mps.driver1}
}
\details{
  The idea of Modular Power Simulations is to provide a generic driver
  function \code{mps.driver()}, to which \dfn{modules} (other functions
  such as \code{\link{ssm.QT}()}, \code{\link{ssm.LM}()}, etc.) are
  passed to specify the simulation model and the test statistic(s) to be
  evaluated.  The driver function handles rote looping over design
  points, and processing of the test statistics from replicate
  simulations to calculate power at different alpha levels and test
  sidednesses, bias, and other summaries.  A modular framework is
  intended to facilitate code re-use, because individual module
  functions can be re-used or modified for different applications.

  The practical use of \code{mps.driver()} is best illustrated by the
  examples and vignettes.  This documentation concentrates on the
  technical specification.
  
  The modularity is implemented by defining two separate functions
  that specify the data generation mechanism, and a further list of
  separate functions that specify test statistic(s).  The data generation
  mechanism is separated into \code{xFun} and \code{yFun}.  \code{xFun}
  will be evaluated only once for each design point, is typically
  deterministic, and may generate independent variables in regression-like
  applications.  \code{yFun} will be evaluated many replicate times for
  each design point, is typically stochastic, and may generate dependent
  variables in regression-like applications.  The return value from
  \code{xFun} is passed as an argument to \code{yFun}; this can be any
  data structure.
  
  Parameter values for the simulation model may be specified for each
  design point, or may be specified globally.  This is achieved by using
  free (unbound) variables inside the functions for individual modules.
  For example, in a function to simulate randomised treatment status:
  \preformatted{simpleTreatment <- function() return(rbinom(sampleSize, 1, 0.5))}
  the variable \code{sampleSize} is a free variable.  To use this
  function with \code{mps.driver()}, there must either be a column in
  the design matrix called \dQuote{sampleSize}, or \code{sampleSize}
  must be defined in the global environment.

  To simulate for the parameter values at a particular the design point,
  the driver function manipulates the environment of the simulation
  functions, essentially by doing:
  \preformatted{design1 <- design[designIdx, , drop = FALSE]
with(as.list(design1), {
  environment(xFun) <- environment()
  environment(yFun) <- environment()
  x1 <- xFun()
  y1 <- yFun(x1)
  ...
})}
  to simulate one one realisation (\code{x1}, \code{y1}) for the
  \code{designIdx}-th design point.

  The return values of \code{xFun} and \code{yFun} can be any data
  structures; they are passed as arguments to a function or functions
  that calculate the test statistics for which power calculations are
  required.  The test statistic(s) are specified using \code{tFun}, which
  should be a named list of functions.  Note if \code{tFun} is a single
  function, it will be made into a list using
  \preformatted{tFun <- list(test = tFun)}
  There is a very specific requirement that each element of the
  list \code{tFun} be a function whose return value can
  be coerced to a double of constant length, with at least a named
  element \dQuote{pval}.  It is \strong{strongly} recommended to use or
  extend the generic function \code{\link{test.extract}} to achieve this,
  e.g. as in
  \preformatted{tFun = list(test1 = function(x, y)
    return(test.extract(lm(phenotype ~ genotypeA, data = y))))}
  See \code{\link{ssm.LM}}, \code{\link{ssm.GLM}},
  \code{\link{ssm.CoxPH}} for further examples.
	  
  \code{mps.driver()} treats each row of the \code{design} argument as a
  separate design point, and generates \code{nrep} replicate simulations.
  Each test function is applied to the same set of replicate
  simulations, and power is evaluated at significance levels \code{alpha}
  and test sidednesses \code{sided}.
}
\value{
  A data frame of results.  There is one row for each combination of
  design point, test function, and significance/sidedness.  That is, the
  return value will have \code{nrow(design)*length(tFun)*length(alpha)}
  rows. 
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
