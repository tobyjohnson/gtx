\name{lm.moments2}
\alias{lm.moments2}
\title{Fit normal linear model using pre-built matrix of second moments.}
\description{
  Regression coefficients for a normal linear model are be calculated,
  using a pre-built sufficient summary statistic matrix that contains
  the second moments between all the variables of (potential) interest.
}
\usage{
lm.moments2(xtx, leftvar, rightvars, n = NULL)
}
\arguments{
  \item{xtx}{a matrix of second moments, typically built using \code{\link{make.moments2}}.}
  \item{leftvar}{name of the response variable (the left hand side of
    the formula).}
  \item{rightvars}{name(s) of the explanatory variables (the right hand
    side of the formula).}
  \item{n}{sample size, only needed if there is not a single intercept
    for all individuals.}
}
\details{
  Non-identifiable variables in \code{rightvars} are removed, with
  preference for keeping variables that occur earlier rather than later
  in \code{rightvars}.
}
\value{
  A list with slots for the effect size estimates, standard errors, and
  a precision matrix.
  
  This function is just a wrapper for calling \code{\link{est.moments2}} with
  \code{scale = NULL}.
}
\examples{
data(mthfrex)
xtx <- make.moments2(mthfr.params, c("SBP", "DBP", "SexC", "Age"), mthfrex)
lm.moments2(xtx, "SBP", c("ONE", "rs6668659_T", "rs4846049_T",
                                "rs1801133_G", "SexC", "Age"))
## Compare against results from lm
## Note have to use coded alleles in original data
lm(SBP ~ rs6668659_G+rs4846049_G+rs1801133_A+Sex+Age, data = mthfrex$data)
## Note in results Sex factor coded differently than SexC
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
