\name{est.moments2}
\alias{est.moments2}
\title{Estimate regression coefficients using quadratic approximation to likelihood function.}
\description{
  To make a quadratic approximation to the likelihood function, the score
  and information are obtained from a pre-built matrix of weighted
  second moments.  This allows a parameter estimate to be obtained by
  one iteration of weighted least squares, or equivalently a score test.
  Typically the weights used to construct the pre-built matrix
  correspond to the MLE under a chosen null model.
}
\usage{
est.moments2(xtwx, leftvar, rightvars, n = NULL, vscale = NULL)
}
\arguments{
  \item{xtwx}{a matrix of weighted second moments, typically built using
    \code{\link{make.moments2}} with the \code{weightvar} argument set.}
  \item{leftvar}{name of the response variable (the left hand side of
    the formula).}
  \item{rightvars}{name(s) of the explanatory variables (the right hand
    side of the formula).}
  \item{n}{sample size, only needed for the normal linear model if there is not a single intercept
    \dQuote{\code{ONE}} for all individuals.}
  \item{vscale}{parameter, set to \code{NULL} for normal linear model and 1 for
    logistic regression.}
}
\details{
  Non-identifiable variables in \code{rightvars} are removed, with
  preference for keeping variables that occur earlier rather than later
  in \code{rightvars}.

  When the \code{vscale} argument is \code{NULL} this function assumes
  that the \code{xtwx} argument was calculated with unit weights and
  therefore that a linear model fit is required with error variance
  estimated from the data.  For this application it is preferred to call
  \code{\link{lm.moments2}}, which is a wrapper for this function with
  \code{vscale=NULL}.

  When the \code{vscale} argument is set equal to 1 this function assumes
  that the \code{xtwx} argument was calculated with weights calculated
  such that a correct likelihood function can be recovered and therefore
  that a generalised linear model fit is required.

  Values other than \code{NULL} or 1 for the \code{vscale} parameter may not be
  what you think.  Do not use other values unless you are absolutely sure what you understand
  what are doing.  See the manuscript for details.
}
\value{
  A list with slots for the effect size estimates, standard errors, and
  a precision matrix.
}
\examples{
data(mthfrex)
mthfrex <- gls.approx.logistic(mthfrex, "HTN", c("SexC", "Age"))
xtwx <- make.moments2(mthfr.params, c("HTNstar", "SexC", "Age"), mthfrex,
                      weightvar = "weight")
est.moments2(xtwx, "HTNstar", c("ONE", "rs6668659_T", "rs4846049_T",
                                "rs1801133_G", "SexC", "Age"), vscale=1)
## Compare against results from glm
## Note have to use coded alleles in original data
glm(HTN ~ rs6668659_G+rs4846049_G+rs1801133_A+Sex+Age,
    family="binomial", data = mthfrex$data)
## Note in results Sex factor coded differently than SexC
## Look at pairwise correlations
cor(subset(mthfrex$data, select = c("rs6668659_G", "rs4846049_G",
                                    "rs1801133_A")))^2
}
\author{
  Toby Johnson \email{T.Johnson@qmul.ac.uk}
}
