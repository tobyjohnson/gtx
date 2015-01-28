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
  \item{xtwx}{an object of class \link{moments2}, typically built using
    \code{\link{make.moments2}} with the \code{weightvar} argument set,
    or a matrix of weighted second moments.}
  \item{leftvar}{name of the response variable (the left hand side of
    the formula).}
  \item{rightvars}{name(s) of the explanatory variables (the right hand
    side of the formula).}
  \item{n}{sample size, only needed for the normal linear model if there is not a single intercept
    \dQuote{\code{ONE}} for all individuals.}
  \item{vscale}{parameter needed if \code{xtwx} is not of class \link{moments2}, set to \code{NULL} for normal linear model and 1 for
    logistic regression.}
}
\details{
  Variables in \code{rightvars} with non-identifiable coefficients are
  removed, with preference for keeping variables that occur earlier
  rather than later in \code{rightvars}.

  When the \code{vscale} attribute of \code{xtwx} (or the \code{vscale}
  function argument) is \code{NULL}, this function assumes
  that the \code{xtwx} argument was calculated with unit weights and
  therefore that a linear model fit is required with error variance
  estimated from the data.  For this application it is preferred to call
  \code{\link{lm.moments2}}, which is a wrapper for this function with
  \code{vscale=NULL}.

  When the \code{vscale} attribute of \code{xtwx} (or the \code{vscale}
  function argument) is set equal to 1, this function assumes that the
  \code{xtwx} argument was calculated with weights calculated such that
  a GLS problem has been correctly set up to approximate a likelihood
  function, and therefore that generalised linear model fit is required.

  Values other than \code{NULL} or 1 for the \code{vscale} parameter may not be
  what you think.  Do not use other values unless you are absolutely sure what you understand
  what are doing.  See the source code for details.
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
myglm <- est.moments2(xtwx, "HTNstar", c("ONE", "rs6668659_T", "rs4846049_T",
                                "rs1801133_G", "SexC", "Age"), vscale=1)
myglm$z <- myglm$betahat/myglm$se
cbind(beta = myglm$betahat, se = myglm$se, z = myglm$z, 
      pval = pnorm(-abs(myglm$z))*2)

## Compare against results from glm
## Note have to use coded alleles used in original data
mycheck <- glm(HTN ~ rs6668659_G+rs4846049_G+rs1801133_A+Sex+Age,
               family="binomial", data = mthfrex$data)
coef(summary(mycheck))
## Note in results Sex factor coded differently than SexC
## Coefficients for covariates used in null model are different,
##     because xtwx approximates around the fitted null model

## Look at pairwise correlations
cor(subset(mthfrex$data, select = c("rs6668659_G", "rs4846049_G",
                                    "rs1801133_A")))^2

## SNP coefficients well approximated (given very high
## inter-SNP correlations) but signs ALL inverted by coded allele flips

## check error less than 10percent
stopifnot(all(-1*myglm$z[2:4]/coef(summary(mycheck))[2:4,3] > 0.9))
stopifnot(all(-1*myglm$z[2:4]/coef(summary(mycheck))[2:4,3] < 1.1))
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
