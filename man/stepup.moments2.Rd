\name{stepup.moments2}
\alias{stepup.moments2}
\title{Stepwise upward model selection using summary statistic matrix.}
\description{
  Terms are added iteratively to a regression model until the
  improvement in fit (judged by the P-value for a partial t-test or
  score test) is no longer signicant at a specified threshold.  The
  method implemented here makes use of a pre-built sufficient summary
  statistic matrix, which contains the (weighted) second moments between
  all the variables that are being assessed for inclusion in the
  regression model.  The calculations are exact for a normal linear
  model and correspond to a score test for a generalised linear model.
}
\usage{
stepup.moments2(xtwx, leftvar, biggest, smallest,
                p.thresh = 0.05, n = NULL, vscale = NULL)
}
\arguments{
  \item{xtwx}{a matrix of (weighted) second moments, typically built using
    \code{\link{make.moments2}}.}
  \item{leftvar}{name of the response variable (the left hand side of
    the formula).}
  \item{biggest}{name(s) of the explanatory variables in the biggest
    model to fit, i.e.\ variables to \emph{consider}
    for inclusion.}
  \item{smallest}{name(s) of the explanatory variables in the smallest
    model to fit, i.e.\ variables that \emph{must} be included.}
  \item{p.thresh}{P-value threshold for proceeding to add a term to the model.}
  \item{n}{sample size, only needed for the normal linear model if there is not a single intercept
    \dQuote{\code{ONE}} for all individuals.}
  \item{vscale}{parameter, set to \code{NULL} for normal linear model and 1 for
    logistic regression.}
}
\details{
  This performs stepwise upward model selection.  Significance of terms
  considered for inclusion is determined using identical calculations to
  \code{\link{lm.moments2}} and \code{\link{est.moments2}}.

  When the \code{vscale} argument is \code{NULL} this function assumes
  that the \code{xtwx} argument was calculated with unit weights and
  therefore that a linear model fit is required with error variance
  estimated from the data.

  When the \code{vscale} argument is set equal to 1 this function assumes
  that the \code{xtwx} argument was calculated with weights calculated
  such that a correct likelihood function can be recovered and therefore
  that a generalised linear model fit is required.

  Values other than \code{NULL} or 1 for the \code{vscale} parameter may not be
  what you think.  Do not use other values unless you are absolutely sure that you understand
  what are doing.  See the manuscript for details.
}
\value{
  The fitted model, as returned by calling
  \code{\link{lm.moments2}} or \code{\link{est.moments2}}.
  This is a list with slots for the effect size estimates, standard errors, and
  a precision matrix.
}
\examples{
data(mthfrex)
xtx <- make.moments2(mthfr.params, c("SBP", "DBP", "SexC", "Age"), mthfrex)
allsnps <- paste(mthfr.params$snp, mthfr.params$coded.allele, sep = "_")
myfit <- stepup.moments2(xtx, "SBP", allsnps, c("ONE", "SexC", "Age"))
## much faster than stepAIC but only steps to bigger models
cbind(beta = myfit$betahat, se = myfit$se, t = myfit$beta/myfit$se)
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
