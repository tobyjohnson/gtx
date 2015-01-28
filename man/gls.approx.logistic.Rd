\name{gls.approx.logistic}
\alias{gls.approx.logistic}
\title{Calculate weights and transformed phenotype so that one iteration
  of generalised least squares approximates a logistic regression.}
\description{Logistic regression models are usually fitted by
  iteratively reweighted generalised least squares (GLS).  This function
  formulates a GLS problem by calculating weights and rescaling the
  response variable such that a logistic regression analysis is
  approximated.  This is equivalent to making a quadratic approximation
  to the likelihood and to performing a score test.}
\usage{
gls.approx.logistic(snpdata, leftvar, rightvars = NULL,
                    outvar = paste(leftvar, "star", sep = ""),
                    weightvar = "weight") 
}
\arguments{
  \item{snpdata}{a list with snpinfo and data, see \link{snpdata}.}
  \item{leftvar}{the name (in \code{snpdata$data}) of the response
    variable for the logistic regression.  Must
    be 0/1 or NA for all individuals.}
  \item{rightvars}{the names (in \code{snpdata$data}) of the explanatory 
    variables for the logistic regression.}
  \item{outvar}{the name (to be added to \code{snpdata$data}) of the
    transformed phenotype that is the response variable in the GLS problem.}
  \item{weightvar}{the name (to be added to \code{snpdata$data}) of the
    weight variable in the GLS problem.}
}
\details{
  An intercept term is \emph{always} included and there is currently no
  way to override this.

  The specified model is printed as a check.  MLE parameter values are
  estimated under the specified model by calling \code{glm}.  These are
  used to compute weights and a transformed phenotype such that one
  iteration of generalised least squares constitutes a score test for
  additional terms that might subsequently be added to the model.
  
  This function works on the (usually phenotypic) columns in the
  \code{$data} slot of \code{snpdata}.  
  }
\value{
  Returns \code{snpdata} with additional columns for the weights and transformed
  phenotype.

  The formula used to find an MLE is printed.
}
\examples{
data(mthfrex)
mthfrex <- gls.approx.logistic(mthfrex, "HTN", c("SexC", "Age"))
xtwx <- make.moments2(mthfr.params, c("HTNstar", "SexC", "Age"), mthfrex,
                      weightvar = "weight")
est.moments2(xtwx, "HTNstar", c("ONE", "rs6668659_T", "rs4846049_T",
                                "rs1801133_G", "SexC", "Age"), vscale=1)
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
