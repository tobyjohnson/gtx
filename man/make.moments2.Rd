\name{make.moments2}
\alias{make.moments2}
\title{Build matrix of second moments from subject-specific data.}
\description{
  Performs input checking, flipping of coded and noncoded alleles,
  removal of individuals with missing data, constructs a
  subject-specific data matrix X, and returns X'X.    
}
\usage{
make.moments2(params, phenolist, snpdata, weightvar)
}
\arguments{
  \item{params}{a data frame, see \link{gtx.params}.}
  \item{phenolist}{a list of phenotypic response variables and
    covariates to be included in addition to genetic variables from \code{params}.}
  \item{snpdata}{a list with snpinfo and data, see \link{snpdata}.}
  \item{weightvar}{name of a variable in \code{snpdata$data} containing
    individual weights when approximating the likelihood for GLMs.}
}
\details{
  After extracting the necessary columns, any rows (individuals) with
  \emph{any} missing data are removed.  
}
\value{
  A numeric matrix of second moments.
}
\examples{
data(mthfrex)
xtx <- make.moments2(mthfr.params, c("SBP", "DBP", "SexC", "Age"), mthfrex)
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
