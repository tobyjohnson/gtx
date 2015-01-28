\name{grs.make.scores}
\alias{grs.make.scores}
\title{Make genetic risk scores from individual-level data.}
\description{
  Constructs the value of one or more genetic risk scores with specified
  parameterisation, from individual-level genotype data.  Automatically
  handles coded allele designation.
}
\usage{
grs.make.scores(params, snpdata, appendage = ".score")
}
\arguments{
  \item{params}{a data frame, see \code{\link{gtx.params}}.}
  \item{snpdata}{a list with snpinfo and data, see \code{\link{snpdata}}.}
  \item{appendage}{text to add to score name to make new column name in \code{snpdata$data}.}
}
\details{
  This function computes subject-specific values of one or more additive
  multi-SNP genetic risk scores using a supplied parameterisation and
  subject-specific genotype and phenotype data.  These risk scores are
  returned as additional columns on the input genotype data frame,
  making it straightforward to use them as explanatory variables for any
  of the model fitting functions available in the R language, such as
  lm, glm, or coxph.
}
\value{
  The snpdata list with additional column(s) added to data.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
