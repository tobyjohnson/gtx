\name{grs.filter.Qrs}
\alias{grs.filter.Qrs}
\title{Filter SNPs for inclusion in genetic risk score using
  heterogeneity test.}
\description{
  Performs a stepwise downward "model selection" in which SNPs are
  iteratively removed from the risk score until the heterogeneity test
  is no longer significant at the specified threshold.
}
\usage{
grs.filter.Qrs(w, b, s, p.thresh = 0.05)
}
\arguments{
  \item{w}{coefficients for the risk score}
  \item{b}{aligned beta coefficients in the testing dataset}
  \item{s}{standard errors}
  \item{p.thresh}{P-value threshold}
}
\details{
  \strong{This function needs better documentation!}
}
\author{
  Toby Johnson \email{T.Johnson@qmul.ac.uk}
}
