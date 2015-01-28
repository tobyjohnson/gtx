\name{multipheno.OBrien}
\alias{multipheno.OBrien}
\title{Multi-phenotype 1 df test for association}
\description{
  For a set of phenotypes, given summary results for association tests
  for each single phenotype, over a large fixed set of marker genotypes
  (e.g. GWAS results), this function calculates a 1 degree of freedom multi-phenotype
  association test, which has a close correspondence to O'Brien's test.
}
\usage{
multipheno.OBrien(z, cor.use = NULL, cor.method = "spearman")
}
\arguments{
  \item{z}{a matrix of association Z statistics, with rows corresponding
    to markers and columns corresponding to phenotypes.}
  \item{cor.use}{a logical vector of length nrow(z), indicating a
    subset of markers to use to calculate the correlation matrix.}
  \item{cor.method}{a method to calculate the correlation matrix.}
}
\details{
  The function is named for the close correspondence with O'Brien's
  test.

  It is the user's responsibility to ensure that the columns of
  \dQuote{\code{z}} are marginally well calibrated, i.e. over a set of
  null markers, each column of \dQuote{\code{z}} should be standard
  normal marginal to the other columns of \dQuote{\code{z}}.
}
\value{
  A list with three elements.  The first element,
  \dQuote{\code{rmatrix}}, is the estimated test statistic correlation
  matrix.  The second and third elements are both of length
  \code{nrow(z)}; 
  \dQuote{\code{OB}} contains the test statistics, which are chi-squared
  with 1 d.f. under the null, and \dQuote{\code{pval}}
  contains the P-values.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
