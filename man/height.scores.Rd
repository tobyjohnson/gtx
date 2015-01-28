\name{height.scores}
\alias{height.scores}
\docType{data}
\title{Genetic risk score for height.}
\description{
  Risk score for height parameterised using GWAS meta-analysis results
  published by Lango Allen et al. (Nature 2010) as part of work by the
  GIANT consortium.
}
\usage{
data(height.scores)
}
\format{
  A data frame suitable for use with other functions in this package,
  see \link{gtx.params}.

  The score is called HEIGHT2010, and has coefficients in height Z score
  units (i.e. phenotypic standard deviations) per coded
  allele.  One phenotypic standard deviation for height is 0.06-0.07m in
  most samples analysed by Lango Allen et al.

  Effect size estimates are taken from the \dQuote{STAGE 2} results in
  Supplementary Table 1 of Lango Allen et al.
}
\source{
  For the publication by Lango Allen et al. from which these data were extracted see \url{http://dx.doi.org/10.1038/nature09410}.
}
\keyword{datasets}