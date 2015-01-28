\name{magic.scores}
\alias{magic.scores}
\docType{data}
\title{Genetic risk scores for glucose/insulin traits.}
\description{
  Risk scores parameterised using GWAS meta-analysis results published
  by the MAGIC consortium (Dupuis et al. Nature Genetics 2010).  There
  are scores for fasting plasma glucose (FPG), fasting insulin (FI),
  HOMA-B and HOMA-IR, all parameterised using genome wide significant
  SNPs reported in Tables 1 and 2 of Dupuis et al. (2010).
  
  These data are actually NOT YET INCLUDED.
}
\usage{magic.scores}
\format{
  A data frame suitable for use with other functions in this package,
  see \link{gtx.params}.
  
  The scores are called FPG2010, FI2010, HOMAB2010, HOMAIR2010.  The
  coefficients are in the same units as the source publication (what exactly?)
}
\source{
  For the publication by Dupuis et al. from which these data were extracted see \url{http://dx.doi.org/10.1038/ng.520}.
}
\keyword{datasets}