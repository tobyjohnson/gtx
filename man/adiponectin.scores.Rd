\name{adiponectin.scores}
\alias{adiponectin.scores}
\docType{data}
\title{Genetic risk score for adiponectin and effects on coronary artery
  disease.}
\description{
  Risk scores parameterised using GWAS meta-analysis results published
  by the ADIPOGen consortium (Dastani et al. PLoS Genetics 2012) for
  adiponectin, and effects on coronary artery disease from the
  CARDIoGRAM consortium and published by
  Dastani et al. (Atherosclerosis 2013).
}
\usage{data(adiponectin.scores)}
\format{
  A data frame suitable for use with other functions in this package,
  see \link{gtx.params}.

  The score using only 17 independent SNPs is called Adiponectin2013,
  and the \dQuote{coef} is defined only for these SNPs.  
  Effect size coefficients for adiponectin (for all 162 SNPs) and for
  CAD (where available) are taken from Tables S1 and S2 of Dastani et al. (2013).
}
\examples{
data(adiponectin.scores)
with(subset(adiponectin.scores, score == "Adiponectin2013"),
  grs.plot(beta_Adiponectin, beta_CAD, se_CAD, text = locus))
}
\source{
  For the publication by Dastani et al. from which these data were extracted see \url{http://dx.doi.org/10.1016/j.atherosclerosis.2013.03.034}.
}
\keyword{datasets}