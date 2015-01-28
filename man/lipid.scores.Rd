\name{lipid.scores}
\alias{lipid.scores}
\docType{data}
\title{Genetic risk scores for serum lipid levels}
\description{
  Risk scores parameterised using GWAS meta-analysis results published
  by the Global Lipid Genetics Consortium (Teslovich et al. Nature
  2010).  There are scores for four different lipid level phenotypes
  (High Density Lipoprotein HDL; Low Density Lipoprotein LDL; Total
  Cholesterol TC; Triglycerides TG), all parameterised using genome wide
  significant SNPs reported in Supplementary Table 1 of Teslovich et
  al. (2010).
}
\usage{lipid.scores}
\format{
  A data frame suitable for use with other functions in this package.
  
  These scores are called HDL2010, LDL2010, TC2010 and TG2010, and have 
  coefficients in mg/dl per coded allele.

  Additional scores based on larger numbers of SNPs are available separately...
  LDL.prune.0.01 and there is an additional column in the data frame,
  p-value, that should be used to threshold on.  The coefficients are
  back-calculates from Z scores and HapMap allele frequency data.
}
\source{
  For the publication by Teslovich et al. from which these
  data were extracted see \url{http://dx.doi.org/10.1038/nature09270}.
  }
\keyword{datasets}