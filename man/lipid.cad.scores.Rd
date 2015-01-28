\name{lipid.cad.scores}
\alias{lipid.cad.scores}
\docType{data}
\title{Genetic risk scores for serum lipid levels and coronary artery disease outcome}
\description{
  Risk scores parameterised using GWAS meta-analysis results published
  by Waterworth et al. (Arteriosclerosis, Thrombosis, and Vascular
  Biology 2010).  There are scores for three different lipid level
  phenotypes (High Density Lipoprotein HDL; Low Density Lipoprotein LDL;
  Triglycerides TG), all parameterised using SNPs reported in previous
  studies, with effect sizes estimated using the data of Waterworth et
  al. (2010 Table 1).  For each SNP the outcome, estimated effect on
  risk of coronary artery disease (CAD), from a case/control collection
  of 9633 cases and 38684 controls, was reported by Waterworth et
  al. (2010 Table 4).
}
\usage{lipid.cad.scores}
\format{
  A data frame suitable for use with other functions in this package.
  
  These scores are called HDL, LDL and TG, and have 
  coefficients in change in natural log transformed lipid level per coded allele.

  The coefficients for outcome effect are in natural log odds ratio per
  coded allele.
}
\source{
  For the publication by Waterworth et al. from which these
  data were extracted see \url{http://dx.doi.org/10.1161/ATVBAHA.109.201020}.
  }
\keyword{datasets}