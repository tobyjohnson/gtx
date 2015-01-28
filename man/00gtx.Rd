\name{gtx}
\alias{gtx}
\title{Genetics ToolboX}
\description{
  This package implements assorted tools for genetic association
  analyses, which is viewed as being entirely an exercise in regressing
  a (possibly multivariate) phenotypic \dQuote{response variable} onto
  one or more \dQuote{explanatory variables} that include genetic
  variables.

  The focus of this package is on implementing (either exactly or
  approximately) regression analyses using summary statistics instead of
  using subject-specific genotype and phenotype data.  So far, functions
  exist to support two applications detailed below: Multi-SNP risk score
  analyses; and multi-SNP conditional regression analyses.  In addition
  there are helper functions for reading and manipulating
  subject-specific genotype and phenotype data, which provide a platform
  for calculating the necessary summary statistics.
  
  The first application is multi-SNP risk score analyses, and the main
  functions provided for analysing summary statistics are
  \code{\link{grs.summary}}, \code{\link{grs.plot}} and
  \code{\link{grs.filter.Qrs}}.  The summary statistics necessary for
  these analyses are single SNP association statistics, which can be
  calculated using a wide variety of existing tools for GWAS analysis
  and meta-analysis.
  
  The second application is multi-SNP conditional or multiple regression
  analyses.  The main functions provided for performing multiple regression using
  summary statistics are \code{\link{combine.moments2}},
  \code{\link{est.moments2}}, \code{\link{lm.moments2}} and \code{\link{stepup.moments2}}.  The summary
  statistics necessary for these analyses can be calculated from
  subject-specific genotype and phenotype data, using the function
  \code{\link{make.moments2}}.
  
  The helper functions for reading and manipulating subject-specific
  genotype and phenotype data provide a convenient interface from R to
  genotype data exported from PLINK, and imputed genotype data generated
  by MACH or by IMPUTE.  The main functions provided are
  \code{\link{read.snpdata.plink}}, \code{\link{read.snpdata.mach}} and
  \code{\link{read.snpdata.impute}}.

}
\author{Toby Johnson \email{T.Johnson@qmul.ac.uk}}
