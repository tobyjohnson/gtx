\name{mthfrex}
\alias{mthfrex}
\alias{mthfr.params}
\docType{data}
\title{Simulated example finemapping genotype and phenotype data.}
\description{
  A simulated example finemapping genotype and phenotype dataset,
  consisting of genotypes for 2000 individuals at 64 SNPs, and
  case/control and continuous phenotypes and covariates.
}
\usage{
data(mthfrex)
}
\format{
  \code{mthfrex} is a list suitable for use with other functions in this package,
  see \link{snpdata}.

  \code{mthfr.params} is a data frame suitable for use with other functions in this package,
  see \link{gtx.params}.
}
\details{
  This is a simulated dataset, whose sole purpose is to illustrate the
  use of functions in this package for multi-SNP regression analyses.  
  
  The dataset is provided so that the usage examples can actually be
  run, without burdening each example with many lines of code to
  generate an analysable dataset.
}
\source{
  The genotypes were simulated using Hudson's mksamples (ms) program,
  assuming an infinite sites standard neutral model with no
  recombination, with 1\% genotypes missing at random.  The simulated
  sites were crudely matched to real SNP names and hg18 map positions
  near the MTHFR locus (in an order to match rankwise with real allele frequencies). Some
  quasi-realistic case/control and continuous phenotypes and covariates
  were simulated conditional on the genotype data with some
  (unrealistically) large effect sizes (many OR>2).
}
