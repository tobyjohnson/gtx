\name{agtstats}
\alias{agtstats}
\docType{data}
\title{Hypertension association statistics for SNPs near the AGT gene.}
\description{
  Hypertension case/control association test statistics (effect size
  estimate and standard error after genomic control) for 21 single
  nucleotide polymorphisms (SNPs) near
  the angiotensinogen (AGT) gene, from a meta-analysis of 25118 subjects
  from 10 cohorts, genotyped using the Illumina HumanCVD BeadChip.
}
\usage{
data(agtstats)
}
\format{
  \code{agtstats} is a data frame with three informative columns:
  \code{dbSBNP} is the name of the SNP, \code{beta} is the meta-analysis
  effect size estimate, and \code{se.GC} is the meta-analysis standard
  error after genomic control (GC).
}
\source{
  Taken directly from Supplemental Table S13 of Johnson et al. (2011),
  see \url{http://dx.doi.org/10.1016/j.ajhg.2011.10.013}.
}
