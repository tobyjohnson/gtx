\name{magic.scores}
\alias{magic.scores}
\docType{data}
\title{Genetic risk scores for glucose/insulin traits.}
\description{
  Risk scores parameterised using GWAS meta-analysis results published
  by the MAGIC consortium (Dupuis et al. Nature Genetics 2010) for
  glucose and insulin traits.
}
\usage{data(magic.scores)}
\format{
  A data frame suitable for use with other functions in this package,
  see \link{gtx.params}.

  The scores are called FG2010 and FI2010, for fasting plasma glucose
  (FG, mmol/l) and fasting insulin (FI, pmol/l), and using the 16 and 2
  genome wide significant loci from Table 1 of Dupuis et al. (2010)
  respectively.  Effect size coefficients are taken from Table 2 of
  Dupuis et al. (2010) and hence use SNPs different from those reported
  in Table 1 at the SLC30A8 and TCF7L2 loci.

  In addition, for all SNPs in both risk scores, there are effect size
  estimates (beta) and standard errors (se) for FG, HOMA-B index of beta
  cell function, FI, HOMA-IR index of insulin resistance, glycated
  hemoglobin (HbA1c, percent), glucose 2 hours after oral glucose
  tolerance test (2hG, mmol/l), insulin 2 hours after oral glucose
  tolerance test (2hI, pmol/l), type 2 diabetes risk (T2D, ln odds,
  converted from odds ratio), body mass index (BMI, kg/m2), diastolic
  blood pressure (DBP, mmHg), systolic blood pressure (SBP, mmHg),
  hypertension risk (HTN, ln odds), and serum lipid concentrations (HDL,
  LDL, TC, TG), all taken from Table 2 and Table 3 of Dupuis et
  al. (2010).
}
\examples{
data(magic.scores)
with(subset(magic.scores, score == "FG2010"),
  grs.plot(beta_FG, beta_TG, se_TG, text = locus))
}
\source{
  For the publication by Dupuis et al. from which these data were extracted see \url{http://dx.doi.org/10.1038/ng.520}.
}
\keyword{datasets}