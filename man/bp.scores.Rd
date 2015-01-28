\name{bp.scores}
\alias{bp.scores}
\docType{data}
\title{Genetic risk scores for blood pressure.}
\description{
  Risk scores parameterised using GWAS meta-analysis results for 29 SNPs
  published by the International Consortium for Blood Pressure GWAS
  (ICBP-GWAS Nature 2011).  There are scores for systolic blood
  pressure (SBP), diastolic blood pressure (DBP), and mean blood
  pressure (MBP).
}
\usage{
data(bp.scores)
}
\format{
  A data frame suitable for use with other functions in this package,
  see \link{gtx.params}.

  There are three scores are called SBP2011, DBP2011 and MBP2011.  Mean
  blood pressure is defined \eqn{(\mbox{SBP}+\mbox{DBP})/2}{(SBP+DBP)/2}
  and is \emph{not} the same as mean arterial pressure.  These are
  parameterised using the \dQuote{all data} effect size estimates from
  Table 1 and Supplementary Table 5 of the ICBP-GWAS manuscript.  There
  are also two scores called SBP2011A and DBP2011A.  There are
  parameterised using the \dQuote{unbiased} effect size estimates from
  Supplementary Appendix A of the ICBP-GWAS manuscript.  The
  coefficients are all in mmHg per coded allele dose after medication
  adjustment (which was +15mmHg SBP and +10mmHg DBP for medicated
  individuals).
}
\examples{
data(bp.scores)
head(subset(bp.scores, score == "SBP2011A"))
sbp <- subset(bp.scores, score == "SBP2011A")
dbp <- subset(bp.scores, score == "DBP2011A")
dbp <- dbp[match(sbp$locus, dbp$locus), ]
plot(sbp$coef/sign(sbp$coef), dbp$coef/sign(sbp$coef),
     xlim = c(0, max(sbp$coef/sign(sbp$coef))),
     ylim = c(0, max(dbp$coef/sign(sbp$coef))),
     xlab = "SBP effect [mmHg]", ylab = "DBP effect [mmHg]",
     las = 1)
}
\source{
  For the publication by ICBP-GWAS from which these data were
  extracted see \url{http://dx.doi.org/10.1038/nature10405}.
}
\keyword{datasets}