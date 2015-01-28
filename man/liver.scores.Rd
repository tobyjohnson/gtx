\name{liver.scores}
\alias{liver.scores}
\docType{data}
\title{Genetic risk scores for liver enzyme traits.}
\description{
  Risk scores parameterised using GWAS meta-analysis results published
  by Chambers et al. (Nature Genetics 2011) for serum ALT, ALP and
  GGT levels.
}
\usage{data(liver.scores)}
\format{
  A data frame suitable for use with other functions in this package,
  see \link{gtx.params}.

  The scores are called ALT2011, ALP2011 and GGT2011, for serum
  concentrations of alanine transaminase (ALT), alkaline phosphatase
  (ALP) and gamma-glutamyl transferase (GGT).  The first twelve columns
  of data were taken from Table 1 of Chambers et al. (2011).  Effect
  sizes (\code{coef}) are on log10(concentration) scale, from inverse
  variance weighted meta-analyses over 17 studies, with analysed sample
  sizes 45596, 56415 and 61089, for ALT, ALP and GGT respectively.  Two
  additional columns of data (\code{beta_replication} and
  \code{se_replication}) were taken from Supplementary Table 18 of
  Chambers et al. (2011), and give the effect size and SE (back
  calculated from the reported P-value) in an independent replication
  sample of size 8112.
}
\examples{
data(liver.scores)
with(subset(liver.scores, score == "ALP2011"), {
  plot(pmin(coded.freq, 1 - coded.freq), coef,
       xlab = "MAF", ylab = "ALP effect",
       xlim = c(0, 0.5), ylim = c(0, 10), las = 1)
  text(pmin(coded.freq, 1 - coded.freq), coef, name,
       pos = 3, cex = 0.5)})
## compare unweighted score used by Chambers et al. 2011
## with weighted score:
with(subset(liver.scores, score == "ALP2011"),
     rbind(unlist(grs.summary(rep(1, length(coef)), beta_replication,
                              se_replication, n = 8112)),
           unlist(grs.summary(coef, beta_replication,
                              se_replication, n = 8112))))
}
\source{
  For the publication by Chambers et al. from which these data were extracted see \url{http://dx.doi.org/10.1038/ng.970}.
}
\keyword{datasets}