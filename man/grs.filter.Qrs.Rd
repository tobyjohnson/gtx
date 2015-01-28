\name{grs.filter.Qrs}
\alias{grs.filter.Qrs}
\title{Filter SNPs for inclusion in genetic risk score using
  heterogeneity test.}
\description{
  Performs a stepwise downward \dQuote{model selection} in which SNPs are
  iteratively removed from the risk score until the heterogeneity test
  is no longer significant at the specified threshold.
}
\usage{
grs.filter.Qrs(w, b, s, p.thresh = 0.05)
}
\arguments{
  \item{w}{coefficients for the risk score}
  \item{b}{aligned beta coefficients in the testing dataset}
  \item{s}{standard errors}
  \item{p.thresh}{P-value threshold}
}
\details{
  When there are m SNPs in the risk score, the heterogeneity
  test is an (m-1) d.f. LRT comparing the 1 d.f. risk score model
  against the unconstrained m d.f. model, as reported by the Qrs
  element of \code{\link{grs.summary}}.  At each iteration, if the
  risk score model has a significant heterogeneity test (at \code{p.thresh}), the 
  SNP is removed that gives the greatest decrease in the heterogeneity
  test statistic.
}
\examples{
  data(magic.scores)
  score1 <- subset(magic.scores, score == "FG2010")
  score1 <- within(score1, okay <- grs.filter.Qrs(coef, beta_TG, se_TG))
  with(score1, {grs.plot(coef, beta_TG, se_TG, locus);
    title(xlab = "FG effect", ylab = "TG effect")})
  with(score1, locus[!okay]) # loci removed
  with(subset(score1, okay), {grs.plot(coef, beta_TG, se_TG, locus);
    title(xlab = "FG effect", ylab = "TG effect")})
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
