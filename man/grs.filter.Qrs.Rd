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
  \strong{A critical bug in this function was corrected in gtx
  package version 0.0.10.  Output from previous versions should be
  treated as incorrect.  Even with the bug fixed, the function
  implements a method with no well defined statistical
  properties.}

  When there are m SNPs in the risk score, the heterogeneity
  test is an (m-1) d.f. LRT comparing the 1 d.f. risk score model
  against the unconstrained m d.f. model, as reported by the Qrs
  element of \code{\link{grs.summary}}.  At each iteration, if the
  risk score model has a significant heterogeneity test (at \code{p.thresh}), the 
  SNP is removed that gives the greatest decrease in the heterogeneity
  test statistic.

  No guarantee is given about the performance of this procedure.
  Intuitively, it is expected to work when the majority of SNPs in the
  risk score only affect the outcome with effects proportional to their
  weights in the score (including the possibility of zero effects), and
  a minority of SNPs affect the outcome with non-proportional effects.
  The working definition of \sQuote{minority} is unclear.

  When using a risk score to make causal inference, the application of
  this procedure (and indeed of any use of the heterogeneity test)
  should not be viewed as a replacement for detailed biological
  knowledge about the mechanisms of action of the causal genetic
  variants tagged by the SNPs used in the risk score.

  Compared with iteratively calculating the heterogeneity test using
  \code{\link{grs.summary}} and removing SNPs based on inspection of
  \code{\link{grs.plot}}, the \code{grs.filter.Qrs} procedure offers
  only convenience and (a modicum of) objectivity, and nothing more.

  Since this heterogeneity test was initially proposed, a number of
  alternative methods have been described that are likely to be
  superior, including the MR-Egger approach of:

  Bowden, J, Davey Smith, G, and Burgess, S (2015) Mendelian
  randomization with invalid instruments: Effect estimation and bias
  detection through Egger regression.  Int. J. Epidemiol. 44:512-525
  \url{https://dx.doi.org/10.1093/ije/dyv080}
}
\value{
  A logical vector of the same length as \code{w}, set TRUE or FALSE
  respectively for SNPs included or excluded at the end of the model
  selection procedure.
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
