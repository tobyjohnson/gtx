\name{grs.plot}
\alias{grs.plot}
\title{Diagnostic plot for genetic risk score calculation from summary statistics.}
\description{
  Each SNP is plotted by coefficient in the risk score (x axis) versus
  estimated effect size for trait of interest in the testing dataset (y axis).  A solid red
  line shows the effect size estimate for the risk score on the trait of
  interest in the testing dataset.
}
\usage{
grs.plot(w, b, s, text = NULL, textpos = NULL, textcex = 0.5,
         alpha = 0.05)
}
\arguments{
  \item{w}{coefficients for the risk score.}
  \item{b}{aligned beta coefficients in the testing dataset.}
  \item{s}{standard errors for \code{b}.}
  \item{text}{optional labels for SNPs}
  \item{textpos}{position for labels for SNPs, as \code{pos} argument to \command{text()}.}
  \item{textcex}{size for lables for SNPs, as \code{cex} argument to \command{text()}.}
  \item{alpha}{determines size of confidence limits plotted.}
}
\examples{
data(t2dex)
library(survival)
mycoxph <- coxph(Surv(FollowupDays,FollowupT2D) ~ Overweight,
                 data = t2dex$data) # fit null model
data(t2d.scores)
assoc1 <- grs.onesnp.apply(t2d.scores, mycoxph) # single SNP association
## risk score fit from single SNPs
grs.plot(t2d.scores$coef, assoc1$beta, assoc1$se, t2d.scores$name)
title(xlab = "risk score weight", ylab = "estimated effect size")
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
