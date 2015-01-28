\name{grs.onesnp.apply}
\alias{grs.onesnp.apply}
\title{Convenience tool to fit a series of single-SNP models.}
\description{
  This function takes as input a fitted model object, typically
  created by the standard R functions lm, glm or coxph,  and fitted to individual-level genotype and
  phenotype data.  The model is augmented by adding a term for each
  single SNP in turn and the refitted.  This provides the single SNP
  summary association statistics needed to apply the summry statistic
  methods.
}
\usage{
grs.onesnp.apply(params, object, coeff.extract.fun = coeff.extract)
}
\arguments{
  \item{params}{a data frame, see \code{\link{gtx.params}}.}
  \item{object}{a fitted model object of class lm, glm, or coxph}
  \item{coeff.extract.fun}{function that extracts Estimate and Std.Err
    from fitted model objects, see \code{\link{coeff.extract}}.}
}
\details{
  By default this uses \code{\link{coeff.extract}} to detect the class of object and handles
  coefficient extraction appropriately.  Supply your own function if you
  have an object that works with update (supply your own update too).
}
\examples{
library(survival)
data(t2d.scores)
data(t2dex)

mycoxph <- coxph(Surv(FollowupDays,FollowupT2D) ~ Overweight,
                 data = t2dex$data) # fit null model
assoc1 <- grs.onesnp.apply(t2d.scores, mycoxph) # single SNP association
## risk score fit from single SNPs
unlist(grs.summary(t2d.scores$coef, assoc1$beta, assoc1$se,
                   n = length(residuals(mycoxph))))

## compare direct analysis of subject-specific data
t2dex <- grs.make.scores(t2d.scores, t2dex)
coxph(Surv(FollowupDays,FollowupT2D) ~ Overweight + T2D2010.score,
      data = t2dex$data)
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
