\name{t2dex}
\alias{t2dex}
\docType{data}
\title{Simulated example genotype and longitudinal phenotype data.}
\description{
  A simulated example genotype and phenotype dataset for genetic risk
  score analyses, consisting of genotypes for 2000 individuals at 31 SNPs, 
  baseline case/control and covariate status, and a survival phenotype.
}
\usage{
data(t2dex)
}
\format{
  \code{t2dex} is a list suitable for use with other functions in this package,
  see \link{snpdata}.
}
\details{
  This is a simulated dataset, whose sole purpose is to illustrate the
  use of functions in this package for genetic risk score analyses.  
  
  The dataset is provided so that the usage examples can actually be
  run, without burdening each example with many lines of code to
  generate an analysable dataset.
  
  The genotype data were simulated assuming exact Hardy-Weinberg and
  linkage equilibrium.  The phenotypes were simulated using a crude
  parametric model of a longitudinal study, intended to be used in
  regression models as if:
  \itemize{
    \item \code{t2dex$data$Age} is age at recruitment into the study.
    \item \code{t2dex$data$Overweight} 0/1 indicator for BMI>25 at
      recruitment.
    \item \code{t2dex$data$T2D} 0/1 indicator for prevalent T2D at
      recruitment.
    \item \code{t2dex$data$FollowupDays} days subject followed up for,
      until incident T2D event or followup stopped.
    \item \code{t2dex$data$FollowupT2D} 0/1 indicator for incident T2D
      event.
  }
}
\source{
  The genotypes were simulated  using allele frequencies from
  \link{t2d.scores}.
}
\examples{
data(t2dex)
summary(subset(t2dex$data, select = c("Age", "Overweight", "T2D",
                                      "FollowupDays", "FollowupT2D")))

library(survival)
plot(survfit(Surv(FollowupDays,FollowupT2D) ~ Overweight,
                 data = t2dex$data), col = c("green", "red"))

mycoxph <- coxph(Surv(FollowupDays,FollowupT2D) ~ Overweight,
                 data = t2dex$data) # fit null model
data(t2d.scores)
assoc1 <- grs.onesnp.apply(t2d.scores, mycoxph) # single SNP association
## risk score fit from single SNPs
unlist(grs.summary(t2d.scores$coef, assoc1$beta, assoc1$se,
                   n = length(residuals(mycoxph))))

## compare direct analysis of subject-specific data
t2dex <- grs.make.scores(t2d.scores, t2dex)
coxph(Surv(FollowupDays,FollowupT2D) ~ Overweight + T2D2010.score,
      data = t2dex$data)
}

