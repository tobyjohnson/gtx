\name{grs.summary}
\alias{grs.summary}
\title{Genetic risk score calculation from summary statistics.}
\description{
  Implements the summary statistic method of Johnson et al. for
  approximating the regression of a response variable onto an additive
  multi-SNP genetic risk score in a given testing dataset, using only single SNP association
  summary statistics.
}
\usage{
grs.summary(w, b, s, n)
}
\arguments{
  \item{w}{coefficients for the risk score.}
  \item{b}{aligned beta coefficients in the testing dataset, of same
    length as \code{w}.}
  \item{s}{standard errors for \code{b}, of same length as \code{w} and \code{b}.}
  \item{n}{sample size of testing dataset.}
}
\details{
  The risk score coefficients \code{w} are the \dQuote{weights} used to construct
  the risk score, for a set of SNPs, in chosen units per dose of the coded allele.
  Typically these are single SNP regression coefficients estimated in a
  \dQuote{discovery} dataset.

  The aligned beta coefficients \code{b} are regression coefficients for
  the response variable of interest, for the same set of SNPs and per
  dose of the same coded allele as used to define \code{w}.  Typically
  these are single SNP regression coefficients estimated in the
  \dQuote{testing} dataset.  The standard errors \code{s} are standard
  errors on \code{b}.

  In applications to causal inference, a common objective is to estimate
  the causal effect of an intermediate trait or biomarker, on a response
  variable or outcome.  In such applications, the \code{w} are the
  estimated effects on the intermediate trait or biomarker, and the
  \code{b} are estimated effects on the response variable or outcome,
  with standard errors \code{s}.
  
  The sample size argument \code{n} is required only to compute the
  (pseudo) variance explained in the testing dataset from the likelihood
  ratio test statistic.

  The method for calculating the regression of the response variable
  onto the risk score was first used for the work of the International
  Consortium for Blood Pressure Genome-Wide Association Studies (2011),
  and described in more detail in Dastani et al. (2012).  The method is
  exact for uncorrelated SNPs and a quadratic log-likelihood, the latter
  being obtained under a normal linear model, or under any regression
  model with a large sample size.
  
  The heterogeneity test is a test of whether the regression
  coefficients for the response variable are proportional to the risk
  score coefficients.  It is described in detail in the
  \dQuote{ashg2012} package vignette.  In applications to causal
  inference, firstly note that the heterogeneity test often lacks power,
  and hence a non-significant heterogeneity test is not evidence of
  clean instruments.  Secondly note that poor fit may be detected either
  when there are pleiotropic effects, or alternatively when one or more
  of the coefficients used to parameterise the risk score have been
  estimated imprecisely or with bias, and therefore a significant
  heterogeneity test is not necessarily evidence of unclean instruments.
  Nonetheless, a significant heterogeneity test may indicate that
  underlying assumptions should be subjected to extra scrutiny before
  any inference is made about causality.
}
\value{
  A named list with the following elements:  \code{m} is the number of
  SNPs used in the risk score.  \code{n} is the input sample size.
  \code{X2m} is the chi squared test statistic for an m d.f. test in the
  testing dataset (all
  SNPs have independent effects).  \code{R2m} is the (pseudo) variance
  explained by the m d.f. model in the testing dataset.  \code{ahat} is
  the estimated coefficent for regressing the response onto the m SNP
  risk score.  \code{aSE} is the standard error.  \code{X2rs} is the chi
  squared test statistic for a 1 d.f. test for the risk score in the
  testing dataset.  \code{R2rs} is the (pseudo) variance
  explained by the risk score model in the testing dataset.  \code{pval}
  is the P-value for the 1 d.f. test.  \code{Qrs} is the (m-1)
  d.f. heterogeneity test statistic.  \code{phet} is the heterogeneity
  test P-value.
}
\references{
  International Consortium for Blood Pressure Genome-Wide Association
  Studies (2011 Nature) \url{http://dx.doi.org/10.1038/nature10405}

  Dastani et al. (2012 PLoS Genetics)
  \url{http://dx.doi.org/10.1371/journal.pgen.1002607}

  Johnson (2012 ASHG poster) see \dQuote{ashg2012} package vignette.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
