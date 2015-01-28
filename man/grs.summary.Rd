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
  \dQuote{testing} dataset.

  The sample size argument \code{n} is required only to compute the
  (pseudo) variance explained in the testing dataset from the likelihood
  ratio test statistic.
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
Cite the paper by T. Johnson.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
