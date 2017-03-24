\name{coloc.fast}
\alias{coloc.fast}
\title{Calculate colocalisation probabilities.}
\description{
  Calculates posterior probabilities for colocalisation of
  two association signals, using the approximate Bayes factor
  method of Giambartolomei et al.
}
\usage{
coloc.fast(data, rounded = 6, priorsd1 = 1, priorsd2 = 1, priorc1 = 1e-4, priorc2 = 1e-4, priorc12 = 1e-5, beta1 = 'beta_1', se1 = 'se_1', beta2 = 'beta_2', se2 = 'se_2')
}
\arguments{
  \item{data}{Data frame of association results.}
  \item{rounded}{Number of digits to round posterior probabilities, or NA for no rounding.}
  \item{priorsd1}{Standard deviation of the prior
    on true effect sizes for phenotype 1.}
  \item{priorsd2}{Standard deviation of the prior
    on true effect sizes for phenotype 2.}
  \item{priorc1}{Prior on a given variant being causal for phenotype 1.}
  \item{priorc2}{Prior on a given variant being causal for phenotype 2.}
  \item{priorc12}{Prior on a given variant being causal for phenotypes 1 and 2.}
  \item{beta1}{Name of effect size for phenotype 1 in data.}
  \item{se1}{Name of standard error for phenotype 1 in data.}
  \item{beta2}{Name of effect size for phenotype 2 in data.}
  \item{se2}{Name of standard error for phenotype 2 in data.}
}
\value{
  A data frame of posterior model probabilities.
}
\details{
  To be described.  Equivalent to Giambartolomei method but uses a
  single framework for quantitative and binary phenotypes, and allows
  user control of the within-model priors on effect sizes (to which the
  results can be surprisingly sensitive).
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
