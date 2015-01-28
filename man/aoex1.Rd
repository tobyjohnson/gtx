\name{aoex1}
\alias{aoex1}
\docType{data}
\title{Example dataset 1}
\description{
A small example dataset, simulated.
}
\usage{data(aoex1)}
\format{
  A data frame with 504 observations on the following 12 variables.
  \describe{
    \item{\code{STUDYID}}{a factor with levels \code{AO000001} \code{AO000002}}
    \item{\code{SUBJID}}{a numeric vector}
    \item{\code{USUBJID}}{a factor}
    \item{\code{TRTGRP}}{a factor with levels \code{DRUG A 10MG} \code{DRUG B 10MG} \code{DRUG B 50MG} \code{PLACEBO}}
    \item{\code{TRTGRPCD}}{a numeric vector}
    \item{\code{SRVDY}}{a numeric vector}
    \item{\code{SRVCFL}}{a factor with levels \code{Censored} \code{Died}}
    \item{\code{SRVCFLCD}}{a numeric vector}
    \item{\code{rs123456}}{a factor with levels \code{A A} \code{T A} \code{T T}}
    \item{\code{rs123456_A}}{a numeric vector}
    \item{\code{bmi}}{a numeric vector}
    \item{\code{AGE}}{a numeric vector}
    \item{\code{SEX}}{a factor with levels \code{Female} \code{Male}}
    \item{\code{SEXCD}}{a numeric vector}
  }
}
\examples{
data(aoex1)
}
\keyword{datasets}
