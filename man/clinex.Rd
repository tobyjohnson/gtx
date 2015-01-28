\name{clinex}
\alias{clinex}
\docType{data}
\title{Simulated example clinical trial data.}
\description{
  A simulated example of clinical trial data, with longitudinal
  laboratory measures and treatment exposure.
}
\usage{
data(clinex)
}
\format{
  \code{clinex} is a list containing three data frames, \code{DEMO},
  \code{LAB} and \code{EXPOSURE}.

  \code{DEMO} is a data frame with 200 rows, one per subject.  The
  columns are the subject identifier \dQuote{USUBJID} and actual
  treatment indicator \dQuote{ATRTGRP}.

  \code{LAB} is a data frame with 1222 rows, describing between 2 and 13
  measurements on each subject.  The columns are the subject identifier
  \dQuote{USUBJID}, a visit identifier \dQuote{VISIT}, the visit date
  \dQuote{LBDT}, the visit days relative to baseline for each subject
  \dQuote{LBACTDY}, and the laboratory measured value \dQuote{LBSTRESN}.

  \code{EXPOSURE} is a data frame with 212 rows, with between 1 and 3
  records on each subject describing treatment exposure.  The columns
  are the subject identifier \dQuote{USUBJID}, a treatment identifier
  \dQuote{ACTTRT}, the exposure start and end dates \dQuote{EXSTDT} and
  \dQuote{EXENDT}, and the exposure start and end days relative to
  baseline for each subject \dQuote{EXACTSDY} and \dQuote{EXACTEDY}.
}
\details{
  This is a simulated dataset, whose sole purpose is to illustrate the
  use of functions in this package.
  
  The dataset is provided so that the usage examples can actually be
  run, without burdening each example with many lines of code to
  generate an analysable dataset.
}
