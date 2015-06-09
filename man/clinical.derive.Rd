\name{clinical.derive}
\alias{clinical.derive}
\title{Derive analysis variables from clinical data}
\description{
  Provides a convenience function
}
\usage{
clinical.derive(datalist, derivations, verbose = TRUE, only)
}
\arguments{
  \item{datalist}{A list of data frames of clinical data}
  \item{derivations}{A data frame of derivation rules}
  \item{verbose}{Whether to print progress messages}
  \item{only}{Character vector of targets to derive}
}
\details{
  Read the code or see the examples.

  The \code{derivations} argument must be a dataframe with the
  following columns (with the following information in each row):
  \dQuote{targets} (a name [or space delimited names] for the derived
  variable[s]); \dQuote{types} (the type [or space delimited types] for
  the derived variable[s]); \dQuote{deps} (a string with space delimited
  names of clinical datasets the derivation depends on); \dQuote{data} (a string
  with an R language statement to extract a data subset); \dQuote{fun}
  (a string with an R language statement of the derivation).
} 
\examples{
\dontrun{
clindata <- clinical.import("path/to/clinical/export/")
data(derivations.standard)
gxvars <- clinical.derive(clindata, derivations.standard)
summary(gxvars)
}
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
