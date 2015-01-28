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
