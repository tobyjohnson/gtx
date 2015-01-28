\name{derivation.add}
\alias{derivation.add}
\title{Add a derivation}
\description{
  Provides a consistent interface for adding
  a derivation.
}
\usage{
derivation.add(derivations, targets, types, deps, data, fun, aept.list, verbose = TRUE)
}
\arguments{
  \item{derivations}{A data frame of derivation rules}
  \item{targets}{Names of derived variables, whitespace delimited}
  \item{types}{Types of derived variables, in sense of \code{storage.mode(target)}, whitespace delimited}
  \item{deps}{Names of dependency data frames, whitespace delimited}
  \item{data}{Parseable text of code to extract data}
  \item{fun}{Parseable text of code to derive variables}
  \item{aept.list}{Not used, for guessing derivation rules}
  \item{verbose}{Whether to print progress messages}
}
\details{
  Read the code or see the examples.
} 
\examples{
\dontrun{
clindata <- derivation.add("path/to/clinical/export/")
data(derivations.standard)
derivations.custom <-
  derivation.add(derivations.standard,
                 targets = 'vsanal.SBP.mean4wk',
                 types = 'double',
                 deps = 'vitals',
                 data = 'subset(vsanal, VSTEST == "Systolic BP~(mmHg)" &
                                        ATTYPE == "On-therapy" &
                                        VSACTDY >= 1 & VSACTDY <= 35)',
                 fun = 'safe(mean, VSSTRESN)')
gxvars <- clinical.derive(clindata, derivations.custom)
summary(gxvars)
}
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
