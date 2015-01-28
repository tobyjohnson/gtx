\name{pgx.trtreat}
\alias{pgx.trtreat}
\title{Derive time in relation to treatment.}
\description{
  Function to derive a variable that classifies time points in relation
  to treatment, as pre-, on-, or post-therapy.
}
\usage{
pgx.trtreat(measure, exposure,
            mtime, winlen, 
            estart = "exstdtR", eend = "exendtR",
            subjid = "USUBJID")
}
\arguments{
  \item{measure}{A data frame containing a time series}
  \item{exposure}{A data frame containing therapy exposure records}
  \item{mtime}{Name of the variable in \code{measure} that describes
    time of measurement}
  \item{winlen}{Length of time beyond last exposure that is in
    on-therapy window}
  \item{estart}{Name of the variable in \code{exposure} that describes
    exposure start time}
  \item{eend}{Name of the variable in \code{exposure} that describes
    exposure end time}
  \item{subjid}{Name of the variable in \code{d}, \code{measure}, and
    \code{exposure}, that identifies individual subjects}
}
\value{
  A data frame of derived values for time in relation to treatment,
  which can be NA, \dQuote{Pre-therapy}, \dQuote{On-therapy}, or
  \dQuote{Post-therapy}.
}
\details{
  This function derives the time in relation to treatment, from clinical
  trial data consisting of: (i) time series; and (ii) therapy exposure records.

  HOW SHOULD WE HANDLE SUBJECTS WITH NO EXPOSURE RECORDS
}
\examples{
data(clinex) # like plain text export of SAS files from CDER/CDISC
             # compliant data
## format dates for R
clinex$LAB <- within(clinex$LAB, {
  lbdtR <- as.Date(LBDT, "\%d\%b\%Y")                
})
clinex$EXPOSURE <- within(clinex$EXPOSURE, {
  exstdtR <- as.Date(EXSTDT, "\%d\%b\%Y")                
  exendtR <- as.Date(EXENDT, "\%d\%b\%Y")                
})

table(pgx.trtreat(clinex$LAB, clinex$EXPOSURE,
            mtime = "lbdtR", winlen = 0), useNA="always")

}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
