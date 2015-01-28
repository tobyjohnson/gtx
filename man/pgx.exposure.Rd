\name{pgx.exposure}
\alias{pgx.exposure}
\title{Derive variables about treatment exposure from longitudinal data.}
\description{
  Function to derive varianles about treatment exposure, from longitudinal data typical of a clinical trial.
}
\usage{
pgx.exposure(d, exposure,
             estart = "exstdtR", eend = "exendtR",
             subjid = "USUBJID")
}
\arguments{
  \item{d}{A data frame containing subject identifiers}
  \item{exposure}{A data frame containing therapy exposure records}
  \item{estart}{Name of the variable in \code{exposure} that describes
    exposure start time}
  \item{eend}{Name of the variable in \code{exposure} that describes
    exposure end time}
  \item{subjid}{Name of the variable in \code{d} and
    \code{exposure}, that identifies individual subjects}
}
\value{
  A data frame of derived variables.
}
\details{
  This function derives total time on treatment, using
  \code{\link{mincover}}.
  
  TO DO:  Calculate first and last days, elapsed time, cumulative dose
}
\examples{
data(clinex) # like plain text export of SAS files
clinex$EXPOSURE <- within(clinex$EXPOSURE, {
  exstdtR <- as.Date(EXSTDT, "\%d\%b\%Y")                
  exendtR <- as.Date(EXENDT, "\%d\%b\%Y")                
})

pgx <- clinex$DEMO
pgx <- merge(pgx,
             pgx.exposure(clinex$DEMO,
                          subset(clinex$EXPOSURE, ACTTRT == "Active")))
head(pgx)
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
