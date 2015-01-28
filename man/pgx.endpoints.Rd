\name{pgx.endpoints}
\alias{pgx.endpoints}
\title{Derive pharmacogenetic (PGx) safety endpoints from longitudinal data.}
\description{
  Function to derive safety endpoints for continuous, ordinal, binary,
  and time-to-event analyses, from longitudinal data typical of a clinical trial.
}
\usage{
pgx.endpoints(d, measure, exposure,
           mtime, mvalue, threshval,
           msign = c("greater", "less"), mstrict = c(TRUE, FALSE), 
           estart = "exstdtR", eend = "exendtR",
           dby = NULL, subjid = "USUBJID")
}
\arguments{
  \item{d}{A data frame containing subject identifiers, and an optional
    stratifying variable}
  \item{measure}{A data frame containing time series trait measurements}
  \item{exposure}{A data frame containing therapy exposure records}
  \item{mtime}{Name of the variable in \code{measure} that describes
    time of measurement}
  \item{mvalue}{Name of the variable in \code{measure} that describes
    value of measurement}
  \item{threshval}{Critical value, for defining events for
    \dQuote{sufficient exposure} to treatment and for survival analysis}
  \item{msign}{Indicator how \code{mvalue} is to be compared to \code{threshval}}
  \item{mstrict}{Indicator how \code{mvalue} is to be compared to \code{threshval}}
  \item{estart}{Name of the variable in \code{exposure} that describes
    exposure start time}
  \item{eend}{Name of the variable in \code{exposure} that describes
    exposure end time}
  \item{dby}{Name of the variable in \code{d} to stratify on for
    computing median exposure times}
  \item{subjid}{Name of the variable in \code{d}, \code{measure}, and
    \code{exposure}, that identifies individual subjects}
}
\value{
  A data frame of derived variables.
}
\details{
  This function derives safety endpoint variables that can be used for
  continuous, ordinal, binary, and time-to-event analyses, from clinical
  trial data consisting of: (i) time series trait measurements,
  typically indicating the presence or severity of symptoms for the
  safety event; and (ii) therapy exposure records.
  
  TO DO:  Allow elapsed time, cumulative dose
}
\examples{
library(survival)
data(clinex) # like plain text export of SAS files from CDER/CDISC
             # compliant data
## format dates for R
clinex$LAB <- within(clinex$LAB, {
  lbdtR <- as.Date(LBDT, "\%d\%b\%Y")                
  grade <- grade(LBSTRESN, c(1, 5, 10), strict = TRUE)
})
clinex$EXPOSURE <- within(clinex$EXPOSURE, {
  exstdtR <- as.Date(EXSTDT, "\%d\%b\%Y")                
  exendtR <- as.Date(EXENDT, "\%d\%b\%Y")                
})

pgx <- clinex$DEMO
pgx <- merge(pgx,
             pgx.endpoints(clinex$DEMO, clinex$LAB, clinex$EXPOSURE,
                        mtime = "lbdtR", mvalue = "LBSTRESN", threshval = 5,
                        msign = "greater", mstrict = TRUE))
pgx <- merge(pgx,
             pgx.endpoints(clinex$DEMO, clinex$LAB, clinex$EXPOSURE,
                        mtime = "lbdtR", mvalue = "grade", threshval = 2,
                        msign = "greater", mstrict = FALSE))
pgx <- within(pgx, {
  LBSTRESN.binary <- ifelse(!LBSTRESN.exposed, NA, 
                            ifelse(LBSTRESN.event == 1, 1,
                                   ifelse(LBSTRESN.max <= 1, 0, NA)))
  ATRTGRP <- as.factor(ATRTGRP)
})

summary(lm(LBSTRESN.max ~ ATRTGRP, data = pgx, subset = LBSTRESN.exposed))
summary(glm(LBSTRESN.binary ~ ATRTGRP, family = binomial, data = pgx))
summary(coxph(Surv(LBSTRESN.time, LBSTRESN.event) ~ ATRTGRP, data = pgx))
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
