\name{Surv2}
\alias{Surv2}
\title{Construct survival object from event and would-be-censored times}
\description{
  Convenience function to construct survival object from time to event
  (NA if event did not occur) and would-be-censored time (used if event
  did not occur).
}
\usage{
Surv2(tevent, tcensor)
}
\arguments{
  \item{tevent}{Vector of times until event, with NA indicating event
    did not occur}
  \item{tcensor}{Vector of would-be-censored times at which censoring
    would occur if event did not occur}
}
\details{
  Right censored survival objects are normally constructed using \code{\link{Surv}},
  specifying the time (until either event or censoring) and an indicator
  variable (censored or event).

  In some applications it may be more natural to construct the same
  survival object using two times: the time until the event (may be NA
  if the event never occurs), and the \dQuote{would-be-censored} time,
  which is the time at which censoring would occur if the event has not
  occured first.  \code{Surv2} provides a convenience wrapper for this
  method of construction.

  The arguments \code{tevent} and \code{tcensor} should be numeric
  vectors of equal length.  Elementwise: An event is determined to have
  occurred if \code{tevent} is nonmissing and less than or equal to
  \code{tcensor}.  The survival time is \code{tevent} if an event has
  occurred, and \code{tcensor} otherwise.  A survival object is created
  with this time and event status with a call to \code{\link{Surv}}.
}
\value{
  An object of class \code{\link{Surv}}.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
\seealso{
  \code{\link{Surv}}
}
\examples{
library(survival)
n <- 500 # sample size
tevent <- ifelse(rbinom(n, 1, 0.9) == 1, rexp(n, 0.1), NA)  # time until event
tcensor <- runif(n, 0, 5) # length of followup
mys <- Surv2(tevent, tcensor)
stopifnot(is.Surv(mys))
stopifnot(!any(is.na(mys)))
## check definition
stopifnot(identical(mys[ , "status"],
                    ifelse(!is.na(tevent) & tevent <= tcensor, 1, 0)))
stopifnot(all(mys[is.na(tevent), "status"] == 0))
}
