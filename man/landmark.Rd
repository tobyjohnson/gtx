\name{landmark}
\alias{landmark}
\title{Cuts survival object using (subject-specific) landmark times}
\description{
  Converts to missing, the survival data for subjects with events of
  censoring before the specified landmark times.
}
\usage{
landmark(s, tstart)
}
\arguments{
  \item{s}{Survival object}
  \item{tstart}{Vector of times at which landmarks occur}
}
\details{
}
\value{
  An object of class \code{\link{Surv}}.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
\seealso{
  \code{\link{Surv}}, dynpred package
}
\examples{
library(survival)
n <- 1000 # sample size
mys <- Surv2(rexp(n, 0.1), runif(n, 0, 5)) # random survival
tlm <- runif(n, 0, 1) # time of landmark
x <- ifelse(mys[ , "time"] >= tlm, rbinom(n, 3, .5), 0)
# x is uncorrelated with mys but zero for subjects with event or
# censor before tlm==tlandmark
drop1(coxph(mys~x), test = "Chisq") # biased
drop1(coxph(landmark(mys, tlm)~x), test = "Chisq") #unbiased
}
