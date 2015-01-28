\name{trial.params}
\alias{trial.params}
\title{Calculate Parameters for a Simple Model of a Clinical Trial}
\description{
  Calculates parameters for a simple model of a clinical trial, with a
  constant rate recruitment period and then a followup period, with
  constant event rate per recruited live subject.
}
\usage{
trial.params(hazardRate = NA, surviveEnd = NA, timeTotal = NA, timeRecruit = NA)
}
\arguments{
  \item{hazardRate}{Instantaneous rate of events}
  \item{surviveEnd}{Proportion of subjects surviving at end of trial}
  \item{timeTotal}{Total duration of the trial}
  \item{timeRecruit}{Time of subject recruitment}
}
\details{
  This function assumes simplified model of a clinical trial (albeit
  less simplified than other models).  Events occur at a fixed
  instantaneous rate (\code{hazardRate}; \eqn{h}{h}) per recruited live
  subject and per unit time, and the expected surviving fraction at the
  end of the trial is \code{surviveEnd} (\eqn{S_\mathrm{End}}{sEnd}).
  The total duration of the trial (from start of recruitment until the
  time at which \eqn{S_\mathrm{End}}{sEnd} is evaulated, is
  \code{timeTotal} (\eqn{t_\mathrm{Total}}{tTotal}).  Recruitment is
  assumed to occur at a constant rate over a period 
  \code{timeRecruit} (\eqn{t_\mathrm{Recruit}}{tRecruit}).  The two
  extreme possibilities are \eqn{t_\mathrm{Recruit}=0}{tRecruit=0},
  where all subjects are recruited instantaneously at the start of the trial, and
  \eqn{t_\mathrm{Recruit}=t_\mathrm{Total}}{tRecruit=tTotal}, where
  recruitment is continuous over the entire duration of the trial.

  For specified values for any three of the four parameters, this
  function calculates the value for the missing paramater.

  The constant instantaneous event rate \eqn{h}{h} implies that for an
  exponential lifetime distribution, the mean lifetime is \eqn{1/h}{1/h}
  and the median lifetime is
  \eqn{\mathrm{qexp}(0.5)/h\simeq0.693/h}{qexp(0.5)/h ~= 0.693/h}.

  Note that with an exponential lifetime model, some parameter
  combinations cannot be solved numerically, e.g. calculating \eqn{h}{h}
  when \eqn{t_\mathrm{Total}=t_\mathrm{Recruit}}{tTotal=tRecruit} and
  \eqn{S_\mathrm{End}}{sEnd} very close to zero.
}
\value{
  A structure with missing values calculated.
}
\examples{
trial.params(hazardRate = c(0.1, 0.3, 0.5, NA),
            surviveEnd = c(NA, 0.4, 0.6, 0.8),
            timeTotal = c(4, NA, 3, 2),
            timeRecruit = c(2, 1.6, NA, 2))
## check missing return value is indeed not achievable:
trial.params(hazardRate = 0.5, surviveEnd = NA,
            timeTotal = 3, timeRecruit = c(0, 3))
## cannot solve all these cases numerically:
trial.params(NA, .00001, 1, c(.9, .95, .99, 1))
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}

