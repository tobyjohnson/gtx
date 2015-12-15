\name{kmplot}
\alias{kmplot}
\alias{kmplot.character}
\alias{kmplot.Surv}
\title{Plot survival functions by discrete categories}
\description{
  Plots survival functions by discrete
  categories (such as genotype).
}
\usage{
\method{kmplot}{character}(object, x, data, ylab, xlab,
                           ylim, xlim, col, lty,
                           legend.location = "topright", legend.cex = 1,
                           legend.median = FALSE, digits = 1,
                           atrisk = FALSE, ...)
\method{kmplot}{Surv}(object, x, data, ylab, xlab,
                      ylim, xlim, col, lty,
                      legend.location = "topright", legend.cex = 1,
                      legend.median = FALSE, digits = 1,
                      atrisk = FALSE, ...)
}
\arguments{
  \item{object}{The name or values of the survival variable}
  \item{x}{The name or values of the discrete category variable}
  \item{data}{A data frame containing values}
  \item{ylab}{A label for the y axis}
  \item{xlab}{A label for the x axis}
  \item{ylim}{Range for the y-axis}
  \item{xlim}{Range for the x-axis}
  \item{col}{A vector of colours}
  \item{lty}{A vector of line types}
  \item{legend.location}{Location for legend}
  \item{legend.cex}{Character expansion (size) for legend}
  \item{legend.median}{Whether to calculate and include median survival
    times in legend}
  \item{digits}{Number of digits for legend data}
  \item{atrisk}{Logical, whether to display at risk numbers}
  \item{...}{Other arguments}
}
\details{
  This function calculates Kaplan--Meier estimates of survival functions
  and draws an annotated plot.  The arguments \code{object} and \code{x}
  specify the survival data and the levels of a discrete factor, either
  directly (as variables of class \code{Surv} and \code{factor}
  respectively), or as names of variables in the data frame \code{data}.

  Several options for annotating the plot are supported, which provide
  improvements over the default \code{print.survfit} method.  A legend
  with the number of subjects in each category is automatically added,
  this may be disabled by setting \code{legend.location="none"}.  The
  legend may be automatically annotated with the median survival times
  and 95\% confidence limits.  Numbers of subjects at risk, at each tick
  on the x-axis, may be added.
}
\value{Returns an invisible null.  The plot is generated as a side effect.}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
\examples{
library(survival)
data(aoex1)
aoex1 <- within(aoex1, srvDays <- Surv(SRVDY, SRVCFLCD))
## show two different ways of calling, either:
with(aoex1, kmplot(srvDays, rs123456, xlab = "Days"))
## or:
kmplot("srvDays", "rs123456", data = aoex1, xlab = "Days")
## show additional annotation
kmplot("srvDays", "rs123456", data = aoex1, xlab = "Days",
       legend.median = TRUE, atrisk = TRUE)
## show correct behaviour when one level has entirely missing survival data
with(aoex1, kmplot(srvDays, STUDYID, xlab = "Days"))
}