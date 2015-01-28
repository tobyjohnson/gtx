\name{kmplot}
\alias{kmplot}
\alias{kmplot.character}
\alias{kmplot.Surv}
\title{Plot a survival variable by discrete categories}
\description{
  Plots a survival variable by discrete
  categories (such as genotype).
}
\usage{
\method{kmplot}{character}(object, x, data, ylab, xlab,
                           ylim, xlim, col, lty,
                           legend.location = "topright", atrisk = FALSE, ...)
\method{kmplot}{Surv}(object, x, data, ylab, xlab,
                      ylim, xlim, col, lty,
                      legend.location = "topright", atrisk = FALSE, ...)
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
  \item{legend.location}{Passed to first argument of call to
    \code{legend}}
  \item{atrisk}{Logical, whether to display at risk numbers}
  \item{...}{Other arguments}
}
\details{
  TO BE WRITTEN
}
\value{Returns an invisible null.  The plot is generated as a side effect.}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
\examples{
library(survival)
data(aoex1)
aoex1 <- within(aoex1, srvDays <- Surv(SRVDY, SRVCFLCD))
kmplot("srvDays", "rs123456", data = aoex1, xlab = "Days")
with(aoex1, kmplot(srvDays, rs123456, xlab = "Days"))
with(aoex1, kmplot(srvDays, STUDYID, xlab = "Days"))
}