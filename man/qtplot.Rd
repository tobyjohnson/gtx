\name{qtplot}
\alias{qtplot}
\alias{qtplot.character}
\alias{qtplot.numeric}
\title{Plot a continuous variable by discrete categories}
\description{
  Plots a continuous variable (such as a quantitative trait) by discrete
  categories (such as genotype).
}
\usage{
\method{qtplot}{character}(object, x, data, ylab, xlab, ylim, col, points.method = "beeswarm", boxplot.layer = "over", ...)
\method{qtplot}{numeric}(object, x, data, ylab, xlab, ylim, col, points.method = "beeswarm", boxplot.layer = "over", ...)
}
\arguments{
  \item{object}{The name or values of the continuous variable}
  \item{x}{The name or values of the discrete category variable}
  \item{data}{A data frame containing values}
  \item{ylab}{A label for the y axis}
  \item{xlab}{A label for the x axis}
  \item{ylim}{y axis range}
  \item{points.method}{Name of method for plotting individual datapoints, either
    \dQuote{stripchart} or \dQuote{beeswarm}}
  \item{boxplot.layer}{Layer to draw the boxplot relative to individual
    datapoints, either \dQuote{over}, \dQuote{under}, or \dQuote{none}}
  \item{col}{A vector of colours}
  \item{...}{Other arguments, ignored}
}
\details{
  Generates a plot displaying numerical data (usually continuous) by
  levels of a factor.

  When displaying such data, it is preferable to avoid the loss of
  information inherent in a classical boxplot, and to display all the
  individual datapoints if at all possible.  However, the boxplot
  summary remains a desirable adjunct.

  The argument \code{points.method} specifies how individual datapoints
  are drawn.  The default, \dQuote{beeswarm}, is deterministic and avoids
  overlap of individual datapoints unless absolutely necessary, but may
  be slow for large numbers of datapoints.  The alternative,
  \dQuote{stripchart}, is fast but is not deterministic and tends to
  overlap datapoints.

  When datapoints are overlapping, it may be possible to visualise this
  using transparency in \code{col}.
}
\value{Returns an invisible null.  The plot is generated as a side effect.}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
\examples{
data(aoex1)
qtplot("bmi", "rs123456", data = aoex1)
with(aoex1, qtplot(bmi, rs123456))
x <- sample(1:2, 1000, replace = TRUE)
y <- rnorm(1000, x, .02)
qtplot(y, x, points.method = "stripchart")
\dontrun{
qtplot(y, x) # takes a long time
}
}