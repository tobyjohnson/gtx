\name{waterfall}
\alias{waterfall}
\alias{waterfall.character}
\alias{waterfall.numeric}
\title{Waterfall plot coloured by discrete categories}
\description{
  Generates a waterfall plot, which summarises the distribution of a continuous variable (such as a quantitative trait) by discrete
  categories (such as genotype).
}
\usage{
\method{waterfall}{character}(object, x, data, ylab, ylim, bylevel, col, ...)
\method{waterfall}{numeric}(object, x, data, ylab, ylim, bylevel, col, ...)
}
\arguments{
  \item{object}{The name or values of the continuous variable}
  \item{x}{The name or values of the discrete category variable}
  \item{data}{A data frame containing values}
  \item{ylab}{A label for the y axis}
  \item{ylim}{Limits for the y axis}
  \item{bylevel}{Whether to group subjects by levels of x}
  \item{col}{A vector of colours}
  \item{...}{Other arguments}
}
\details{
  Waterfall plots are widely used in clinical oncology, to display the
  distribution of a specific quantitative trait:  best on-study change
  in tumour burden, relative to baseline.  The plot is typically
  coloured by levels of a factor, such as treatment arm or categorical
  best response.  In principle the plot could
  be used for any quantitative trait, but the display is most visually
  appealing when the trait takes both positive and negative values.
}
\value{Returns an invisible null.  The plot is generated as a side effect.}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
\examples{
data(aoex1)
par(mfrow = c(1, 2))
with(aoex1, waterfall(bmi-mean(bmi, na.rm = TRUE), rs123456, bylevel = TRUE))
with(aoex1, waterfall(bmi-mean(bmi, na.rm = TRUE), rs123456, bylevel = FALSE))
}