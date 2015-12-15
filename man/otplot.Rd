\name{otplot}
\alias{otplot}
\alias{otplot.character}
\alias{otplot.factor}
\alias{otplot.numeric}
\title{Plot a factor variable by discrete categories}
\description{
  Plots a factor variable (such as binary yes/no, case/control, or an ordinal trait) by discrete
  categories (such as genotype).
}
\usage{
\method{otplot}{character}(object, x, data, ylab, xlab, col, col.mix = c("white", "black"), style = "percent", ycrop = FALSE, yext = 1.25)
\method{otplot}{factor}(object, x, data, ylab, xlab, col, col.mix = c("white", "black"), style = "percent", ycrop = FALSE, yext = 1.25)
\method{otplot}{numeric}(object, x, data, ylab, xlab, col, col.mix = c("white", "black"), style = "percent", ycrop = FALSE, yext = 1.25)
}
\arguments{
  \item{object}{The name or values of the factor variable}
  \item{x}{The name or values of the discrete category variable}
  \item{data}{A data frame containing values}
  \item{ylab}{A label for the y axis}
  \item{xlab}{A label for the x axis}
  \item{col}{A vector of colours for levels of \code{x}}
  \item{col.mix}{A pair of colours to mix with \code{col} for levels of \code{object}}
  \item{style}{A character vector, see details}
  \item{ycrop}{Logical, whether to crop y axis to data range}
  \item{yext}{Extension of y axis to accomodate legend}
}
\details{
  Generates a plot displaying factor data by levels of a factor.

  Currently, all supported styles are variations on a barplot, with bar
  heights showing either the count, fraction, or percentage of the
  response factor \code{object}.

  Colour mixing is not currently as flexible as it could be.

  Legend placement heuristics could be improved.
}
\value{Returns an invisible null.  The plot is generated as a side effect.}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
\examples{
  data(aoex1)
  aoex1 <- within(aoex1, {
                  bmicat <- cut(aoex1$bmi, 
                                c(-Inf, 25, 30, Inf))})
  ## Simple example
  with(aoex1, otplot(bmicat, rs123456, 
         ylab = "BMI category", xlab = "rs123456 genotype"))
  ## Illustration of different plot styles,
  ## using alternative function calling style
  oldpar <- par(mfcol = c(2, 3))
  for (style1 in c("percent", "percentall", "fraction", "fractionall", "count")) {
    otplot("bmicat", "rs123456", aoex1,
           ylab = "BMI category", xlab = "rs123456 genotype",
           style = style1)
    title(main = paste('style = "', style1, '"', sep = ""))
  }
  par(oldpar)
  ## Numerics are coerced to factors
  with(aoex1, otplot(rs123456_A, rs123456, 
         ylab = "Numeric genotype", xlab = "rs123456 genotype"))
  ## Illustration of alternative colour scheme      
  otplot("bmicat", "rs123456", aoex1,
         ylab = "BMI category", xlab = "rs123456 genotype",
         col = rgb(.3, 0, .15),
         col.mix = c(rgb(1, 1, .9), rgb(0, 0, .4)))
}
