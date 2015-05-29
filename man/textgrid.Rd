\name{textgrid}
\alias{textgrid}
\title{Add matrix of scaled text to a plot}
\description{
  Draws the strings given in a dataframe or matrix, with a font size
  scaled to ensure the text fits inside a specified region.
}
\usage{
textgrid(tmat, x0 = 0, y0 = 0, x1 = x0 + 1, y1 = y0 + 1, colsep = 3, rowsep = 1)
}
\arguments{
  \item{tmat}{dataframe or matrix of character strings}
  \item{x0}{coordinate of left edge of region}
  \item{y0}{coordinate of bottom edge of region}
  \item{x1}{coordinate of right edge of region}
  \item{y1}{coordinate of top edge of region}
  \item{colsep}{width (in multiples of width of M) of space between columns}
  \item{rowsep}{height (in multiples of height of M) of space between rows}
}
\value{
  Returns a list with the following elements.  \code{cex} provides the
  character scaling used.  \code{xpos} provides the coordinates for the
  left edge of each column.  \code{ypos} provides the coordinates for
  the top of each row.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
