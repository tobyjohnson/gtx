\name{textgrid}
\alias{textgrid}
\title{Add matrix of scaled text to a plot}
\description{
  Draws the strings given in a dataframe or matrix, with a font size
  scaled to ensure the text fits inside a specified region.
}
\usage{
textgrid(tmat, x0 = 0, y0 = 0, x1 = x0 + 1, y1 = y0 + 1, xalign = "l",
xjust = 0.5, yjust = 0, colsep = 3, rowsep = 1, include.rownames,
include.colnames, draw = TRUE)
}
\arguments{
  \item{tmat}{dataframe or matrix of character strings}
  \item{x0}{coordinate of left edge of region}
  \item{y0}{coordinate of bottom edge of region}
  \item{x1}{coordinate of right edge of region}
  \item{y1}{coordinate of top edge of region}
  \item{xalign}{character vector of \sQuote{l}, \sQuote{c}, or \sQuote{r} specifying
    alignment of text within columns}
  \item{xjust}{justification within \code{x0}-\code{x1}}
  \item{yjust}{justification within \code{y0}-\code{y1}}
  \item{colsep}{width (in multiples of width of M) of space between columns}
  \item{rowsep}{height (in multiples of height of M) of space between rows}
  \item{include.rownames}{whether to draw rownames of \code{tmat}}
  \item{include.colnames}{whether to draw colnames of \code{tmat}}
  \item{draw}{logical vector}
}
\value{
  Returns a list with the following elements.  \code{cex} provides the
  character scaling used.  \code{xpos} provides the coordinates for the
  left edge of each column.  \code{ypos} provides the coordinates for
  the top of each row.

  Draws the contents of \code{tmat} as a side effect.  
}
\details{
  To be written
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
