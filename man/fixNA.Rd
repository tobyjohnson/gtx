\name{fixNA}
\alias{fixNA}
\title{Fix NA Table Row and Column Names}
\description{
  Fixes NA table row and column names by converting to a string, so that
  tables produced with the useNA option can be printed using xtable.
}
\usage{
fixNA(tt, na.string = "NA")
}
\arguments{
  \item{tt}{a table}
  \item{na.string}{a string}
}
\details{
  Only intended for printing tables using xtable, so only works for 1-
  and 2-dimensional tables.
}
\value{
  The table with NA row and column names replaced with \code{na.string}.
}
\examples{
\dontrun{
x <- factor(c("A", "B", "C", NA)[rbinom(100, 3, .5) + 1])
y <- factor(c("foo", "bar", NA)[rbinom(100, 2, .5) + 1])
xtable(fixNA(table(x, useNA = "always")))
xtable(fixNA(table(y, x, useNA = "always")))
}
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
