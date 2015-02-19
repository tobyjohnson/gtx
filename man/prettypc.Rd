\name{prettypc}
\alias{prettypc}
\title{Pretty print percentages}
\description{
  Pretty print percentages.
}
\usage{
prettypc(x, digits = 0)
}
\arguments{
  \item{x}{a vector of integer counts}
  \item{digits}{precision for percentages}
}
\value{
  A character vector of percentages, rounded to \code{digits}
  precision.
}
\details{
  When \code{digits}>0, the pretty printing right pads with
  zeros.  Percentages that are 0 (or 100) to \code{digits} precision,
  but are not exactly 0 (or 100), are displayed like \dQuote{<1} (or
  \dQuote{>99}; as appropriate for \code{digits}).
}
\examples{
  counts <- c(2000, 150, 50, 10, 1, 0)
  prettypc(counts, digits = 0)
  prettypc(counts, digits = 1)
  prettypc(c(2000, 1))
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
