\name{sanitise.whitespace}
\alias{sanitise.whitespace}
\title{Remove leading and trailing spaces and convert multiple to
  single spaces}
\description{
  Convenience function for a series of gsub() calls.
}
\usage{
sanitise.whitespace(tt)
}
\arguments{
  \item{tt}{a character vector.}
}
\value{
  A character vector with leading and trailing spaces removed, and
  double spaces converted to single spaces.
}
\examples{
sanitise.whitespace(c("  A to G", "A to  G", "A to G", "A    to G   "))
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
