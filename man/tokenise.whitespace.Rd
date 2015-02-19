\name{tokenise.whitespace}
\alias{tokenise.whitespace}
\title{Extract whitespace separated tokens}
\description{
  Converts whitespace separated tokens to a vector with one token per element.
}
\usage{
tokenise.whitespace(tt)
}
\arguments{
  \item{tt}{a character vector.}
}
\value{
  A character vector with one token per element.
}
\details{
  This is a simple convenience function for \code{unname},
  \code{strsplit}, \code{unlist} followed by removal of any zero length
  strings.

  This is used inside some other gtx functions where such as
  \code{\link{clinical.derive}}.
}
\examples{
tokenise.whitespace(c("pop demo", "pop exposure lab ae", "pop demo onctte"))
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
