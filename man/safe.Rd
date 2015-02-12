\name{safe}
\alias{safe}
\title{Safely apply a function}
\description{
  Convenience function to safely apply a function
  when the argument may have length zero.
}
\usage{
safe(FUN, x, ...)
}
\arguments{
  \item{FUN}{a function to apply}
  \item{x}{first argument for \code{FUN}}
  \item{...}{additional arguments passed to \code{FUN}}
}
\value{
  The result of \code{FUN(x, ...)} after preprocessing \code{x}, or NA
  if \code{x} has length zero after preprocessing.
}
\examples{
\dontrun{
x <- 1:5
y <- rep("A", 5)
x[y == "B"]
min(x[y == "B"]) ## produces warning and returns Inf
safe(min, x[y == "B"]) ## returns NA
}
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
