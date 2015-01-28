\name{contrasting.rainbow}
\alias{contrasting.rainbow}
\title{Rainbow of colours permuted to maximise contrast.}
\description{
  Returns a rainbow of colours, permuted so that adjacent elements
  constrast each other as much as possible.
}
\usage{
contrasting.rainbow(x, ...)
}
\arguments{
  \item{x}{number of colours in the rainbow.}
  \item{...}{other arguments to be passed to \code{rainbow}.}
}
\details{
  The permutation is chosen by finding the largest coprime of \code{x}
  that is not greater than \code{x}/2, and using that coprime as an
  increment size to cycle through the normal rainbow.
}
\value{
  A vector of colours of length \code{x}.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
\examples{
## contrasting colours suitable for Manhattan plot for 22 autosomes
plot(1:22, rep(1, 22), pch = 22, cex = 2, ann = FALSE, 
     yaxt = "n", bg = contrasting.rainbow(22))
}