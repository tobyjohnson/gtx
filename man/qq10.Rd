\name{qq10}
\alias{qq10}
\alias{qq10.new}
\alias{qq10.envelope}
\alias{qq10.points}
\title{QQ plots on -log10 scale for p-values}
\description{
  Plots a QQ plot on a -log10 scale, appropriate for P-values.
  Component functions help build up more complex QQ plots.
}
\usage{
qq10(p, pmin = NULL, alpha = 0.01, ...)
qq10.new(pmin)
qq10.envelope(n, pmin, alpha = 0.01, col = "grey")
qq10.points(p, ...)
}
\arguments{
  \item{p}{A vector of observed p-values}
  \item{pmin}{The smallest observed or expected p-value to include in
    the plot}
  \item{alpha}{The alpha level for the envelope for the null distribution} 
  \item{n}{The number of p-values for the envelope for the null distribution}
  \item{col}{The colour for shading the envelope for the null distribution}
  \item{...}{Additional arguments passed to points for plotting the p-values}
}
\details{
  The plot is drawn in the current screen on the current device.

  If the default \code{pmin = NULL} is used, \code{qq10} sets
  \code{pmin} to the smallest observed or expected p-value.  If there
  are observed p-values much smaller than expected, this may not produce
  a nice plot, and the user should specify a value for \code{pmin} and
  note that the y axis will be truncated.
  
  The 100(1-alpha)\% envelope for the null distribution is approximate,
  and is based on asymptotics for the distribution of the expected
  P-value, pointwise for each value of the observed P-value.

  \code{qq10.points.fast} uses an experimental heuristic to thin points
  in the overlapping part of the QQ plot to speed up plotting.  Because
  this involves subsetting the points to be plotted, additional
  graphical parameters must be subsetted to match and hence must be
  passed explicitly.
}
\value{
  Nothing should be returned.
}
\examples{
## single command for simple QQ plot
qq10(sort(runif(100)), col = rainbow(100))
## example of building a plot to compare samples from two distributions
p1 <- runif(1000)
p2 <- rbeta(1000, shape1 = 0.6, shape2 = 1)
qq10.new(min(c(p1, p2, 1/2000)))
qq10.envelope(1000, min(c(p1, p2, 1/2000)))
qq10.points(p1, pch = 21, bg = "yellow")
qq10.points(p2, pch = 21, bg = "red")
legend("topleft", pch = rep(21, 2), pt.bg = c("red", "yellow"),
       legend = c(expression(beta(0.6,1)), expression(beta(1,1))))
box()
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
