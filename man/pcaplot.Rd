\name{pcaplot}
\alias{pcaplot}
\title{Principal component plot coloured by factor}
\description{
  Plots pairwise plots for three principal components, with symbols
  coloured by a factor.
}
\usage{
pcaplot(pcs, f)
}
\arguments{
  \item{pcs}{A data frame or matrix with columns called PC1, PC2, and PC3}
  \item{f}{A factor with length matching the rows of pcs}
}
\details{
  The plot is drawn in 4 panels using \code{\link{split.screen}}, and
  therefore works inside \code{\link{gtxpipe}}'s plotting functions (but
  is incompatible with \code{par(mfrow)}, \code{\link{layout}} or other
  mechanisms).

  Two features of the plot are intended to give a more accurate visual
  impression of the density of different symbol types, in areas where
  they overlap.  First, symbols are plotted using dark but highly
  transparent colours.  Second, symbols that appear with higher
  frequencies in the plot are deliberately interleaved in the plotting
  order (so that all symbols of one type are not plotted \sQuote{on top}
  of another type), but symbols that appear with lower frequencies are
  placed towards the end of the plotting order (so tend to be plotted
  \sQuote{on top} for maximum visibility).

  The legend is plotted with font size scaling to force the labels to
  fit into the available space.
}
\value{
  Returns invisible null.
}
\examples{
pcaplot(pcs = cbind(PC1 = rnorm(500) + 7*rbinom(500, 1, .3), 
                    PC2 = rnorm(500),
                    PC3 = rnorm(500) + 6*rbinom(500, 1, .4)),
        f = c(rep("Class A", 230), rep("Class B", 200),
              rep("Class C", 25), rep("Class D", 20), 
              rep("Class E", 10), rep("Class F", 5), 
              rep("Class G", 4), rep("Class H", 3), 
              rep("Class J", 2), rep("Class K", 1))) 
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
