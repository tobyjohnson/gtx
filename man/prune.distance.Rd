\name{prune.distance}
\alias{prune.distance}
\title{Prune markers by physical position}
\description{
  Selects subset of markers so that all excluded markers are within
  threshold physical distance of included markers.
}
\usage{
prune.distance(chr, pos, pval, win.size = 500000)
}
\arguments{
  \item{chr}{A vector of chromosome locations}
  \item{pos}{A vector of position (on chromosome) locations}
  \item{pval}{A vector of p-values to order by}
  \item{win.size}{Threshold distance}
}
\details{
  The pruning is done so that prune is TRUE only if there is no smaller
  P-value within \code{win.size} distance.
}
\value{
  A logical vector.
}
\examples{
tophits <- data.frame(chr = floor(runif(1000, 1, 22.99)),
                      pos = floor(runif(1000, 1, 1000000.99)),
                      pval = runif(1000))
tophits <- tophits[order(tophits$pval), ]
tophits <- within(tophits, prune <- prune.distance(chr, pos, pval))
head(tophits, 20)
with(tophits, table(chr[prune]))
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
