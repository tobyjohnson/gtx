\name{manhattan}
\alias{manhattan}
\title{Manhattan plots for P-values}
\description{
  Plots a \sQuote{Manhattan} plot, appropriate for P-values from a
  genome-wide analysis.
}
\usage{
manhattan(p, SNP, chr, pos, ...)
}
\arguments{
  \item{p}{A vector of P-values}
  \item{SNP}{A vector of SNP names in chr:pos[:allele_allele] format}
  \item{chr}{A vector of chromosome locations}
  \item{pos}{A vector of positions within chromosomes}
  \item{...}{Additional arguments passed to points for plotting the p-values}
}
\details{
  The plot is drawn in the current screen on the current device.  The
  spacing between chromosomes is hard coded to a value that looks nice
  for human chromosome lengths.
}
\value{
  Nothing should be returned.
}
\examples{
chr <- rep(c(1:22, "X"), each = 5000)
pos <- do.call(c, lapply(c(2e8/sqrt(1:22), 2e8),
                         function(chrlen)
                         return(sort(chrlen*runif(5000)))))
SNP <- paste(chr, pos, sep = ":")
pvals <- runif(115000)
manhattan(pvals, SNP, pch = 20)
## FIXME blue colour repeated for X
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
