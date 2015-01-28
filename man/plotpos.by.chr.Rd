\name{plotpos.by.chr}
\alias{plotpos.by.chr}
\alias{plotcol.by.chr}
\alias{axis.by.chr}
\title{Utilities for Manhatten plots}
\description{
  Constructs vectors of plotting positions and colours, for points that
  have positions on chromosomes.
}
\usage{
plotpos.by.chr(chr, pos, gap = 5e7,
               chrset = c(1:22, "XY", "X", "Y", "M"))
plotcol.by.chr(chr, cols = NULL, col.missing = "black",
               chrset = c(1:22, "XY", "X", "Y", "M"))
axis.by.chr(chr, plotpos, side = 1, lines = 2)
}
\arguments{
  \item{chr}{A vector of chromosome identities}
  \item{pos}{A vector of positions within chromosomes}
  \item{gap}{Gap size between chromosomes, in same units as \code{pos}}
  \item{chrset}{Set of all chromosomes}
  \item{cols}{Colours for chromosomes, recycled as necessary}
  \item{col.missing}{Colour to use for chromosomes not in \code{chrset}}
  \item{plotpos}{A vector of plotting positions}
  \item{side}{Side to plot the axis}
  \item{lines}{Number of lines to distribute chromosome labels over}
}
\details{
  These functions all convert \code{chr} to upper case and remove leading
  \dQuote{CHR} if present, before doing anything else.

  If argument \code{cols} to \code{plotcol.by.chr} is left NULL, it is
  filled with a suitable call to \code{\link{contrasting.rainbow}}.
  For an alternative popular style of Manhatten plot, try \code{cols = c("cyan", "grey")}.
 
  The example illustrates how to use these functions to help draw a Manhatten plot.
}
\value{
  \code{plotpos.by.chr} returns a vector of plotting positions of length
  equal to the input \code{chr} and \code{pos}.

  \code{plotcol.by.chr} returns a vector of plotting colours of length
  equal to the input \code{chr}.

  \code{axis.by.chr} draws an axis and returns the positions the ticks
  were drawn at.
}
\examples{
## generate dummy GWAS results
gwas <- data.frame(chr = paste("chr", rep(1:22, each = 500), sep = ""),
                   pos = do.call(c, lapply(1:22, function(chr)
return(sort(ceiling(runif(500, 1, 1e8/chr)))))),
                   pval = runif(11000),
                   stringsAsFactors = FALSE)
## example drawing a Manhatten plot
gwas$plotpos <- plotpos.by.chr(gwas$chr, gwas$pos)
gwas$plotcol <- plotcol.by.chr(gwas$chr)
with(gwas, plot(plotpos, -log10(pval), col = plotcol,
                ann = FALSE, xaxt = "n", las = 1))
title(ylab = expression(-log[10](italic(P))))
axis.by.chr(gwas$chr, gwas$plotpos)
mtext("Genomic position by chromosome", side = 1, line = 3)
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
