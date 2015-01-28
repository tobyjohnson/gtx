\name{prune.genes}
\alias{prune.genes}
\title{Prune markers by physical position}
\description{
  Selects subset of markers so that all excluded markers are within
  threshold physical genes of included markers.
}
\usage{
prune.genes(chr, pos, genetable, genes, win.size = 10000)
}
\arguments{
  \item{chr}{Chromosome}
  \item{pos}{Position}
  \item{genetable}{A data frame of a UCSC genome browser table}
  \item{genes}{Genes}
  \item{win.size}{Distance threshold}
}
\details{
  To write.

  \code{genetable} could be the refgene or GENCODE table from UCSC.  It
  must have columns \dQuote{chrom}, \dQuote{txStart}, \dQuote{txEnd},
  \dQuote{name2}.
}
\value{
  A vector of true/false.
}
\examples{
data(gencode14.UGT1A1)
markers <- data.frame(CHROM = paste("chr", c(rep(1, 1000), rep(2, 3000)), sep = ""),
                      POS = sort(runif(4000, 233000000, 235000000)),
		      pval = runif(4000, 1e-5, 1))
## within 5kb of either of two genes
markers$nearGene <- prune.genes(markers$CHROM, markers$POS, gencode14.UGT1A1, c("DNAJB3", "HJURP"), 5000)
with(markers, table(CHROM, nearGene))

leftrightpos <- c(234568893, 234781945)
plot.new()
plot.window(leftrightpos, c(-5, 5))
abline(h = 0, col = "grey")
gene.draw(2, leftrightpos[1], leftrightpos[2], gencode14.UGT1A1, genesep = 5000)
axis(1, at = pretty(leftrightpos * 1e-6)*1e6, labels = pretty(leftrightpos * 1e-6))
title(xlab = "chr2 genomic position (Mb)")
axis(2, at = 0:5, las = 1)
box()

with(markers, points(POS, -log10(pval), pch = 21, bg = ifelse(nearGene, "red", "grey")))
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
