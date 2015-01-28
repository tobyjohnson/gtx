\name{gene}
\alias{gene.annotate}
\alias{gene.nearest}
\alias{gene.draw}
\title{Annotate genomic positions with gene transcript information}
\description{
  Annotates genomic positions (chromosome and position) with information
  about nearby transcripts from a UCSC genome browser table such as GENCODE.
}
\usage{
gene.annotate(chr, pos, genetable, win.size = 10000)
gene.nearest(chr, pos, genetable)
gene.draw(chr, leftpos, rightpos,
          genetable, nodraw = NULL, 
          genesep = 10000, hlines.min = NULL, yhi = -1, ylo = -5,
          exony = 0.05, genecex = 1)
}
\arguments{
  \item{chr}{A vector of chromosome numbers}
  \item{pos}{A vector of chromosomal base pair positions (one-based)}
  \item{genetable}{A data frame of a UCSC genome browser table}
  \item{win.size}{Size of window to include transcripts for}
  \item{leftpos}{Leftmost position to draw transcripts for}
  \item{rightpos}{Rightmost position to draw transcripts for}
  \item{nodraw}{A list of gene names not to draw}
  \item{genesep}{Desired space between genes drawn on same horizontal line}
  \item{hlines.min}{Minimum number of horizontal lines to use}
  \item{yhi}{Y coordinate of highest line}
  \item{ylo}{Y coordinate of lowest line}
  \item{exony}{Y coordinate height for exons}
  \item{genecex}{cex parameter for gene names}
}
\details{
  \code{genetable} could be the refgene or GENCODE table from UCSC.  It
  must have columns \dQuote{chrom}, \dQuote{txStart}, \dQuote{txEnd},
  \dQuote{name2}, and for \code{gene.draw} it must also have columns
  \dQuote{exonStarts} and \dQuote{exonEnds}.  Note that all positions in
  \code{genetable} are assumed to be zero-based, but function argument
  \code{pos} is assumed to be one-based.
  
  The algorithm for \code{gene.nearest} is peculiar and may change in
  future versions.  If the queried position is inside one or more
  transcripts, their names are concatenated (comma separated).
  Otherwise, the names of the nearest transcripts to the left and right
  are concatendated (hyphen separated).  This may produce curious output
  for transcript names that themselves contain a hyphen.
    
  \code{gene.draw} draws representations of all genes in the specified
  region.  For each gene (determined by unique values of
  \dQuote{name2}), all transcripts are drawn in an overlapping style,
  with yellow blocks for exons and black lines for introns, and
  \dQuote{name2} written below.  Different genes are distributed over a number of lines,
  chosen adaptively to avoid horizontal overlap.  In large or gene dense
  regions, this may result in overlap of the gene names with the  vertical overlap. 
}
\value{
  Character string, or plot.
}
\examples{
\dontrun{
## genome wide annotation 
genetable <- read.table(gzfile("path/to/UCSC-genome-tables/hg19.GENCODE14.gz"),
                          header = TRUE, comment.char = "", sep = "\t", as.is = TRUE)
}

data(gencode14.UGT1A1)
leftrightpos <- c(234568893, 234781945)
plot.new()
plot.window(leftrightpos, c(-5, 5))
abline(h = 0, col = "grey")
gene.draw(2, leftrightpos[1], leftrightpos[2], gencode14.UGT1A1, genesep = 5000)
axis(1, at = pretty(leftrightpos * 1e-6)*1e6, labels = pretty(leftrightpos * 1e-6))
title(xlab = "chr2 genomic position (Mb)")
axis(2, at = 0:5, las = 1)
box()

querypos <- c(234680000, 234685000)
points(querypos, rep(0, 2), pch = 23)
gene.nearest(c("chr2", "chr2"), querypos, gencode14.UGT1A1)

gene.annotate("chr2", querypos[1], gencode14.UGT1A1) # genes within 10kb
gene.annotate("chr2", querypos[1], gencode14.UGT1A1, win.size = 0) # genes within 0kb
gene.annotate(2, querypos[1], gencode14.UGT1A1, win.size = 0) # same
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
