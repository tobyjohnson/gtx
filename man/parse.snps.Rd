\name{parse.snps}
\alias{parse.snps}
\title{Parse text representation of a SNP embedded in flanking sequences.}
\description{
  Extracts the first occurence of a pattern like [A/T], along with the
  5' and 3' flanking sequences.
}
\usage{
parse.snps(seqs, seqnames = NULL)
}
\arguments{
  \item{seqs}{a (vector of) sequences.}
  \item{seqnames}{optional names.}
}
\value{
  A data frame with seven columns, called seq5p, len5p, poly, alleleA,
  alleleB, seq3p, and len3p.
}
\examples{
data(snps.BRCA1)
parse.snps(snps.BRCA1$SourceSeq[1:10])
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
