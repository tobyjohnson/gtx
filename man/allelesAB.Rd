\name{allelesAB}
\alias{allelesAB}
\title{Paste together (vectors of) A and B alleles after sorting alphabetically.}
\description{
  A tool to summarise and compare alleles when we do not care what order
  they are reported in, i.e. when A/G and G/A are considered the same.
}
\usage{
allelesAB(A1, A2, sep = "/")
}
\arguments{
  \item{A1}{a vector of allele names.}
  \item{A2}{a vector of allele names, of same length as \code{A1}.}
  \item{sep}{string to use as a separator.}
}
\value{
  A vector of combined allele names; elementwise \code{A1} and \code{A2}
  are pasted together in alphabetically increasing order.
}
\examples{
  data(t2d.scores)
  table(allelesAB(t2d.scores$coded.allele, t2d.scores$noncoded.allele))
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}