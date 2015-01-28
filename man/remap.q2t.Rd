\name{remap.q2t}
\alias{remap.q2t}
\title{Remap coordinates from BLAT query sequence to BLAT target sequence.}
\description{
  For a single BLAT match of a query sequence against a target sequence,
  a vector of nucleotide positions in the query sequence are remapped to
  the corresponding nucleotide positions in the target sequence.
}
\usage{
remap.q2t(qpos, blatres)
}
\arguments{
  \item{qpos}{a vector of 0-offset positions in the query sequence.}
  \item{blatres}{a single match from BLAT.}
}
\value{
  A vector of remapped positions in the target sequence.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
