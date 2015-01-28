\name{hapmap.read.haplotypes}
\alias{hapmap.read.haplotypes}
\title{Read hapmap haplotypes.}
\description{
  Reads HapMap haplotypes and legend and stores in format expected by functions in gtx package.
}
\usage{
hapmap.read.haplotypes(chr, path, sample = "CEU", release = 22)
}
\arguments{
  \item{chr}{chromosome to read}
  \item{path}{path to local mirror of ftp://ftp.ncbi.nlm.nih.gov/hapmap}
  \item{sample}{a HapMap sample, e.g. CEU or YRI}
  \item{release}{a HapMap release version, currently must be 21 or 22}
}
\details{
  You need a local mirror of HapMap data.  To build this use e.g. a
  shell command\cr
  \code{rsync ftp://ftp.ncbi.nlm.nih.gov/hapmap ~}\cr
  then use an argument \code{path="~/hapmap"} for this function call.

  Technically the local mirror only needs the contents of\cr
  \url{ftp://ftp.ncbi.nlm.nih.gov/hapmap/phasing/2006-07_phaseII/phased/}\cr
  and\cr
  \url{ftp://ftp.ncbi.nlm.nih.gov/hapmap/phasing/2007-08_rel22/phased/}

  Other gtx functions expect to be passed a list with 22 elements; one
  for each chromosome.  See the examples.
}
\examples{
\dontrun{
## chromosomes 1-22
hapmap <- lapply(1:22, hapmap.read.haplotypes, path = "~/hapmap")

## selected chromosomes only
hapmap <- lapply(1:22, function() return(NULL))
for (chr in c(1, 11, 17)) {
  hapmap[[chr]] <- hapmap.read.haplotypes(chr, path = "~/hapmap")
}
}
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}