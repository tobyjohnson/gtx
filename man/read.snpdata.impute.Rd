\name{read.snpdata.impute}
\alias{read.snpdata.impute}
\title{Read genotype dosages in the format output by IMPUTE.}
\description{
  Does what it says on the tin!
}
\usage{
  read.snpdata.impute(samplefile, genofile, phenotypes = NULL)
}
\arguments{
  \item{samplefile}{filename for samples, assumed 000Ps format.}
  \item{genofile}{filename for genotypes, assumed IMPUTE verbose format.}
  \item{phenotypes}{if not null, a data frame of phenotypes to be merged with the
    genotypes; first two columns will be used to match against first two
  columns of the samplefile.}
}
\details{
  This function needs better documentation.

  This function will be slow for large input files.  Best to use gtool
  or grep/awk out the relevant lines (SNPs) into a smaller file first.
}
\value{
  Returns a list with snpinfo and data slots, see \code{\link{snpdata}}.
}
