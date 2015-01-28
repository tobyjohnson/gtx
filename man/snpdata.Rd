\name{snpdata}
\alias{is.snpdata}
\alias{as.snpdata}
\alias{summary.snpdata}
\alias{snpdata}
\title{Class for SNP genotype and phenotype data.}
\description{
  Objects of this class contain subject-specific SNP genotype or allele
  dosage data along with subject-specific phenotypes data, and SNP
  coded/noncoded allele information.
}
\usage{
is.snpdata(object)
as.snpdata(object)
\method{summary}{snpdata}(object, \dots)
}
\arguments{
  \item{object}{an object to test or coerce to snpdata class.}
  \item{\dots}{other arguments to \code{summary} are ignored.}
}
\details{
  An object of class \code{snpdata} is a list with at least two elements, called
  \dQuote{\code{snpinfo}} and \dQuote{\code{data}}.  \code{snpdata$snpinfo} must be a data frame
  with (at least) columns called \dQuote{\code{snp}}, \dQuote{\code{coded.allele}} and \dQuote{\code{noncoded.allele}}.
  \code{snpdata$data} must be a data frame with columns for each SNP in
  \code{snpdata$snpinfo} with column names formed by pasting together the
  \code{snp} and \code{coded.allele} with a separating underscore character,
  e.g. \dQuote{\code{rs12345_A}}.
}
\examples{
\dontrun{
align.snpdata.coding(params, snpdata)
grs.make.scores(params, snpdata)
grs.onesnp.apply(params, object)
make.moments2(params, phenolist, snpdata)
}
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
\seealso{
  Widely used formats for measured or
  imputed SNP genotype data are read by the \code{\link{read.snpdata.plink}},
  \code{\link{read.snpdata.mach}} and \code{\link{read.snpdata.impute}} functions,
  which all return data in the format described above.
}