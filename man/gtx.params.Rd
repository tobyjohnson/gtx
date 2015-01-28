\name{gtx.params}
\alias{gtx.params}
\title{Parameter format for multi-SNP analyses}
\description{
  Format used by functions in the gtx package for parameterisations for
  multi-SNP analyses.  The parameterisation describes the SNPs to be
  used, choice of coded and noncoded alleles, and potentially allele
  frequency and effect size information.
}
\details{
  \code{params} must be a data frame with (at least) columns called \dQuote{\code{snp}}, \dQuote{\code{coded.allele}} and
  \dQuote{\code{noncoded.allele}}.  For some functions, additional columns called
  \dQuote{\code{coded.freq}} and \dQuote{\code{coef}} are also required.
  
  Suitable parameterisations in the required format are provided for
  some genetic risk score analyses by \code{\link{cad.scores}}, \code{\link{t2d.scores}} etc.
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
