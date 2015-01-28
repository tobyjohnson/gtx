\name{align.snpdata.coding}
\alias{align.snpdata.coding}
\title{Update genotype coding when there are coded allele designation flips.}
\description{
  The input parameterisation specify a desired coded allele for each
  SNP.  This function examines the coded and noncoded alleles used in
  the input genotype data, and for each SNP where the the input genotype
  data are encoded as the dose of the opposite (desired noncoded allele)
  allele, an additional column is added to the output genotype data with
  the dose of the desired coded allele.
}
\usage{
align.snpdata.coding(params, snpdata, ploidy = 2,
                     missing.snp = "fail")
}
\arguments{
  \item{params}{a data frame, see \link{gtx.params}.}
  \item{snpdata}{a list with snpinfo and data, see \link{snpdata}.}
  \item{ploidy}{if dosage for the noncoded allele is \var{x}, the dosage
    for the coded allele is calculated as \code{ploidy}-\var{x}.}
  \item{missing.snp}{character, either "fail" or "okay".}
}
\details{
  The PLINK convention of calling the coded allele \dQuote{0} for
  monomorphic SNPs is handled transparently, by assuming that the absent
  allele in the input genotype data matches whatever allele in the
  desired parameterisation does not match the present allele in the
  input genotype data.  This behaviour should not cause inadvertent
  strand flips.
  
  You should not need to call this function, unless you are intending to
  call \code{\link{grs.onesnp.apply}} without calling
  \code{\link{grs.make.scores}} first.  Note that
  \code{\link{grs.onesnp.apply}} has no way to check whether columns for
  desired coded alleles are present and may return NA for codes it
  cannot find.
  
  The \code{ploidy} argument defaults to 2, but should be set to 1 if
  the input genotype data are haplotypes (either phased or male X or Y
  chromosome).
  
  The \code{missing.snp} argument controls how to handle SNPs in the
  desired paramterisation that are not present in the input genotype
  data.  If "okay" then SNPs listed in the desired parameterisation but
  not present in the input genotype data are assumed to have dosage zero
  for all individuals.
  
  This function is one of the main computational bottlenecks and should
  be aggresively optimised in future releases.
}
\value{
  List with \code{$params} and \code{$snpdata} slots, contain the
  input arguments with additional columns.  The input \code{params} has
  an extra
  column \code{data.coded.freq} and the input \code{snpdata$data} has
  extra column(s) for doses of the specified coded alleles.
}
\examples{
data(mthfrex)
"rs1537514_G" \%in\% names(mthfrex$data) # FALSE
mthfrex <- align.snpdata.coding(mthfr.params, mthfrex)$snpdata  
"rs1537514_G" \%in\% names(mthfrex$data) # TRUE
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}