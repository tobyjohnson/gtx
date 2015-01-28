\name{read.snpdata.plink}
\alias{read.snpdata.plink}
\title{Read genotype dosages in the format output by PLINK.}
\description{
  Reads genotype dosages as output by PLINK using the \code{--recodeA}
  option, and combines these with information about coded \emph{and
    noncoded} alleles from information output by PLINK using the
  \code{--freq} option.
}
\usage{
read.snpdata.plink(fileroot, tol.af = 0.01, phenotypes = NULL)
}
\arguments{
  \item{fileroot}{a filename root to which .frq and .raw extensions will be
    appended}
  \item{tol.af}{a tolerence for checking allele frequencies between the
    two files read}
  \item{phenotypes}{if not null, a data frame of phenotypes to be merged with the
    genotypes; must contain columns called FID and IID}
}
\details{
  Run PLINK twice, once with \code{--freq} and once with
  \code{--recodeA}, and with \emph{otherwise identical options} for e.g.
    individual and SNP inclusions and exclusions, in order to make
    suitable files for this function to read.
}
\value{
  Returns a list with snpinfo and data slots, see \code{\link{snpdata}}.
}
\examples{
\dontrun{
mypheno <- read.table("MYPHENO.dat", header = TRUE,
                      as.is = TRUE, na.strings = c("NA", "-9"))
system("plink --bfile MYBED --freq --out MYEXPORT")
system("plink --bfile MYBED --recodeA --out MYEXPORT")
mydata <- read.snpdata.plink("MYEXPORT", phenotypes = mypheno)
}
}
\references{
  Information about PLINK and its input and output data formats can be found at
  \url{http://pngu.mgh.harvard.edu/~purcell/plink}.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}