\name{snphwe}
\alias{snphwe}
\title{Exact test of Hardy--Weinberg.}
\description{
  Computes a p-value for the exact test of Hardy--Weinberg proportions,
  for genotype counts for a biallelic locus, as described by Wiggington
  et al. (2005)
}
\usage{
snphwe(g)
}
\arguments{
  \item{g}{Vector of genotypes.}
}
\details{
  This function is a wrapper for \code{\link{snphweCounts}}, which uses
  code by Wiggington et al. to compute a P-value for an exact test of
  Hardy--Weinberg proportions.
  
  This function requires that when \code{g} is sorted, the order is
  hom/het/hom.  For example \code{g} can be additively coded genotypes
  (0, 1, 2), or text coded genotypes like AA, AC, CC.  For text coded
  genotypes the heterozygote MUST ALWAYS be coded the same way, e.g. a
  mixture of AC and CA are not allowed.
} 
\value{
  The p-value.
}
\examples{
snphwe(rbinom(100, 2, 0.2))
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
