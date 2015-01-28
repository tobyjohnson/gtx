\name{snphweCounts}
\alias{snphweCounts}
\title{Exact test of Hardy--Weinberg.}
\description{
  Computes a p-value for the exact test of Hardy--Weinberg proportions,
  for genotype counts for a biallelic locus, as described by Wiggington
  et al. (2005).
}
\usage{
snphweCounts(obs_hets, obs_hom1, obs_hom2)
}
\arguments{
  \item{obs_hets}{Observed count of the heterozygote genotype.}
  \item{obs_hom1}{Observed count of one homozygote genotype.}
  \item{obs_hom2}{Observed count of the other homozygote genotype.}
}
\details{
  This is the original code of Wiggington et al. (see references below),
  with minor modifications.  The original (and perhaps counterintuitive)
  order of arguments has been preserved. It is better to use the
  \code{\link{snphwe}} wrapper function.
} 
\value{
  The p-value.
}
\examples{
snphwe(rbinom(100, 2, 0.2))
}
\references{
  The method and algorithm used are described by Wiggington, Cutler and
  Abecasis (2005; A Note on Exact Tests of Hardy-Weinberg
  Equilibrium. American Journal of Human Genetics. 76(5):887-893)
  available at \url{http://www.ncbi.nih.gov/pmc/articles/PMC1199378}.
  The code of Wiggington et al. that is used here was downloaded from
  \url{http://www.sph.umich.edu/csg/abecasis/Exact}.  With appropriate
  citation, this code is freely available for use and can be
  incorporated into other programs.
}