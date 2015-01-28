\name{combine.moments2}
\alias{combine.moments2}
\title{Builds a matrix of second moments for a set of individuals, from the
  matrices of second moments for all constituent subsets of individuals.}
\description{
  Given matrices of non-central second moments between a set of
  variables, calculated for two or more subsets of individuals, the
  matrix of non-central second moments between that set of variables, in
  the combined set of individuals, is obtained simply by addition.  This
  function provides a convenience interface that deals with different row and
  column orderings in the input matrices, and also allows a subset of
  variables to be treated as pertaining to specific subsets of
  individuals only.
}
\usage{
combine.moments2(xtxlist, fixed)
}
\arguments{
  \item{xtxlist}{a named list of matrices of second moments.}
  \item{fixed}{names of variables for which regression coefficients are
    considered fixed across all studies.}
}
\value{
  A numeric matrix of second moments.
}
\details{
  The names of the elements of \code{xtxlist} are used as identifiers
  for the subsets of individuals for which each matrix has been
  calculated.
  
  Variables whose names are included in the \code{fixed} argument must
  be included in \emph{all} matrices in \code{xtxlist}, and are assumed
  to have identical definitions in all individuals.  The intended
  application is that a regression model will be fitted with a
  coefficient for each of these variables that takes a single fixed
  value for all individuals.

  Any variables included in matrices in \code{xtxlist} whose names are
  not included in the \code{fixed} argument are assumed to have
  potentially different definitions in the subsets of individuals
  corresponding to the matrices in \code{xtxlist}.  An example would be
  ancestry principal component covariates, that might have the same
  names (\dQuote{PC1}, \dQuote{PC2}, ...) but do not have the same
  definitions in different subsets of individuals.  The intended
  application is that a regression model will be fitted with a
  subset-specific coefficient for each of these variables, and the
  variables are therefore renamed by prepending identifiers constucted
  using the names of the elements of \code{xtxlist}.

  The combined matrix has no \dQuote{ONE} intercept term, and instead
  has \code{length(xtxlist)} subset-specific intercepts (named by
  prepending identifiers constucted using the names of the elements of
  \code{xtxlist}).  This allows the the full rank of the combined matrix
  to be maintained while treating all subsets symetrically.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
\examples{
data(mthfrex)
## artifical example with two datasets obtained by
## splitting mthfrex by the HTN variable

xtx.hi <- make.moments2(mthfr.params, c("SBP", "SexC", "Age"),
                        as.snpdata(list(snpinfo = mthfrex$snpinfo,
                                   data = subset(mthfrex$data, HTN == 1))))

xtx.lo <- make.moments2(mthfr.params, c("SBP", "SexC", "Age"),
                        as.snpdata(list(snpinfo = mthfrex$snpinfo,
                                   data = subset(mthfrex$data, HTN == 0))))

## make list of X'X matrices
xtx.list <- list(hi = xtx.hi, lo = xtx.lo)

## combine for outcome SBP and fixed effects for all SNPs
## other variables SexC and Age will be treated as study-specific
fixed <- paste(mthfr.params$snp, mthfr.params$coded.allele, sep = "_")
xtx.comb <- combine.moments2(xtx.list, c("SBP", fixed))

## fit regression model
n.comb <- sum(diag(xtx.comb)[1:length(xtx.list)])
lm.moments2(xtx.comb, "SBP", c("hi_ONE", "lo_ONE", "rs4846052_T"), n = n.comb)

## equivalent regression directly using subject-specific data
coeff.extract(lm(SBP ~ HTN + rs4846052_C, data = mthfrex$data))
}