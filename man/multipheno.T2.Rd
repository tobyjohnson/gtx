\name{multipheno.T2}
\alias{multipheno.T2}
\title{Multi-phenotype test for association}
\description{
  For a set of phenotypes, given summary results for association tests
  for each single phenotype, over a large fixed set of marker genotypes
  (e.g. GWAS results), this function calculates a multi-phenotype
  association test for each marker that is equivalent to using the
  subject-specific data to perform a multivariate analysis of variance 
  for each marker.
}
\usage{
multipheno.T2(z, cor.use = NULL, cor.method = "spearman")
}
\arguments{
  \item{z}{a matrix of association Z statistics, with rows corresponding
    to markers and columns corresponding to phenotypes.}
  \item{cor.use}{a logical vector of length nrow(z), indicating a
    subset of markers to use to calculate the correlation matrix.}
  \item{cor.method}{a method to calculate the correlation matrix.}
}
\details{
  The function is named for the close correspondence with Hotelling's T2
  test.

  It is the user's responsibility to ensure that the columns of
  \dQuote{\code{z}} are marginally well calibrated, i.e. over a set of
  null markers, each column of \dQuote{\code{z}} should be standard
  normal marginal to the other columns of \dQuote{\code{z}}.

  Rank deficient situations are currently not handled.
}
\value{
  A list with three elements.  The first element,
  \dQuote{\code{rmatrix}}, is the estimated test statistic correlation
  matrix.  The second and third elements are both of length
  \code{nrow(z)}; 
  \dQuote{\code{T2}} contains the test statistics, which are chi-squared
  with \code{ncol(z)} d.f. under the null, and \dQuote{\code{pval}}
  contains the P-values.
}
\examples{
\dontrun{
nsubj <- 400 # number of subjects
nsnp <- 4000 # intended number of SNPs in GWAS

snps <- replicate(nsnp, rbinom(nsubj, 2, rbeta(1, 1.2, 1.2)))
## simulate with random allele frequencies
snps <- snps[ , apply(snps, 2, var) > 0]
nsnp <- ncol(snps) # number after filtering monomorphic

beta <- matrix(rnorm(30, 0, 0.1), ncol = 3)
## matrix of effects of 10 snps on 3 phenotypes

y1 <-  rnorm(nsubj)
y2 <- .1*y1 + rnorm(nsubj)
y3 <- .1*y1 + .3*y2 + rnorm(nsubj) # simulate correlated errors
y <- cbind(y1, y2, y3) + snps[ , 1:10] \%*\% beta
## wlog the truely associated snps are 1:10
rm(y1, y2, y3)

astats <- function(ii) {
  with(list(snp = snps[ , ii]),
       c(coef(summary(lm(y[ , 1] ~ snp)))["snp", 3],
         coef(summary(lm(y[ , 2] ~ snp)))["snp", 3],
         coef(summary(lm(y[ , 3] ~ snp)))["snp", 3],
         summary(manova(y ~ snp))$stats["snp", 6]))
}
system.time(gwas <- t(sapply(1:nsnp, astats)))
## columns 1-3 contain single phenotype Z statistics
## column 4 contains manova P-value
pv <- multipheno.T2(gwas[ , 1:3])$pval

plot(-log10(gwas[ , 4]), -log10(pv), type = "n",
     xlab = "MANOVA -log10(P)", ylab = "Summary statistic -log10(P)", las = 1)
points(-log10(gwas[-(1:10), 4]), -log10(pv[-(1:10)]))
points(-log10(gwas[1:10, 4]), -log10(pv[1:10]), cex = 2, pch = 21, bg = "red")
legend("topleft", pch = c(1, 21), pt.cex = c(1, 2), pt.bg = c("white", "red"),
       legend = c("null SNPs", "associated SNPs"))
## nb approximation gets better as nsnp becomes large, even with small nsubj
}
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
