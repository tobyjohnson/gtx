\name{blockassoc}
\alias{blockassoc}
\alias{blockstats}
\alias{blockstats.lm}
\alias{blockstats.glm}
\alias{blockstats.coxph}
\alias{blockstats.clm}
\title{Association tests for blocks of genotypes (one by one)}
\description{
  Calculations association test statistics for each genetic variant in a
  block.
}
\usage{
blockassoc(qcall, data, minimac,
usubjid = getOption("gtx.usubjid", "USUBJID"), 
threshold.MAF = 0, threshold.Rsq = 0, threshold.pass = NULL,
message.begin = "blockassoc", out.signif = 6, use.compiled = FALSE)
\method{blockstats}{lm}(m1, m0, coefname = "GENOTYPE")
\method{blockstats}{glm}(m1, m0, coefname = "GENOTYPE")
\method{blockstats}{coxph}(m1, m0, coefname = "GENOTYPE")
\method{blockstats}{clm}(m1, m0, coefname = "GENOTYPE")
}
\arguments{
  \item{qcall}{A quoted call to a model fitting function}
  \item{data}{Data frame (phenotype data) in which qcall is to be evaluated}
  \item{minimac}{Path to minimac inputed genotype info and dosages}
  \item{usubjid}{Name of the unique subject id variable in \code{data}}
  \item{threshold.MAF}{Minimum minor allele frequency (MAF) for inclusion in output}
  \item{threshold.Rsq}{Minimum imputation quality (Rsq) for inclusion in output}
  \item{threshold.pass}{Character vector of SNP names to include in
    output independent of MAF or Rsq}
  \item{message.begin}{Optional text to begin each progress message with}
  \item{out.signif}{Number of significant digits for output}
  \item{use.compiled}{Flag whether to use fast C++ code for special cases}
  \item{m1}{Fitted model object with genotype term}
  \item{m0}{Fitted model object without genotype term}
  \item{coefname}{Name of the coefficient for genotype in \code{m1}}
}
\details{
  TO BE WRITTEN

  TO DO - Add arguments for genotypes to be an existing data frame,
  PLINK binary or --recodeA output, or VCF
}
\value{Returns a data frame with info and association statistics.}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
\examples{
\dontrun{
library(ordinal)
library(survival)
library(MASS)

write.table(blockassoc(quote(lm(simCount ~ sexMale + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10)),
                       data = read.table("phenotypes/simCount.pheno", header = TRUE, as.is = TRUE),
                       minimac = "imputed-20120314/chr22chunk6",
                       threshold.MAF = 0.01, threshold.Rsq = 0.5),
            file = "test-lm.out", row.names = FALSE, quote = FALSE)

write.table(blockassoc(quote(glm(simCount ~ sexMale + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
family = poisson)),
                       data = read.table("phenotypes/simCount.pheno", header = TRUE, as.is = TRUE),
                       minimac = "imputed-20120314/chr22chunk6",
                       threshold.MAF = 0.01, threshold.Rsq = 0.5),
            file = "test-glm.out", row.names = FALSE, quote = FALSE)

write.table(blockassoc(quote(clm(factor(simCount) ~ sexMale + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10)),
                       data = read.table("phenotypes/simCount.pheno", header = TRUE, as.is = TRUE),
                       minimac = "imputed-20120314/chr22chunk6",
                       threshold.MAF = 0.01, threshold.Rsq = 0.5),
            file = "test-clm.out", row.names = FALSE, quote = FALSE)

write.table(blockassoc(quote(coxph(Surv(simLife, simSurv) ~ sexMale + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10)),
                       data = read.table("phenotypes/simLife.pheno", header = TRUE, as.is = TRUE),
                       minimac = "imputed-20120314/chr22chunk6",
                       threshold.MAF = 0.01, threshold.Rsq = 0.5),
            file = "test-coxph.out", row.names = FALSE, quote = FALSE)
}
}