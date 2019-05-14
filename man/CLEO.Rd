\name{CLEO}
\alias{CLEO}
\title{Conditional Leave-Each-Out (CLEO) analyses}
\description{
CLEO finemapping is based on an approach originally described by Giambartolomei et al. in their 2014 paper on colocalization.  Although Giambartolomei et al. only considered colocalization, CLEO generalizes the underlying idea.  CLEO provides a unifying framework for several problems concerning inference about causal variants and causal mechanisms.  One application of CLEO is similar in intent to other Bayesian finemapping approaches, notably FINEMAP.  Specifically, CLEO provides a Bayesian posterior probability, for any specific configuration of variants in a region, for that configuraiton to be causal for an observed association signal for a given disease or phenotype.  A different application of CLEO is to provide colocalization posterior probabilities when there are multiple causal variants in a region.

CLEO takes its name from "Conditional Leave Each Out", which is only one step (the middle step) of a three-step process to compute these posterior probabilities.

The advantages of CLEO over other finemapping approaches are:
\enumerate{
\item When there is more than one causal variant, arbitrary labels are assigned to the signals from each causal variant, which greatly assist in interpreting the results.
\item Besides computing posterior probabilities of causality for a given disease or phenotype, the approach extends to other finemapping-related questions such as colocalization between two traits or diseases (including gene expression traits) in the presence of multiple causal variants.
\item The computations are separated into three stages.  The first two are computationally intensive but do not depend on the Bayesian priors, which are needed for the final on-the-fly stage.  This means that results can quickly be (re-)generated under different priors e.g. the prior distribution on effect sizes, the prior probabilities on different variant types (coding, noncoding, etc) being causal, etc.
\item Compared to FINEMAP, the priors used are more flexible, are more intuitive to specify as probabilites, and lead to more intuitively acceptable results in terms of number of causal variants at a typical locus.
}

The three steps of a complete finemapping analysis using CLEO are: 
\enumerate{
\item Model selection.  Identify an upper bound on the number of causal variants or signals, and identify a "representative" index variant for each signal.
\item Conditional Leave-Each-Out (CLEO) analyses to derive "pure signal" beta and SE.
\item ABF calculation using the "pure signal" beta and SE.
}

Currently, step 1 is performed by running GCTA for each analysis/chromosome combination, in model selection mode, with a threshold P-value of 5e-06.  Step 2 is performed by running GCTA once for each index variant (for each analysis/chromosome), conditioning on all other index variants for that analysis/chromosome.  Separate scripts are used to coordinate and run these analyses at scale and to write the results to \code{TABLE gwas_results_joint} and \code{TABLE gwas_results_cond}.  Apart from checking whether results exist for a given analysis (using \code{\link{gtxanalyses}(..., has_cleo = TRUE)}), most users will not need to be concerned with the details of steps 1 and 2.

Step 3 is performed on the fly by several functions including \code{\link{regionplot}}, \code{\link{fm_cleo}}, \code{\link{coloc}} and \code{\link{multicoloc}}.

More details and empirical comparison with other methods and on simulated data to be added here...
}
\references{
Giambartolomei C, Vukcevic D, Schadt EE, Franke L, Hingorani AD, Wallace C, et al. (2014) Bayesian Test for Colocalisation between Pairs of Genetic Association Studies Using Summary Statistics. PLoS Genet 10(5): e1004383. \url{https://doi.org/10.1371/journal.pgen.1004383}
}
