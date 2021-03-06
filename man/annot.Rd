% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annot.R
\name{annot}
\alias{annot}
\title{Query analysis, variant, gene etc. annotation}
\usage{
annot(analysis, chrom, pos, ref, alt, rs, hgnc, ensg, max_variants = 30L,
  max_genes = 3L, dbc = getOption("gtx.dbConnection", NULL))
}
\arguments{
\item{analysis}{Analysis identifier}

\item{chrom}{Chromosome identifier (string)}

\item{pos}{Position}

\item{ref}{Reference allele}

\item{alt}{Alternate allele}

\item{rs}{dbSNP identifier (e.g. "rs123456")}

\item{hgnc}{HGNC gene symbol}

\item{ensg}{Ensembl gene identifier (e.g. "ENSG00000123456")}

\item{max_variants}{Maximum number of variants to return}

\item{max_genes}{Maximum number of genes to return}

\item{dbc}{Database connection (see \code{\link{gtxconnect}})}
}
\value{
Query results
}
\description{
Common queries of the annotation tables of the GWAS summary statistics.
}
\details{
All the main function arguments are optional, but at least one argument 
must be given in order to generate any results.
}
\author{
Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
