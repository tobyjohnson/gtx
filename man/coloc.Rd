% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coloc.R
\name{coloc}
\alias{coloc}
\title{Colocalization analysis.}
\usage{
coloc(analysis1, analysis2, signal1, signal2, chrom, pos_start, pos_end,
  pos, hgncid, ensemblid, rs, surround = 5e+05, entity, entity1, entity2,
  style = "Z", priorsd1 = 1, priorsd2 = 1, priorc1 = 1e-04,
  priorc2 = 1e-04, priorc12 = 1e-05, join_type = "inner",
  dbc = getOption("gtx.dbConnection", NULL))
}
\arguments{
\item{analysis1}{The key value for the first GWAS analysis to analyze}

\item{analysis2}{The key value for the second GWAS analysis to analyze}

\item{chrom}{Character specifying chromosome}

\item{pos_start}{Start position of region}

\item{pos_end}{End position of region}

\item{hgncid}{HGNC identifier of gene to define region around}

\item{ensemblid}{ENSEMBL gene identifier to define region around}

\item{surround}{Distance around gene to include in region. Default: 500000}

\item{entity}{Identifier for an entity, for analyses of multiple entities}

\item{entity1}{Identifier for an entity, for analysis}

\item{entity2}{Identifier for an entity, for analysis2}

\item{style}{Character specifying plot style. Default: 'Z'}

\item{dbc}{Database connection. Default: getOption("gtx.dbConnection", NULL)}
}
\value{
a data frame containing the result of the
 colocalization analysis, see \code{\link{coloc.compute}} for details.
 The plot is generated as a side effect.
}
\description{
Colocalization analysis using summary statistics from database.
}
\details{
This high level function conducts a colocalization analysis, using
summary statistics for association with two traits, across a region of
the genome.  The \code{analysis1} and \code{analysis2} arguments specify
which analyses to use summary statistics from.
Where one or both analyses 
have summary statistics for multiple entities (e.g. from eQTL or
pQTL analyses), the desired entities must be specified (see below).

The \code{style} argument can be set to \code{'Z'} to plot Z
statistics for the two analyses, and/or \code{'beta'} to plot beta (effect
size) statistics for the two analyses.  \dQuote{One-sided} plots, where
the ref/alt alleles are flipped so that beta is always positive for
\code{analysis1}, are provided as styles \code{'Z1'} and \code{'beta1'}.
The style \code{'none'} suppresses plotting altogether.  In all plots,
the x and y axes are used for analysis1 and analysis2 respectively.

To help interpretation, 95% credible sets are calculated
separately for each analysis.  Variants in both credible sets are plotted
as diamonds, and variants only in one credible set are plotted as triangles
(down for the x axis analysis1; up for the y axis analysis 2).  Additionally,
two channel shading is used to indicate posterior probability of
causality for analysis1 (red), analysis2 (green), or both (yellow).

Note that when using a \code{hgncid} or \code{ensemblid} gene
identifier to specify the region from which to use summary statistics,
the default \code{surround=500000} will \emph{not} include the full
cis eQTL region as usually specified.
 
The region of interest can be specified in several different ways.
The region can be supplied as physical coordinates using the arguments
\code{chrom}, \code{pos_start} and \code{pos_end}.  Alternatively, the
region can be centered on a gene of interest, using either the
\code{hgncid} or \code{emsemblid} argument, and the size of region
around the gene can be modified using the \code{surround} argument.
Note that the primary purpose of gene-identifying arguments
\code{hgncid} or \code{ensemblid} is to specify the genomic region of
interest (and thus the set of the variants to analyse).  It is only a
secondary purpose that the entity for eQTL or pQTL analyses will be
inferred from \code{hgncid} or \code{ensemblid}, if no explicit
\code{entity} argument is given.

Entities are used to distinguish genomic features, where a single set
analysis includes summary statistics, for each variant, for
associations with one or more entities.  E.g. in an eQTL analysis,
each transcript or gene is an entity, and a single typical variant
will have summary statistics for associations with multiple
transcripts or genes.  If either of the analyses specified by
\code{analysis1} and \code{analysis2} have results separated by
entity, then the arguments \code{entity1} and \code{entity2} are used
to specify the desired entity from each.  If either \code{entity1} or
\code{entity2} is missing, the argument \code{entity} is used
instead.  (This mechanism facilitates e.g. colocalization between analyses
for the same transcript between two different eQTL datasets.)  If the
argument \code{entity} is also missing, the function attempts to infer
a suitable entity from the \code{hgncid} or \code{ensemblid}
arguments.  (This leads to sensible default behaviour, and facilitates
the most common use case of centering the genomic region of interest
on the entity being analysed in an eQTL or pQTL dataset.)
}
\author{
Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
