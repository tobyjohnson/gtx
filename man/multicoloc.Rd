\name{multicoloc}
\alias{multicoloc}
\title{Multiple colocalization analyses}
\description{
  Multiple colocalization analyses.
}
\usage{
multicoloc(analysis1, analysis2,
           chrom, pos_start, pos_end, pos, 
           hgncid, ensemblid, rs, surround = 0,
           hard_clip = FALSE,
           style = 'heatplot',
           thresh_analysis = 0.1, thresh_entity = 0.1, 
           dbc = getOption("gtx.dbConnection", NULL))
multicoloc.data(analysis1, analysis2,
                chrom, pos_start, pos_end, pos, 
                hgncid, ensemblid, rs, surround = 0,
                hard_clip = FALSE, 
                dbc = getOption("gtx.dbConnection", NULL))
}
\arguments{
  \item{analysis1}{The key value(s) for GWAS analysis/es to analyze}
  \item{analysis2}{The key value for the second GWAS analysis to analyze}
  \item{chrom}{Argument passed to \code{\link{gtxregion}()}}
  \item{pos_start}{Argument passed to \code{\link{gtxregion}()}}
  \item{pos_end}{Argument passed to \code{\link{gtxregion}()}}
  \item{pos}{Argument passed to \code{\link{gtxregion}()}}
  \item{hgncid}{Argument passed to \code{\link{gtxregion}()}}
  \item{ensemblid}{Argument passed to \code{\link{gtxregion}()}}
  \item{surround}{Argument passed to \code{\link{gtxregion}()}}
  \item{hard_clip}{Logical, see details}
  \item{style}{Character specifying plot style(s)}
  \item{thresh_analysis}{Probability threshold for inclusion in plots}
  \item{thresh_entity}{Probability threshold for inclusion in plots}
  \item{dbc}{Database connection}
}
\details{
  \code{multicoloc()} is an entry point for multiple colocalization
  analyses.  It supports the most common use case, to colocalize
  association signals from one or more analyses of gene
  expression/protein levels (specified by \code{analysis1}), each of
  which includes association statistics for multiple entities (genes or
  proteins), against an association signal from a single analysis
  (typically a disease or clinical phenotype, specified by
  \code{analysis2}).  For this use case, \code{multicoloc()} is
  typically more convenient and (much) faster than looping over multiple
  calls to \code{\link{coloc}()}.

  \code{multicoloc()} offers a choice of two different algorithms for
  controlling the genomic region from which summary statistics are used
  for colocalization analyses, controlled by the argument
  \code{hard_clip}.  The default, \code{hard_clip=FALSE}, uses the full
  set of available summary statistics for the entity/ies analyzed from
  each analysis included in \code{analysis1}.  In this mode, the genomic
  range arguments \code{chrom}, \code{pos_start}, \code{pos_end},
  \code{hgncid} etc. are only used (via \code{\link{gtxregion()}}) to
  determine the set of entities to be analyzed.  Typically, this results
  in different entities being analyzed for colocalization using
  different (albeit overlapping) regions of the association signal from
  \code{analysis2}.  The alternative \code{hard_clip=TRUE}, uses only
  summary statistics within the genomic range specified (via the
  arguments passed to \code{\link{gtxregion()}}).  Typically, different
  entities will be analyzed for colocalization using the same or similar
  regions of the association signal from \code{analysis2}, depending on
  how the genomic range overlaps the summary statistics available for
  each entity.  The exact algorithms used in each mode are detailed
  below.  (And can be visualized using plot style= in a forthcoming
  release.)

  When \code{hard_clip=FALSE}, the algorithm used by \code{multicoloc()}
  first determines a \dQuote{seed region} using the genomic region
  arguments, as interpreted by \code{\link{gtxregion()}}.  Next, a set
  of entities is determined from the summary statistics for all analyses
  included in \code{analysis1}, consisting of all entites with summary
  statistics overlapping this \dQuote{seed region}.  (Better
  implementation of overlap is forthcoming).  Finally, an
  \dQuote{expanded region} is determined, that includes all available
  summary statistics for all of these entities.  This \dQuote{expanded
  region} is then used for each colocalization analyses, for each entity
  within each analysis within \code{analysis1}, against
  \code{analysis2}.  Notes and Warnings: This algorithm only makes sense
  if the summary statistics are restricted to localized regions around
  each entity, such as cis- regions for eQTL analyses.  Typically,
  different entities will be evaluated for colocalization using
  different regions of summary statistics for \code{analysis2}.  Because
  the set of entities is determined by aggregating over all analyses in
  \code{analysis1}, unexpected results may be produced if a given entity
  has summary statistics at very different genomic positions in
  different analyses.  The set of entities is combined across
  \dQuote{Seed regions} specified using only the index variant from a
  GWAS signal (e.g. using \code{pos} or \code{rs} with the default
  \code{surround=0}) will not guarantee to select all entities with
  summary statistics for cis- regions spanning such a single base pair
  \dQuote{seed} region, if some entities are missing summary statistics
  for the variant in that \dQuote{seed} region.  [This last issue will
  be fixed in a forthcoming update.]

  When \code{hard_clip=FALSE}, the algorithm used by \code{multicoloc()}
  is simply to select all summary statistics within the genomic region
  arguments, as interpreted by \code{\link{gtxregion()}}.  The typical
  use case is to set this genomic region as the extent of the
  \sQuote{significant} part of the association signal for
  \code{analysis2}.  The \code{hard_clip=FALSE} mode is (currently) not
  the default option, because in initial exploratory analyses it is
  unusual to precisely specify this region, and because we believe the
  number of \sQuote{false positive} colocalizations is reduced by
  including the whole cis- eQTL region (assuming that the strongest
  disease signal in the region \sQuote{should} be aligned with the
  strongest cis- eQTL signal).  Notes and Warnings: In general, a given
  entity may have summary statistics that only partially overlap the
  genomic region specified, which may have unexpected consequences.  In
  a future release it will be possible to automatically subset to
  entities that overlap the genomic region specified by more than a
  chosen percentage.  When using a \code{hgncid} or \code{ensemblid}
  gene identifier to specify the region from which to use summary
  statistics, the default \code{surround=0} will \emph{not} include the
  full cis eQTL region.

  In a future release the output of multicoloc will be a long skinny
  dataframe with the full colocalization results (all priors, bfs and
  posteriors, numbers of variants and min and max positions used).
}
\value{
  \code{multicoloc} returns a data frame containing the result of the
  colocalization analyses, see \code{\link{coloc.fast}} for details.
  The plot is generated as a side effect.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
		  
