\name{regionplot}
\alias{regionplot}
\alias{regionplot.region}
\alias{regionplot.new}
\title{Regional association plot}
\description{
  Draw regional association plot.
}
\usage{
regionplot(phenotype,
           chrom, pos_start, pos_end,
           hgncid, ensemblid, surround = 500000,
           style = 'signals', 
           protein_coding_only = TRUE,
           dbc = getOption("gtx.dbConnection", NULL))
regionplot.region(chrom, pos_start, pos_end,
                  hgncid, ensemblid, surround = 500000, 
                  dbc = getOption("gtx.dbConnection", NULL))
regionplot.new(chrom, pos_start, pos_end,
               hgncid, ensemblid, surround = 500000, 
               pmin = 1e-10,
	       protein_coding_only = TRUE,   
               dbc = getOption("gtx.dbConnection", NULL))
regionplot.genelayout(chrom, pos_start, pos_end,
                      ymax, cex = 0.75, 
		      protein_coding_only = TRUE, 
                      dbc = getOption("gtx.dbConnection", NULL))
regionplot.genedraw(gl)
regionplot.points(pos, pval,
                  pch = 21,
                  bg = rgb(.67, .67, .67, .5),
                  col = rgb(.33, .33, .33, .5),
                  cex = 1, 
                  suppressWarning = FALSE)
regionplot.recombination(chrom, pos_start, pos_end,
                         yoff = -.5, 
                         dbc = getOption("gtx.dbConnection", NULL))
}
\arguments{
  \item{phenotype}{The key value for the GWAS analysis to plot results from}
  \item{chrom}{Character specifying chromosome}
  \item{pos_start}{Start position of region}
  \item{pos_end}{End position of region}
  \item{hgncid}{HGNC identifier of gene to define region around}
  \item{ensemblid}{ENSEMBL gene identifier to define region around}
  \item{surround}{Distance around gene to include in region}
  \item{style}{Character specifying plot style}
  \item{protein_coding_only}{Whether to restrict annotation to protein
    coding genes only}
  \item{pmin}{Smallest P-value to size plot region for}
  \item{ymax}{Largest y-axis value to plan gene layout for}
  \item{cex}{Relative size for gene name labels}
  \item{gl}{Object returned by \code{regionplot.genelayout}}
  \item{yoff}{Y axis offset for baseline of recombination rate}
  \item{dbc}{Database connection}
}
\details{
  The highest level function, \code{regionplot}, draws a regional
  association plot, by querying the required data via a database
  connection (see \code{\link{gtxdbcheck}}).  Some lower level functions
  provide component parts to assist in generating modified plots
  (e.g. multiple analysis results in the same plot, or user specified
  plot colours or symbols).

  The region to plot results over can be specified in several different
  ways.  The region can be supplied as physical coordinates using the
  arguments \code{chrom}, \code{pos_start} and \code{pos_end}.
  Alternatively, the region can be centered on a gene of interest, using
  either the \code{hgncid} or \code{emsemblid} argument, and the size of
  region around the gene can be modified using the \code{surround}
  argument.
}
\value{
  \code{regionplot} returns an invisible null.  The plot is generated as
  a side effect.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
