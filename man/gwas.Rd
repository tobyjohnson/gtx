\name{gwas}
\alias{gwas}
\title{GWAS summary}
\description{
  Summarize and visualise GWAS results.
}
\usage{
gwas(analysis,
     style = 'manhattan',
     pval_thresh = 5e-08, maf_ge, rsq_ge,
     gene_annotate = TRUE,
     manhattan_thresh = 5e-08, 
     manhattan_ymax = 30, 
     manhattan_col = c('#064F7C', '#6D97BD'),
     manhattan_interspace = 50e6,
     manhattan_fastbox = 2, 
     dbc = getOption("gtx.dbConnection", NULL))
}
\arguments{
  \item{analysis}{The key value for the GWAS analysis to summarize}
  \item{style}{Character specifying plot style(s)}
  \item{pval_thresh}{P-value threshold for top hits table}
  \item{maf_ge}{MAF greater-or-equal threshold}
  \item{rsq_ge}{Imputation r-squared greater-or-equal threshold}
  \item{gene_annotate}{Logical, whether to annotate top hits with gene names}
  \item{manhattan_thresh}{P-value for threshold line on Manhattan plot}
  \item{manhattan_ymax}{Y axis maximum, -log10(p) scale, for Manhattan plot}
  \item{manhattan_col}{Chromosome colour cycle for Manhattan plot}
  \item{manhattan_interspace}{Chromosome interspacing for Manhattan plot}
  \item{manhattan_fastbox}{Box height for fast Manhattan plot}
  \item{dbc}{Database connection}
}
\details{
  This high level function summarizes and visualizes results
  (i.e. complete summary statistics) from a single GWAS analysis.

  Currently, a distance-pruned table of top hits, and an optional
  Manhattan plot, are generated.

  The distance pruning and gene annotation is currently very slow and
  needs to be optimized.

  The Manhattan plot uses an approximation to speed up both the data
  query and the plot rendering, where boxes are drawn instead of
  individual points, for y-axis values between 0 and
  \code{manhattan_fastbox}.

  Y-axis points are truncated, but this needs to be made clear in the
  plot.
}
\value{
  \code{gwas} returns a data frame containing the distance pruned top
  hits.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
		  