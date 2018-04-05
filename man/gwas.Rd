\name{gwas}
\alias{gwas}
\title{GWAS summary}
\description{
  Summarize and visualise GWAS results.
}
\usage{
gwas(analysis,
     style = c('manhattan', 'qqplot'), 
     pval_thresh = 5e-08, maf_ge, rsq_ge,
     gene_annotate = TRUE,
     plot_ymax = 30,
     manhattan_thresh = 5e-08,
     manhattan_col = c('#064F7C', '#6D97BD'),
     manhattan_interspace = 50e6,
     qqplot_col = '#064F7C',
     qqplot_alpha = 0.01,
     plot_fastbox = 2, 
     dbc = getOption("gtx.dbConnection", NULL))
}
\arguments{
  \item{analysis}{The key value for the GWAS analysis to summarize}
  \item{style}{Character vector specifying plot style(s)}
  \item{pval_thresh}{P-value threshold for top hits table}
  \item{maf_ge}{MAF greater-or-equal threshold}
  \item{rsq_ge}{Imputation r-squared greater-or-equal threshold}
  \item{gene_annotate}{Logical, whether to annotate top hits with gene names}
  \item{plot_ymax}{Y axis maximum, -log10(p) scale, for plots}
  \item{manhattan_thresh}{P-value for threshold line on Manhattan plot}
  \item{manhattan_col}{Chromosome colour cycle for Manhattan plot}
  \item{manhattan_interspace}{Chromosome interspacing for Manhattan plot}
  \item{qqplot_col}{Colour for QQ plot}
  \item{qqplot_alpha}{Null envelope alpha value for QQ plot}
  \item{plot_fastbox}{Box height for fast approximate plotting}
  \item{dbc}{Database connection}
}
\details{
  This high level function summarizes and visualizes results
  (i.e. complete summary statistics) from a single GWAS analysis.

  Currently, a distance-pruned table of top hits, and optional
  Manhattan and QQ plots, are generated.

  The distance pruning and gene annotation is currently very slow and
  needs to be optimized.

  The Manhattan and QQ plots use an approximation to speed up both the data
  query and the plot rendering, where boxes or lines are drawn instead of
  individual points, for y-axis values between 0 and
  \code{plot_fastbox}.  Faint grey hatching is shown over this part of
  the plot to indicate that the true data are not being plotted.

  Y-axis points are truncated, but this needs to be made more clear in the
  plot.
}
\value{
  \code{gwas} returns a data frame containing the distance pruned top
  hits.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
		  