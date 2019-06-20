# gwas2() commented arguments from gwas(), may be implemented in future

#' GWAS pairwise comparison 
#' 
#' Compare results from two Genome Wide Association Studies.
#' 
#' Generates a plot showing (flipped) Z scores for all variants significant in either 
#' analysis.  Distance-based pruning is applied separately to each analysis, and 
#' index variants for first/second/both analyses are shown as red/green/yellow 
#' down-triangles/up-triangles/squares.  Using each set of index variants a GRS 
#' is constructed and tested for association in the other analysis.
#' 
#' @param analysis1 Analysis identifier for first GWAS (plot x-axis)
#' @param analysis2 Analysis identifier for second GWAS (plot y-axis)
#' @param pval1_thresh Significance threshold for first GWAS. Default: 5e-08
#' @param pval2_thresh Significance threshold for second GWAS. Default: 5e-08
#' @param prune_dist Distance-based pruning threshold (for both).  Default: 500000
#' @param dbc Database connection. Default: getOption("gtx.dbConnection", NULL)
#' @return Invisible dataframe of significant variants, with prune1/2 columns added
#' 
#' @author Toby Johnson \email{Toby.x.Johnson@gsk.com}
#' @export
gwas2 <- function(analysis1, analysis2, 
                  pval1_thresh = 5e-08, pval2_thresh = 5e-08, 
#                 maf_ge, rsq_ge, emac_ge, case_emac_ge, 
#                 gene_annotate = TRUE,
#                 plot_ymax = 30,
                  prune_dist = 500000L,
#                 manhattan_thresh = 5e-08,
#                 manhattan_col = c('#064F7C', '#6D97BD'),
#                 manhattan_interspace = 50e6,
#                 qqplot_col = '#064F7C',
#                 qqplot_alpha = 0.01,
#                 plot_fastbox = 2, 
#                 zrok = FALSE,
                  dbc = getOption("gtx.dbConnection", NULL)) {
  gtxdbcheck(dbc)

  # note, gtxanalysisdb() not being used as this idea is deprecated...
  res <- sqlWrapper(dbc,
                  sprintf('SELECT t1.chrom, t1.pos, t1.ref, t1.alt,
                             t1.beta AS beta1, t1.se AS se1, t1.pval AS pval1,
                             t2.beta AS beta2, t2.se AS se2, t2.pval AS pval2
                             FROM gwas_results AS t1
                             JOIN gwas_results AS t2
                             USING (chrom, pos, ref, alt)
                             WHERE (t1.pval <= %s OR t2.pval <= %s) 
                                 AND t1.pval IS NOT NULL AND t2.pval IS NOT NULL
                                 AND %s AND %s
                             ORDER BY chrom, pos, ref, alt',
                          sanitize1(pval1_thresh, type = 'double'),
                          sanitize1(pval2_thresh, type = 'double'),
                          gtx:::gtxwhat(analysis1 = analysis1, tablename = 't1'),
                          gtx:::gtxwhat(analysis1 = analysis2, tablename = 't2')),
           uniq = FALSE)
  res <- data.table::data.table(res)

  # compute flip +/-1 so that effect on analysis1 always positive
  bs = sign(res$beta1)
  bz = which(bs == 0L)
  if (length(bz) > 0) {
    warning('Randomly flipping signs of beta2 for ', length(bz), ' variant(s) with beta1==0.')
    bs[bz] <- ifelse(runif(length(bz)) < 0.5, -1L, 1L)
  }
  stopifnot(all(bs == -1L | bs == 1L, na.rm = TRUE))
  res[ , flip := bs]

  # prune by distance for highlighting in plot and GRS MR analyses
  # see UB-254, (1) we should prune by max(z2) or max(abs(z))
  #   because arbitrary tie-breaking when pval=0 gives visually alarming results
  #  (2) we should not have to rename a column to pval here, but pass column 
  #   name to prune.distance
  res_thresh1 <- which(res$pval1 <= pval1_thresh)
  res[ , pval := pval1]
  res_sigs1 <- prune.distance(res[res_thresh1, ], surround = prune_dist, sorted = TRUE)
  res[ , prune1 := FALSE]
  res[res_thresh1[res_sigs1$row], prune1 := TRUE]
  
  res_thresh2 <- which(res$pval2 <= pval1_thresh)
  res[ , pval := pval2]
  res_sigs2 <- prune.distance(res[res_thresh2, ], surround = prune_dist, sorted = TRUE)
  res[ , prune2 := FALSE]
  res[res_thresh2[res_sigs2$row], prune2 := TRUE]

  pdesc1 <- gtxanalysis_label(analysis = analysis1, nlabel = FALSE, dbc = dbc)
  pdesc2 <- gtxanalysis_label(analysis = analysis2, nlabel = FALSE, dbc = dbc)
  
  plot.new()
  plot.window(c(0, max(c(qnorm(pval1_thresh/2., lower.tail = FALSE), with(res, flip*beta1/se1)), na.rm = TRUE)),
              range(c(0, qnorm(pval2_thresh/2.), qnorm(pval2_thresh/2., lower.tail = FALSE), with(res, flip*beta2/se2)), na.rm = TRUE))
  # forces (0,0) and pval thresholds to be included in plot
  
  polyx <- c(0, rep(qnorm(pval1_thresh/2., lower.tail = FALSE), 2), 0)
  polyy <- c(-1, -1, 1, 1)*qnorm(pval2_thresh/2., lower.tail = FALSE)
  polygon(polyx, polyy, 
          density = NA, col = rgb(.83, .83, .83, 1), border = rgb(.83, .83, .83, 1), lwd = 9*.5)
  polygon(polyx, polyy, 
          density = 10, angle = 67.5, col = rgb(.75, .75, .75, .5), border = NA)

  # both=22; prune1 only=25; prune2 only=24; neither=21 [border] or 19 [no border]
  with(res[(!prune1 & !prune2)], points(flip*beta1/se1, 
               flip*beta2/se2,
               pch = 19, 
               col = rgb(.67, .67, .67, .1)))
  with(res[(prune1 & !prune2)], points(flip*beta1/se1, 
                  flip*beta2/se2,
     pch = 25, 
     bg = rgb(1, 0, 0, .5),
     col = rgb(.33, .33, .33, .5)))
  with(res[(prune2 & !prune1)], points(flip*beta1/se1, 
                    flip*beta2/se2,
                    pch = 24, 
                    bg = rgb(0, 1, 0, .5),
                    col = rgb(.33, .33, .33, .5)))
  with(res[(prune2 & prune1)], points(flip*beta1/se1, 
                                       flip*beta2/se2,
                                       pch = 22, 
                                       bg = rgb(1, 1, 0, .5),
                                       col = rgb(.33, .33, .33, .5)))
  abline(h = 0)
  abline(v = 0)
  axis(1);axis(2, las = 1); box()
  mtext.fit(xlab = paste(pdesc1, 'association Z score'),
            ylab = paste(pdesc2, 'association Z score'))
  mtext(paste0('Y~grs(X) p=', signif(with(res[(prune1)], grs.summary(beta1, beta2, se2, n = 1000))$pval, 3)),
        side = 3, line = 0)
  mtext(paste0('X~grs(Y) p=', signif(with(res[(prune2)], grs.summary(beta2, beta1, se1, n = 1000))$pval, 3)),
        side = 3, line = 1)
  
  return(invisible(res)) 
}
  # FIXME get analysis metadata, only needed n for for R2rs calculation...
  

