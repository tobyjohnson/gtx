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

  bs = sign(res$beta)
  bz = which(bs == 0L)
  if (length(bz) > 0) {
    warning('Randomly flipping signs of beta2 for ', length(bz), ' variant(s) with beta1==0.')
    bs[bz] <- ifelse(runif(length(bz)) < 0.5, -1L, 1L)
  }
  stopifnot(all(bs == -1L | bs == 1L, na.rm = TRUE))
  res[ , flip := bs]

  res_thresh1 <- which(res$pval1 <= pval1_thresh)
  res[ , pval := pval1]
  res_sigs1 <- prune.distance(res[res_thresh1, ], surround = prune_dist, sorted = TRUE)
  res[ , prune1 := FALSE]
  res[res_thresh1[res_sigs1$row], prune1 := TRUE]
  res1 <- res[res_thresh1[res_sigs1$row], ] # after distance-based pruning for wrt pval for analysis1

  pdesc1 <- gtxanalysis_label(analysis = analysis1, nlabel = FALSE, dbc = dbc)
  pdesc2 <- gtxanalysis_label(analysis = analysis2, nlabel = FALSE, dbc = dbc)
  
  plot.new()
  plot.window(c(0, max(c(qnorm(pval1_thresh/2., lower.tail = FALSE), with(res, flip*beta/se)), na.rm = TRUE)),
              range(c(0, qnorm(pval2_thresh/2.), qnorm(pval2_thresh/2., lower.tail = FALSE), with(res, flip*beta2/se2)), na.rm = TRUE))
  # forces (0,0) and pval thresholds to be included in plot
  
  polygon(c(0, rep(qnorm(pval1_thresh/2., lower.tail = FALSE), 2), 0), 
          c(-1, -1, 1, 1)*qnorm(pval2_thresh/2., lower.tail = FALSE), 
          density = 20, angle = 45, col = rgb(.5, .5, .5, 1), border = NA)
  
    with(res, points(flip*beta/se, 
               flip*beta2/se2,
               pch = 19, 
               col = rgb(.67, .67, .67, .1)))
  with(res1, points(flip*beta/se, 
                  flip*beta2/se2,
     pch = 21, 
     bg = rgb(1, 0, 0, .5),
     col = rgb(.33, .33, .33, .5)))
  abline(h = 0)
  abline(v = 0)
  axis(1);axis(2, las = 1); box()
  mtext.fit(main = paste0('GRS p=', signif(with(res1, grs.summary(beta, beta2, se2, n = 1000))$pval, 3)),
            xlab = paste(pdesc1, 'association Z score'),
            ylab = paste(pdesc2, 'association Z score'))

   return(invisible(NULL)) 
}
  # FIXME get analysis metadata, only needed n for for R2rs calculation...
  

