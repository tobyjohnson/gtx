## Colocalisation analysis, implementing method of Giambartolomei et al. 2014
coloc.fast <- function(data, rounded = 6,
                       priorsd1 = 1, priorsd2 = 1, priorc1 = 1e-4, priorc2 = 1e-4, priorc12 = 1e-5, 
                       beta1 = 'beta1', se1 = 'se1', beta2 = 'beta2', se2 = 'se2') {
  abf1 <- abf.Wakefield(data[[beta1]], data[[se1]], priorsd1, log = TRUE)
  abf2 <- abf.Wakefield(data[[beta2]], data[[se2]], priorsd2, log = TRUE)
  w <- which(!is.na(abf1) & !is.na(abf2))
  ## Fill in zeros not filter out, thanks Karl and Karsten!!!!
  nv <- length(w)
  abf1 <- norm1(c(0, abf1[w]), log = TRUE)
  abf2 <- norm1(c(0, abf2[w]), log = TRUE)
  res <- data.frame(hypothesis = paste0("H", 0:4),
                    label = c("No association",
                      "One variant associated with phenotype 1 only",
                      "One variant associated with phenotype 2 only",
                      "Two variants separately associated with phenotypes 1 and 2",
                      "One variant associated with phenotypes 1 and 2"),
                    prior = norm1(c(1, priorc1*nv, priorc2*nv, priorc1*priorc2*nv*(nv - 1), priorc12*nv)),
                    bf = c(abf1[1]*abf2[1], 
                      sum(abf1[-1])*abf2[1]/nv, 
                      abf1[1]*sum(abf2[-1])/nv, 
                      (sum(abf1[-1])*sum(abf2[-1]) - sum(abf1[-1]*abf2[-1]))/(nv*(nv - 1)), 
                      sum(abf1[-1]*abf2[-1])/nv))
  res$bf <- res$bf/res$bf[1]
  res$posterior <- norm1(res$prior*res$bf)
  if (is.finite(rounded)) {
    res$posterior = round(res$posterior, rounded)
  }
  return(list(results = res, nvariants = length(w)))
}

## TO DO
## Function rename rationalization
##   regionplot.region -> gtxregion
## Sanitation for required length 1
## Wrapper for SQL queries that must return data frame of 1 or >=1 rows



coloc <- function(analysis1, analysis2,
                  chrom, pos_start, pos_end,
                  hgncid, ensemblid, surround = 500000,
                  entity, entity1, entity2,
                  style = 'Z', 
                  dbc = getOption("gtx.dbConnection", NULL)) {
  gtxdbcheck(dbc)

  ## Determine genomic region from arguments
  xregion <- gtxregion(chrom, pos_start, pos_end,
                               hgncid, ensemblid, surround,
                               dbc = dbc)
  chrom = xregion$chrom
  pos_start = xregion$pos_start
  pos_end = xregion$pos_end

  ## substitute generic entity for entity1 and entity2 if needed
  if (missing(entity1) && !missing(entity)) entity1 <- entity
  if (missing(entity2) && !missing(entity)) entity2 <- entity

  ## Determine entity, if required, for each analysis
  xentity1 <- gtxentity(analysis1, entity = entity1, hgncid = hgncid, ensemblid = ensemblid)
  xentity2 <- gtxentity(analysis2, entity = entity2, hgncid = hgncid, ensemblid = ensemblid)

  ## Get association statistics
  res <- sqlWrapper(dbc, 
                    sprintf('SELECT 
                                 t1.beta AS beta1, t1.se AS se1, 
                                 t2.beta AS beta2, t2.se AS se2 
                             FROM 
                                 (SELECT
                                      chrom, pos, ref, alt, beta, se, pval 
                                  FROM %s.gwas_results 
                                  WHERE
                                      analysis=\'%s\' 
                                      AND %s 
                                      %s 
                                 ) AS t1 
                                 FULL JOIN 
                                 (SELECT 
                                      chrom, pos, ref, alt, beta, se, pval 
                                  FROM %s.gwas_results 
                                  WHERE 
                                      analysis=\'%s\' 
                                      AND %s 
                                      %s
                                 ) AS t2
                                 USING (chrom, pos, ref, alt);',
                            sanitize(gtxanalysisdb(analysis1), type = 'alphanum'), # may not require sanitation
                            sanitize(analysis1, type = 'alphanum'),
                            gtxwhere(chrom, pos_ge = pos_start, pos_le = pos_end),
                            if (!is.null(xentity1)) sprintf(' AND feature=\'%s\'', xentity1$entity) else '', # FIXME will change to entity
                            sanitize(gtxanalysisdb(analysis2), type = 'alphanum'), # may not require sanitation
                            sanitize(analysis2, type = 'alphanum'),
                            gtxwhere(chrom, pos_ge = pos_start, pos_le = pos_end),
                            if (!is.null(xentity2)) sprintf(' AND feature=\'%s\'', xentity2$entity) else ''  # FIXME will change to entity
                            ),
                    uniq = FALSE) # expect >=1 rows

  gtxlog('coloc query returned ', nrow(res), ' rows')

  resc <- coloc.fast(res)

  pdesc1 <- sqlWrapper(dbc,
                       sprintf('SELECT description FROM analyses WHERE analysis = \'%s\';',
                               sanitize(analysis1, type = 'alphanum')))$description
  if (!is.null(xentity1)) pdesc1 <- paste(xentity1$entity_label, pdesc1)
  pdesc2 <- sqlWrapper(dbc,
                       sprintf('SELECT description FROM analyses WHERE analysis = \'%s\';',
                               sanitize(analysis2, type = 'alphanum')))$description
  if (!is.null(xentity2)) pdesc2 <- paste(xentity2$entity_label, pdesc2)
 
  if (identical(style, 'Z')) {
      ## would like to draw one sided plot, but unclear what to do when sign(beta1)==0. HMMMM FIXME
      with(res, {
          plot(beta1/se1,
               beta2/se2,
               pch = 21, bg = rgb(.67, .67, .67, .5), col = rgb(.33, .33, .33, .5), cex = 1,
               ann = FALSE)
	  abline(h = 0)
	  abline(v = 0)
          mtext.fit(main = paste0('H', c('0', 'x', 'y', 'x,y', 'xy'), '=', round(resc$results$posterior*100), '%', collapse = ', '),
                    xlab = paste(pdesc1, 'association Z score'),
                    ylab = paste(pdesc2, 'association Z score'))
	  mtext(paste0('colocalization at chr', chrom, ':', pos_start, '-', pos_end), 3, 3)
      })
  } else if (identical(style, 'beta')) {
      ## would like to draw one sided plot, but unclear what to do when sign(beta1)==0. HMMMM FIXME
      with(res, {
          plot(beta1,
               beta2,
               pch = 21, bg = rgb(.67, .67, .67, .5), col = rgb(.33, .33, .33, .5), cex = 1,
               ann = FALSE)
	  abline(h = 0)
	  abline(v = 0)
          mtext.fit(main = paste0('H', c('0', 'x', 'y', 'x,y', 'xy'), '=', round(resc$results$posterior*100), '%', collapse = ', '),
                    xlab = paste(pdesc1, 'association effect size'),
                    ylab = paste(pdesc2, 'association effect size'))
	  mtext(paste0('colocalization at chr', chrom, ':', pos_start, '-', pos_end), 3, 3)
      })
  } else if (identical(style, 'none')) {
      # do not draw plot
  } else {
      stop('Unrecognized style; must be Z, beta, or none')
  }

  return(resc)
}
