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
                    bf = if (nv > 0) c(abf1[1]*abf2[1], 
                      sum(abf1[-1])*abf2[1]/nv, 
                      abf1[1]*sum(abf2[-1])/nv, 
                      (sum(abf1[-1])*sum(abf2[-1]) - sum(abf1[-1]*abf2[-1]))/(nv*(nv - 1)), 
                      sum(abf1[-1]*abf2[-1])/nv) else rep(NA, 5))
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
                  chrom, pos_start, pos_end, pos, 
                  hgncid, ensemblid, rs, surround = 500000,
                  entity, entity1, entity2,
                  style = 'Z', 
                  dbc = getOption("gtx.dbConnection", NULL)) {
  gtxdbcheck(dbc)

  ## Determine genomic region from arguments
  xregion <- gtxregion(chrom = chrom, pos_start = pos_start, pos_end = pos_end, pos = pos, 
                       hgncid = hgncid, ensemblid = ensemblid, rs = rs, surround = surround,
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

## analysis1 must have entities, analysis2 must not
## note that surround=0 is a sensible default because
## of the region expansion performed within this function
multicoloc.data <- function(analysis1, analysis2,
                            chrom, pos_start, pos_end, pos, 
                            hgncid, ensemblid, rs, surround = 0,
##                            entity, entity1, entity2,
##                            style = 'Z', 
                            hard_clip = FALSE, 
                            dbc = getOption("gtx.dbConnection", NULL)) {
  gtxdbcheck(dbc)
  ## note there is a *niche* case where analysis2 would have an entity
  ## e.g. coloc ABC123 expression against expression of all other nearby entities

  
  ## Determine genomic region from arguments
  xregion <- gtxregion(chrom = chrom, pos_start = pos_start, pos_end = pos_end, pos = pos, 
                       hgncid = hgncid, ensemblid = ensemblid, rs = rs, surround = surround,
                       dbc = dbc)
  chrom = xregion$chrom
  pos_start = xregion$pos_start
  pos_end = xregion$pos_end

  ## Currently only works if analysis1 all in the same db table
  db1 <- sanitize1(unique(sapply(analysis1, gtxanalysisdb)), type = 'alphanum')
  
#  ## substitute generic entity for entity1 and entity2 if needed
#  if (missing(entity1) && !missing(entity)) entity1 <- entity
#  if (missing(entity2) && !missing(entity)) entity2 <- entity

#  ## Determine entity, if required, for each analysis
#  xentity1 <- gtxentity(analysis1, entity = entity1, hgncid = hgncid, ensemblid = ensemblid)
#  xentity2 <- gtxentity(analysis2, entity = entity2, hgncid = hgncid, ensemblid = ensemblid)

  ## We want to include complete association statistics for all entities that
  ## partly or fully overlap the query region.  (This cannot be done
  ## by adding a buffer region since entities are of variable size, and have
  ## association statistics that extend by unknown intervals around the
  ## entity [GTEx 1Mb cis-regions are not guaranteed].)
  ##
  ## Therefore we run a series of queries, to
  ## 1. Find all entities with >=1 association statistic within the overlap query region
  ## 2. Find the interval that includes all association statistics for those entities
  ## 3. Do a coloc query for this region (with WHERE ... AND entity=)
  ##
  ## Note that for efficiency regions we want query 3. to be by a defined physical region
  ## rather than directly selecting on WHERE entity IN ...
  
  gtxlog('Query region is chr', chrom, ':', pos_start, '-', pos_end, 
         ' (', prettyNum(pos_end - pos_start, big.mark = ',', scientific = FALSE), ' bp)')

  eq <- sqlWrapper(dbc, 
                    sprintf('SELECT 
                                 DISTINCT feature
                             FROM %s.gwas_results
                             WHERE
                                 %s AND %s ;',
                            db1, 
                            gtxwhat(analysis = analysis1), 
                            gtxwhere(chrom = chrom, pos_ge = pos_start, pos_le = pos_end)),
                    uniq = FALSE)$feature
  ## FIXME this may return zero rows, should handle gracefully
  gtxlog('Query region includes association statistics for ', length(eq), ' entities')

  if (!hard_clip) {
      ep <- sqlWrapper(dbc, 
                       sprintf('SELECT 
                                 min(pos) as minpos, max(pos) as maxpos
                             FROM %s.gwas_results
                             WHERE
                                 %s AND %s ;',
                             db1, 
                             gtxwhat(analysis = analysis1),
                             gtxwhere(chrom = chrom, entity = eq)))
      pos_start <- ep$minpos
      pos_end <- ep$maxpos
      gtxlog('Expanded region is chr', chrom, ':', pos_start, '-', pos_end,
             ' (', prettyNum(pos_end - pos_start, big.mark = ',', scientific = FALSE), ' bp)')
  }

  ## We use a (INNER) JOIN and silently drop rows that don't match
  ## FIXME if hard_clip can make this query run faster by not selecting on entity
  t0 <- as.double(Sys.time())
  res <- sqlWrapper(dbc,
                    sprintf('SELECT 
                                 t1.analysis AS analysis1, t1.entity AS entity1,
                                 t1.beta AS beta1, t1.se AS se1, 
                                 t2.beta AS beta2, t2.se AS se2 
                             FROM 
                                 (SELECT
                                      chrom, pos, ref, alt, analysis, feature AS entity, beta, se
                                  FROM %s.gwas_results 
                                  WHERE
                                      %s AND %s AND %s
                                 ) AS t1 
                             JOIN 
                                 (SELECT 
                                      chrom, pos, ref, alt, beta, se
                                  FROM %s.gwas_results 
                                  WHERE 
                                      analysis=\'%s\' 
                                      AND %s 
                                 ) AS t2
                             USING (chrom, pos, ref, alt);',
                            db1, 
                            gtxwhat(analysis = analysis1),
                            gtxwhere(chrom = chrom, pos_ge = pos_start, pos_le = pos_end),
                            gtxwhere(chrom = chrom, entity = eq), 
                            sanitize(gtxanalysisdb(analysis2), type = 'alphanum'), # may not require sanitation
                            sanitize(analysis2, type = 'alphanum'),
                            gtxwhere(chrom = chrom, pos_ge = pos_start, pos_le = pos_end)
                            ),
                    uniq = FALSE) # expect >=1 rows
  t1 <- as.double(Sys.time())
  gtxlog('Results query returned ', nrow(res), ' rows in ', round(t1 - t0, 3), 's.')
  return(res)
}

multicoloc <- function(analysis1, analysis2,
                       chrom, pos_start, pos_end, pos, 
                       hgncid, ensemblid, rs, surround = 0,
                       style = 'heatplot', 
                       dbc = getOption("gtx.dbConnection", NULL)) {
  gtxdbcheck(dbc)

  ## get summary stats
  ss <- multicoloc.data(analysis1 = analysis1, analysis2 = analysis2,
                         chrom = chrom, pos_start = pos_start, pos_end = pos_end, pos = pos, 
                         hgncid = hgncid, ensemblid = ensemblid, rs = rs,
                         surround = surround,
                         dbc = dbc)

  res <- sqlWrapper(dbc, 
                    sprintf('SELECT * FROM genes WHERE %s ORDER BY pos_start;',
                            gtxwhere(ensemblid = unique(ss$entity))), # FIXME not guaranteed entity_type is ENSG
                    uniq = FALSE)
  res$entity <- res$ensemblid # FIXME not guaranteed entity_type
  analyses <- unique(ss$analysis1)

  t0 <- as.double(Sys.time())  
  for (this_analysis in analyses) {
      resc <- do.call(rbind,
                      lapply(unique(res$entity), function(this_entity) {
                          return(subset(coloc.fast(subset(ss, analysis1 == this_analysis & entity1 == this_entity))$results,
                                        hypothesis == 'H4')$posterior)
                      }))
      colnames(resc) <- paste0('Hxy', '_', this_analysis)
      res <- cbind(res, resc)
  }
  t1 <- as.double(Sys.time())
  gtxlog('Colocalization analyzed for ', sum(!is.na(res[ , paste0('Hxy', '_', analyses)])), ' pairs in ', round(t1 - t0, 3), 's')

  if (identical(style, 'none')) {
      ## do nothing
  } else if (identical(style, 'heatplot')) {
      zmat <- as.matrix(res[ , paste0('Hxy', '_', analyses)])
      colnames(zmat) <- analyses
      rownames(zmat) <- with(res, ifelse(hgncid != '', as.character(hgncid), as.character(ensemblid))) # FIXME will this work for all entity types
      thresh_z <- .1*max(zmat, na.rm = TRUE) # threshold
      zmat <- zmat[ , order(apply(zmat, 2, function(x) if (any(x >= thresh_z, na.rm = TRUE)) max(x, na.rm = TRUE) else NA), na.last = NA), drop = FALSE]
      multicoloc.plot(zmat)
  } else {
      stop('unknown style [ ', style, ' ]')
  }
  
  return(res)
}

## Input, a matrix of z values with analysis as column names and entity as row names 
multicoloc.plot <- function(zmat, 
                            dbc = getOption("gtx.dbConnection", NULL)) {
    gtxdbcheck(dbc)
    
    ## Query plot labels for analyses
    label_y <- sqlWrapper(dbc, 
                         sprintf('SELECT analysis, label FROM analyses WHERE %s',
                                 gtxwhat(analysis = colnames(zmat))),
                         uniq = FALSE)
    label_y <- label_y$label[match(colnames(zmat), label_y$analysis)]
    ## label_y <- ifelse(!is.na(label_y), label_y, colnames(zmat)) # fall back to analysis if label lookup failed
    
    plot.new()
    x_labelmax <- .4 # max fraction of total x space to use for analysis descriptions
    y_linesep <- 2. # spacing desired as multiple of strheight()
    
    cex_ylab <- 1.
    while (TRUE) {
        y_used <- sum(strheight(label_y, cex = cex_ylab)*y_linesep)
        x_used <- max(strwidth(label_y, cex = cex_ylab))
        if (y_used <= 1. && x_used <= x_labelmax) break
        cex_ylab <- cex_ylab*min(1./y_used, x_labelmax/x_used)
    }
    x_labeluse <- x_used
    
    cex_values <- 1.
    while (TRUE) {
        y_used <- strheight('000', cex = cex_values)*ncol(zmat)*y_linesep
        x_used <- strwidth('000', cex = cex_values)*nrow(zmat)
        if (y_used <= 1. && x_used <= (1. - x_labeluse)) break
        cex_values <- cex_values*min(1./y_used, (1. - x_labeluse)/x_used)
    }

    plot.window(c(-x_labeluse/(1. - x_labeluse), 1.)*(nrow(zmat) + .5), c(.5, ncol(zmat) + .5))
    abline(v = 0)
    
    image(x = 1:nrow(zmat), y = 1:ncol(zmat),
          z = zmat,
          zlim = c(0, 1),
          col = rgb(1, 100:0/100, 100:0/100,),
          add = TRUE) # should add options for different colour scalings

    text(0, 1:ncol(zmat), label_y, pos = 2, cex = cex_ylab)
    for (idx in 1:ncol(zmat)) {
        zvals <- as.integer(round(zmat[, idx]*100))
        text(1:nrow(zmat), idx, ifelse(!is.na(zvals), sprintf('%02i', zvals), ''), cex = cex_values)
    }
    axis(1, at = 1:nrow(zmat),
         labels = rownames(zmat),
         las = 2, cex.axis = .5, font = 3)
    box()

    return(invisible(NULL))
}
