## Colocalisation analysis, implementing method of Giambartolomei et al. 2014
##
## changed API to pass in vectors instead of dataframe
#' @export
coloc.fast <- function(beta1, se1, beta2, se2, 
                       priorsd1 = 1, priorsd2 = 1, priorc1 = 1e-4, priorc2 = 1e-4, priorc12 = 1e-5, 
                       rounded = 6) {
  abf1 <- abf.Wakefield(beta1, se1, priorsd1, log = TRUE)
  abf2 <- abf.Wakefield(beta2, se2, priorsd2, log = TRUE)
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
  res$bf <- res$bf/max(res$bf) # was: res$bf/res$bf[1] # caused division by zero for strongly colocalized signals
  res$posterior <- norm1(res$prior*res$bf)
  ## compute model averaged effect size ratios
  mw <- abf1[-1]*abf2[-1] # model weights
  alpha12 <- sum(beta1[w]/beta2[w]*mw)/sum(mw)
  alpha21 <- sum(beta2[w]/beta1[w]*mw)/sum(mw)
  if (is.finite(rounded)) {
    res$posterior = round(res$posterior, rounded)
  }
  return(list(results = res, nvariants = length(w), alpha12 = alpha12, alpha21 = alpha21))
}

coloc.compute <- function(data,
                          priorsd1 = 1, priorsd2 = 1, priorc1 = 1e-4, priorc2 = 1e-4, priorc12 = 1e-5,
                          join_type = 'inner',
                          summary_only = FALSE) {
    ## coloc_compute replaces lnabf that are NA by -Inf , hance zero probability.
    ## Can be used two ways:
    ##   1. Set to *match* the type of join used in the preceding db query, so that rows
    ##      missing from one side of the query are handled appropriately.
    ##   2. Use a broader db query e.g. full/outer join, then set this option to
    ##      use the results *as if* a narrower query had been used.  E.g. previous
    ##      behaviour was always to do full/outer join db query and then drop NAs in
    ##      coloc.fast, hence as if a normal/inner join query had been done.
    stopifnot(is.data.frame(data))
    stopifnot(all(c('beta1', 'se1', 'beta2', 'se2') %in% names(data)))
    stopifnot(join_type %in% c('inner', 'left', 'right', 'outer'))
    data <- within(data, {
        lnabf1 <- abf.Wakefield(beta1, se1, priorsd1, log = TRUE)
        lnabf2 <- abf.Wakefield(beta2, se2, priorsd2, log = TRUE)
        inc1 <- if (join_type == 'inner' | join_type == 'left') !is.na(lnabf1) else TRUE
        inc2 <- if (join_type == 'inner' | join_type == 'right') !is.na(lnabf2) else TRUE
        inc <- inc1 & inc2
    })
    data <- data[which(data$inc), ]
    data <- within(data, {
        inc1 <- NULL
        inc2 <- NULL
        inc <- NULL
    })        
    abf1 <- norm1(c(0, data$lnabf1), log = TRUE)
    abf2 <- norm1(c(0, data$lnabf2), log = TRUE)
    abf1[is.na(abf1)] <- 0
    abf2[is.na(abf2)] <- 0
    data$abf1 <- abf1[-1]
    data$abf2 <- abf2[-1]
    attr(data, 'nullabf1') <- abf1[1]
    attr(data, 'nullabf2') <- abf2[1]

    nv <- nrow(data)
    ## compute colocalization model probabilities
    prior <- norm1(c(1, priorc1*nv, priorc2*nv, priorc1*priorc2*nv*(nv - 1), priorc12*nv))
    names(prior) <- c('H0', 'H1', 'H2', 'H1,2', 'H12')
    bf = if (nv > 0) c(abf1[1]*abf2[1], 
                       sum(abf1[-1])*abf2[1]/nv, 
                       abf1[1]*sum(abf2[-1])/nv, 
                       (sum(abf1[-1])*sum(abf2[-1]) - sum(abf1[-1]*abf2[-1]))/(nv*(nv - 1)), 
                       sum(abf1[-1]*abf2[-1])/nv) else rep(NA, 5)
    bf <- bf/max(bf)
    names(bf) <- names(prior)
    posterior <- norm1(prior*bf)
    names(posterior) <- names(prior)

    ## compute model averaged effect size ratios
    mw <- data$abf1*data$abf2 # model weights
    alpha12 <- sum(mw*data$beta1/data$beta2)/sum(mw)
    alpha21 <- sum(mw*data$beta2/data$beta1)/sum(mw)

    if (summary_only) {
        return(list(nv = nv, prior = prior, bf = bf, posterior = posterior, alpha12 = alpha12, alpha21 = alpha21))
    } else {
        attr(data, 'coloc') <- list(nv = nv, prior = prior, bf = bf, posterior = posterior, alpha12 = alpha12, alpha21 = alpha21)
        ## should we include the assumed priors?
        return(data)
    }

    ## should never:
    return(invisible(NULL))
}

## TO DO
## Function rename rationalization
##   regionplot.region -> gtxregion
## Sanitation for required length 1
## Wrapper for SQL queries that must return data frame of 1 or >=1 rows

#' @export
coloc.data <- function(analysis1, analysis2,
                  chrom, pos_start, pos_end, pos, 
                  hgncid, ensemblid, rs, surround = 500000,
                  entity, entity1, entity2,
                  dbc = getOption("gtx.dbConnection", NULL), ...) {
  gtxdbcheck(dbc)

  ## Determine genomic region from arguments
  xregion <- gtxregion(chrom = chrom, pos_start = pos_start, pos_end = pos_end, pos = pos, 
                       hgncid = hgncid, ensemblid = ensemblid, rs = rs, surround = surround,
                       dbc = dbc)

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
                                  FROM %sgwas_results 
                                  WHERE
                                      %s AND %s %s
                                 ) AS t1 
                                 FULL JOIN 
                                 (SELECT 
                                      chrom, pos, ref, alt, beta, se, pval 
                                  FROM %sgwas_results 
                                  WHERE 
                                      %s AND %s %s
                                 ) AS t2
                                 USING (chrom, pos, ref, alt);',
                            gtxanalysisdb(analysis1), 
                            gtxwhat(analysis1 = analysis1), # analysis1= argument allows only one analysis
                            gtxwhere(chrom = xregion$chrom, pos_ge = xregion$pos_start, pos_le = xregion$pos_end),
                            if (!is.null(xentity1)) sprintf(' AND feature=\'%s\'', xentity1$entity) else '', # FIXME will change to entity
                            gtxanalysisdb(analysis2), 
                            gtxwhat(analysis1 = analysis2), # analysis1= argument allows only one analysis, different to arguments analysis1 and analysis2
                            gtxwhere(chrom = xregion$chrom, pos_ge = xregion$pos_start, pos_le = xregion$pos_end),
                            if (!is.null(xentity2)) sprintf(' AND feature=\'%s\'', xentity2$entity) else ''  # FIXME will change to entity
                            ),
                    uniq = FALSE) # expect >=1 rows

  attr(res, 'region') <- xregion
  attr(res, 'entity1') <- xentity1 # if NULL then no attr is set
  attr(res, 'entity2') <- xentity2 # if NULL then no attr is set
  
  return(res)
}


#'
#'
#' Colocalization analysis.
#'
#' Colocalization analysis.
#'
#' This high level function conducts a colocalization analysis, using
#' summary statistics for association with two traits, across a region of
#' the genome.  The two sets of summary statistics are specified using
#' the \code{analysis1} and \code{analysis2} arguments.  Where one or
#' both contain summary statistics for multiple entities (e.g. from eQTL or
#' pQTL analyses), the desired entities must be specified (see below).
#'
#' Note that when using a \code{hgncid} or \code{ensemblid} gene
#' identifier to specify the region from which to use summary statistics,
#' the default \code{surround=500000} will \emph{not} include the full
#' cis eQTL region as usually specified.
#'  
#' The region of interest can be specified in several different ways.
#' The region can be supplied as physical coordinates using the arguments
#' \code{chrom}, \code{pos_start} and \code{pos_end}.  Alternatively, the
#' region can be centered on a gene of interest, using either the
#' \code{hgncid} or \code{emsemblid} argument, and the size of region
#' around the gene can be modified using the \code{surround} argument.
#' Note that the primary purpose of gene-identifying arguments
#' \code{hgncid} or \code{ensemblid} is to specify the genomic region of
#' interest (and thus the set of the variants to analyse).  It is only a
#' secondary purpose that the entity for eQTL or pQTL analyses will be
#' inferred from \code{hgncid} or \code{ensemblid}, if no explicit
#' \code{entity} argument is given.
#'
#' Entities are used to distinguish genomic features, where a single set
#' analysis includes summary statistics, for each variant, for
#' associations with one or more entities.  E.g. in an eQTL analysis,
#' each transcript or gene is an entity, and a single typical variant
#' will have summary statistics for associations with multiple
#' transcripts or genes.  If either of the analyses specified by
#' \code{analysis1} and \code{analysis2} have results separated by
#' entity, then the arguments \code{entity1} and \code{entity2} are used
#' to specify the desired entity from each.  If either \code{entity1} or
#' \code{entity2} is missing, the argument \code{entity} is used
#' instead.  (This mechanism facilitates e.g. colocalization between analyses
#' for the same transcript between two different eQTL datasets.)  If the
#' argument \code{entity} is also missing, the function attempts to infer
#' a suitable entity from the \code{hgncid} or \code{ensemblid}
#' arguments.  (This leads to sensible default behaviour, and facilitates
#' the most common use case of centering the genomic region of interest
#' on the entity being analysed in an eQTL or pQTL dataset.)
#'
#' The \code{style} argument can be set to \sQuote{Z} to plot Z
#' statistics for the two analyses, \sQuote{beta} to plot beta (effect
#' size) statistics for the two analyses, or \sQuote{none} to suppress
#' plotting altogether.
#'
#' @param analysis1 The key value for the first GWAS analysis to analyze
#' @param analysis2 The key value for the second GWAS analysis to analyze
#' @param chrom Character specifying chromosome
#' @param pos_start Start position of region
#' @param pos_end End position of region
#' @param hgncid HGNC identifier of gene to define region around
#' @param ensemblid ENSEMBL gene identifier to define region around
#' @param surround Distance around gene to include in region. Default: 500000
#' @param entity Identifier for an entity, for analyses of multiple entities
#' @param entity1 Identifier for an entity, for analysis
#' @param entity2 Identifier for an entity, for analysis2
#' @param style Character specifying plot style. Default: 'Z'
#' @param dbc Database connection. Default: getOption("gtx.dbConnection", NULL)
#' 
#' @return
#'  a data frame containing the result of the
#'  colocalization analysis, see \code{\link{coloc.fast}} for details.
#'  The plot is generated as a side effect.
#'
#' @author Toby Johnson \email{Toby.x.Johnson@gsk.com}
#' @export

coloc <- function(analysis1, analysis2,
                  chrom, pos_start, pos_end, pos, 
                  hgncid, ensemblid, rs, surround = 500000,
                  entity, entity1, entity2,
                  style = 'Z', 
                  priorsd1 = 1, priorsd2 = 1, priorc1 = 1e-4, priorc2 = 1e-4, priorc12 = 1e-5,
                  join_type = 'inner',
                  dbc = getOption("gtx.dbConnection", NULL)) {
  gtxdbcheck(dbc)

  ## query variant-level summary stats
  res <- coloc.data(analysis1 = analysis1, analysis2 = analysis2,
                    chrom = chrom, pos_start = pos_start, pos_end = pos_end, pos = pos,
                    hgncid = hgncid, ensemblid = ensemblid, rs = rs, surround = surround,
                    entity = entity, entity1 = entity1, entity2 = entity2,
                    dbc = dbc)  
  gtxlog('coloc.data query returned ', nrow(res), ' rows')

  pdesc1 <- gtxanalysis_label(analysis = analysis1, entity = attr(res, 'entity1'), nlabel = FALSE, dbc = dbc)
  pdesc2 <- gtxanalysis_label(analysis = analysis2, entity = attr(res, 'entity2'), nlabel = FALSE, dbc = dbc)
  
  ## compute coloc probabilities
  resc <- coloc.compute(res,
                        priorsd1 = priorsd1, priorsd2 = priorsd2,
                        priorc1 = priorc1, priorc2 = priorc2, priorc12 = priorc12,
                        join_type = join_type)

  ## compute extra variables and symbol/colour for plotting
  if ('Z' %in% style || 'Z1' %in% style) {
      # compute Z statistics if needed for plots
      resc <- within(resc, { z1 <- beta1/se1 ; z2 <- beta2/se2 })
  }
  if ('beta1' %in% style || 'Z1' %in% style) {
      # prepare for sign flipping if needed for half-scatter plots
      bs = sign(resc$beta1)
      bz = which(bs == 0L)
      if (length(bz) > 0) {
          warning('Randomly flipping signs of beta2 for ', length(bz), ' variant(s) with beta1==0.')
          bs[bz] <- ifelse(runif(length(bz)) < 0.5, -1L, 1L)
      }
      stopifnot(all(bs == -1L | bs == 1L, na.rm = TRUE))
      resc$flip <- bs
  }

  ## tidy this up, put twochannel in gtxbase; more carefully name the plot style variables
  twochannel <- function(x, y, grey = 0.67, alpha = 0.5) {
      r <- pmin(pmax(1+x*y-y, 0), 1)
      g <- pmin(pmax(1+x*y-x, 0), 1)
      b <- pmin(pmax(1+x*y-y-x, 0), 1)
      return(rgb(grey*r, grey*g, grey*b, alpha = alpha))
  }
  p1 <- gtx:::norm1(c(attr(resc, 'nullabf1'), priorc1*resc$abf1))
  p2 <- gtx:::norm1(c(attr(resc, 'nullabf2'), priorc2*resc$abf2))
  cs1 <- gtx:::credset(p1)[-1]
  cs2 <- gtx:::credset(p2)[-1]
  bg <- twochannel(x = p1[-1]/max(p1),
                   y = p2[-1]/max(p2))
  cex <- 0.75 + 0.5*p1[-1] + 0.5*p2[-1]
  pch = ifelse(cs1, ifelse(cs2, 22, 25), ifelse(cs2, 24, 21))
  
  if ('Z' %in% style) {
      with(resc, {
          plot(beta1/se1,
               beta2/se2,
               xlim = range(c(0, beta1/se1), na.rm = TRUE),
               ylim = range(c(0, beta2/se2), na.rm = TRUE), # forces (0,0) to be included in plot
               cex = cex, pch = pch, bg = bg, col = rgb(.33, .33, .33, .5), 
               ann = FALSE)
	  abline(h = 0)
	  abline(v = 0)
          mtext.fit(main = paste0('H', c('0', 'x', 'y', 'x,y', 'xy'), '=',
                                  round(attr(resc, 'coloc')$posterior*100), '%', collapse = ', '),
                    xlab = paste(pdesc1, 'association Z score'),
                    ylab = paste(pdesc2, 'association Z score'))
	  mtext(paste('colocalization at', attr(resc, 'region')$label), 3, 3) # should force to fit plot area
      })
  }
  if ('Z1' %in% style) {
      with(resc, {
          plot(flip*beta1/se1,
               flip*beta2/se2,
               xlim = c(0, max(flip*beta1/se1, na.rm = TRUE)),
               ylim = range(c(0, flip*beta2/se2), na.rm = TRUE), 
               cex = cex, pch = pch, bg = bg, col = rgb(.33, .33, .33, .5), 
               ann = FALSE)
	  abline(h = 0)
	  abline(v = 0)
          mtext.fit(main = paste0('H', c('0', 'x', 'y', 'x,y', 'xy'), '=',
                                  round(attr(resc, 'coloc')$posterior*100), '%', collapse = ', '),
                    xlab = paste(pdesc1, 'association Z score'),
                    ylab = paste(pdesc2, 'association Z score'))
	  mtext(paste('colocalization at', attr(resc, 'region')$label), 3, 3) # should force to fit plot area
      })
  }
  if ('beta' %in% style) {
      with(resc, {
          plot(beta1,
               beta2,
               xlim = range(c(0, beta1), na.rm = TRUE),
               ylim = range(c(0, beta2), na.rm = TRUE), # forces (0,0) to be included in plot
               cex = cex, pch = pch, bg = bg, col = rgb(.33, .33, .33, .5), 
               ann = FALSE)
	  abline(h = 0)
	  abline(v = 0)
          mtext.fit(main = paste0('H', c('0', 'x', 'y', 'x,y', 'xy'), '=',
                                  round(attr(resc, 'coloc')$posterior*100), '%', collapse = ', '),
                    xlab = paste(pdesc1, 'association effect size'),
                    ylab = paste(pdesc2, 'association effect size'))
          mtext(paste('colocalization at', attr(resc, 'region')$label), 3, 3) # should force to fit plot area
      })
  }
  if ('beta1' %in% style) {
      with(resc, {
          plot(flip*beta1,
               flip*beta2,
               xlim = c(0, max(flip*beta1, na.rm = TRUE)),
               ylim = range(c(0, flip*beta2), na.rm = TRUE), 
               cex = cex, pch = pch, bg = bg, col = rgb(.33, .33, .33, .5), 
               ann = FALSE)
	  abline(h = 0)
	  abline(v = 0)
          mtext.fit(main = paste0('H', c('0', 'x', 'y', 'x,y', 'xy'), '=',
                                  round(attr(resc, 'coloc')$posterior*100), '%', collapse = ', '),
                    xlab = paste(pdesc1, 'association effect size'),
                    ylab = paste(pdesc2, 'association effect size'))
          mtext(paste('colocalization at', attr(resc, 'region')$label), 3, 3) # should force to fit plot area
      })
  }

  return(resc) # in future, will return invisible res with resc as an attribute
}


## analysis1 must have entities, analysis2 must not
## note that surround=0 is a sensible default because
## of the region expansion performed within this function
#' @export
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
         ' (', prettyNum(pos_end - pos_start + 1, big.mark = ',', scientific = FALSE), ' bp)')

  eq <- sqlWrapper(dbc, 
                    sprintf('SELECT 
                                 DISTINCT feature
                             FROM %sgwas_results
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
                             FROM %sgwas_results
                             WHERE
                                 %s AND %s ;',
                             db1, 
                             gtxwhat(analysis = analysis1),
                             gtxwhere(chrom = chrom, entity = eq)))
      pos_start <- ep$minpos
      pos_end <- ep$maxpos
      gtxlog('Expanded region is chr', chrom, ':', pos_start, '-', pos_end,
             ' (', prettyNum(pos_end - pos_start + 1, big.mark = ',', scientific = FALSE), ' bp)')
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
                                  FROM %sgwas_results 
                                  WHERE
                                      %s AND %s AND %s
                                 ) AS t1 
                             JOIN 
                                 (SELECT 
                                      chrom, pos, ref, alt, beta, se
                                  FROM %sgwas_results 
                                  WHERE 
                                      %s AND %s 
                                 ) AS t2
                             USING (chrom, pos, ref, alt);',
                            db1, 
                            gtxwhat(analysis = analysis1),
                            gtxwhere(chrom = chrom, pos_ge = pos_start, pos_le = pos_end),
                            gtxwhere(chrom = chrom, entity = eq), 
                            gtxanalysisdb(analysis2),
                            gtxwhat(analysis1 = analysis2), # analysis1= argument allows only one analysis
                            gtxwhere(chrom = chrom, pos_ge = pos_start, pos_le = pos_end)
                            ),
                    uniq = FALSE) # expect >=1 rows
  t1 <- as.double(Sys.time())
  gtxlog('Results query returned ', nrow(res), ' rows in ', round(t1 - t0, 3), 's.')
  return(res)
}

#' @export
multicoloc <- function(analysis1, analysis2,
                       chrom, pos_start, pos_end, pos, 
                       hgncid, ensemblid, rs, surround = 0,
                       hard_clip = FALSE, style = 'heatplot',
                       thresh_analysis = 0.1, thresh_entity = 0.1, 
                       dbc = getOption("gtx.dbConnection", NULL), ...) {
  gtxdbcheck(dbc)

  ## get summary stats
  ss <- multicoloc.data(analysis1 = analysis1, analysis2 = analysis2,
                         chrom = chrom, pos_start = pos_start, pos_end = pos_end, pos = pos, 
                         hgncid = hgncid, ensemblid = ensemblid, rs = rs,
                         surround = surround, hard_clip = hard_clip, 
                         dbc = dbc)
  ss <- data.table(ss) # could inline
  
  res <- sqlWrapper(getOption('gtx.dbConnection_cache_genes', dbc), 
                    sprintf('SELECT * FROM genes WHERE %s ORDER BY pos_start;',
                            gtxwhere(ensemblid = unique(ss$entity))), # FIXME not guaranteed entity_type is ENSG
                    uniq = FALSE)
  res$entity <- res$ensemblid # FIXME not guaranteed entity_type
  analyses <- unique(ss$analysis1)

  t0 <- as.double(Sys.time())
  resc1 <- ss[ ,
              subset(coloc.fast(beta1, se1, beta2, se2, ...)$results, hypothesis == 'H4')$posterior,
              by = .(analysis1, entity1)]
  names(resc1)[3] <- 'Hxy'
  resc <- reshape(resc1, direction = 'wide', idvar = 'entity1', timevar = 'analysis1')
  res <- cbind(res, resc[match(res$entity, resc$entity1), ])
  #for (this_analysis in analyses) {
  #    resc <- do.call(rbind,
  #                    lapply(unique(res$entity), function(this_entity) {
  #                        return(subset(coloc.fast(subset(ss, analysis1 == this_analysis & entity1 == this_entity))$results,
  #                                      hypothesis == 'H4')$posterior)
  #                    }))
  #    colnames(resc) <- paste0('Hxy', '_', this_analysis)
  #    res <- cbind(res, resc)
  #}
  t1 <- as.double(Sys.time())
  gtxlog('Colocalization analyzed for ', sum(!is.na(res[ , paste0('Hxy', '.', analyses)])), ' pairs in ', round(t1 - t0, 3), 's')

  if (identical(style, 'none')) {
      ## do nothing
  } else if (identical(style, 'heatplot')) {
      zmat <- as.matrix(res[ , paste0('Hxy', '.', analyses)])
      colnames(zmat) <- analyses
      rownames(zmat) <- with(res, ifelse(hgncid != '', as.character(hgncid), as.character(ensemblid))) # FIXME will this work for all entity types
      ## thresh_analysis <- thresh_analysis*max(zmat, na.rm = TRUE) # threshold could be relative instead of absolute
      ## thresh_entity <- thresh_entity*max(zmat, na.rm = TRUE) # threshold, ditto
      ## FIXME, if nothing passes thresholds, should adapt more gracefully
      zmat.rsel <- apply(zmat, 1, function(x) return(any(x >= thresh_entity, na.rm = TRUE) && any(!is.na(x))))
      if (all(!zmat.rsel)) message('No entities to plot, because no Hxy >= thresh_entity=', thresh_entity)    
      zmat.cord <- order(apply(zmat, 2, function(x) if (any(x >= thresh_analysis, na.rm = TRUE)) max(x, na.rm = TRUE) else NA), na.last = NA)
      if (identical(length(zmat.cord), 0L)) message('No analyses to plot, because no Hxy >= thresh_analysis=', thresh_analysis)
      zmat <- zmat[zmat.rsel, zmat.cord, drop = FALSE]
      if (any(zmat.rsel) && length(zmat.cord) > 0L) {
          multicoloc.plot(zmat, dbc = dbc)
      } else {
          message('Skipping plotting')
      }
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

#' @export
multicoloc.kbs <- function(analysis1, analysis2,
                       chrom, pos_start, pos_end, pos, 
                       hgncid, ensemblid, rs, surround = 0,
                       hard_clip = FALSE, style = 'heatplot',
                       thresh_analysis = 0.1, thresh_entity = 0.1, 
                       dbc = getOption("gtx.dbConnection", NULL), ...) {
  gtxdbcheck(dbc)

  ## get summary stats
  ss <- multicoloc.data(analysis1 = analysis1, analysis2 = analysis2,
                         chrom = chrom, pos_start = pos_start, pos_end = pos_end, pos = pos, 
                         hgncid = hgncid, ensemblid = ensemblid, rs = rs,
                         surround = surround, hard_clip = hard_clip, 
                         dbc = dbc) %>% 
          as.data.table()
  
  res <- sqlWrapper(getOption('gtx.dbConnection_cache_genes', dbc), 
                    sprintf('SELECT * FROM genes WHERE %s ORDER BY pos_start;',
                            gtxwhere(ensemblid = unique(ss$entity))), # FIXME not guaranteed entity_type is ENSG
                    uniq = FALSE)
  res$entity <- res$ensemblid # FIXME not guaranteed entity_type

  ret <- ss[ ,
              coloc.fast(beta1, se1, beta2, se2, ...)$results,
              by = .(analysis1, entity1)] %>% 
         dplyr::select(-label, -prior, -bf) %>% 
         left_join(., 
                  res %>% 
                    rename(entity1 = entity) %>% 
                    select(entity1, hgncid), 
                  by = "entity1") %>% 
        dplyr::select(-posterior, -hypothesis, everything()) %>%
        rename(tissue = analysis1) %>%
        rename(ensembl_id = entity1) %>%
        rename(hgnc_id = hgncid)
  
  return(ret)
}
