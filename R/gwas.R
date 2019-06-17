#define s4 class for gwas
setClass('gwas',
         representation(unPrunedData = 'data.table',
                        prunedData = 'data.table',
                        commonPlotArgs = 'list',
                        manhattanPlotArgs = 'list',
                        qqplotPlotArgs = 'list'),
         prototype(unPrunedData = NULL,
                   prunedData = NULL,
                   commonPlotArgs = NULL,
                   manhattanPlotArgs = NULL,
                   qqplotPlotArgs = NULL))

validateGwasInputs <- function(connectionType){
  checkType(connectionType)
}

#get the sql query string
#NOTE: We might be able to generalize this for all the calls
getSQLQuery <- function(analysis, pval_thresh, maf_ge, rsq_ge, emac_ge, case_emac_ge, dbc){
  queryText <- sprintf('SELECT chrom, pos, ref, alt, pval, rsq, freq, beta, se
                        FROM %sgwas_results
                        WHERE %s AND %s AND pval IS NOT NULL ORDER BY chrom, pos;', 
          gtxanalysisdb(analysis),
          gtxwhat(analysis1 = analysis),
          gtxfilter(pval_le = pval_thresh,
                    maf_ge = maf_ge, rsq_ge = rsq_ge,
                    emac_ge = emac_ge, case_emac_ge = case_emac_ge, 
                    analysis = analysis,
                    dbc = dbc))
  
  return(queryText)
}


#' Get Unpruned Data from the database
#' 
#' This functione is part of the gwas plotting functionality.
#' It retrives data that will ater be pruned
#' 
#' @inheritParams gwas
#' 
#' @return A gwas S4 object
#' 
#' @export
getUnprunedData <- function(connectionType = 'SQL', analysis, pval_thresh,
                            maf_ge, rsq_ge, emac_ge, case_emac_ge, dbc){

  t0 <- as.double(Sys.time())
  
  validateGwasInputsMsg <- validateGwasInputs(connectionType)
  
  callArgs <- switch(connectionType,
                     'SQL' = getSQLQuery(analysis, pval_thresh, maf_ge, rsq_ge,
                                         emac_ge, case_emac_ge, dbc),
                     stop('Unkown connectionType'))
  
  unprunedData <- getDataFromDB(connectionType = connectionType,
                                connectionArguments = list(dbc,
                                                           callArgs,
                                                           uniq = FALSE,
                                                           zrok = zrok))
    
  unprunedData <- data.table::data.table(unprunedData) # in future, getDataFromDB will return data.table objects always
  t1 <- as.double(Sys.time())
  gtx_info('Significant results query returned {nrow(res)} rows in {round(t1 - t0, 3)}s.')
  
  gwasObj <- new('gwas')
  gwasObj@unPrunedData <- unPrunedData
    
  return(gwasObj)
    
}

getEmptyUnprunedData <- function(){
  emptyTabel <- data.table::data.table(signal = NA, chrom = NA, pos_start = NA,
                                       pos_end = NA, num_variants = NA,
                                       min_pval = NA, pos_index = NA,
                                       ref_index = NA, alt_index  = NA,
                                       pval_index = NA, rsq_index = NA,
                                       freq_index = NA, beta_index = NA,
                                       se_index = NA)
  
  return(emptyTable)
  
}


checkPruneDataParams <- function(unPrunedData, gene_annotate){
  validColumnNames <- c("signal", "chrom", "pos_start", "pos_end",
                        "num_variants", "min_pval", "pos_index", "ref_index",
                        "alt_index", "pval_index", "rsq_index", "freq_index",
                        "beta_index", "se_index")
  if (class(unPrunedData)[1] != 'data.table'){
    stop('"unPrunedData" must be a data.table')
  }
  if (!all(validColumnNames %in% colnames (unPrunedData))){
    stop(paste0('unPrunedData must contain the following columns: ',
                paste(validColumnNames, collapse = ", ")))
  }
  if (nrow(unPrunedData) == 0){
    warning('"unPrunedData" has 0 rows.')
  }
  
  if (!is.logical(gene_annotate)){
    stop('gene_annotate must be logical.')
  }
}


#' Function to prune data
#' 
#' @inheritParams gwas
#' @param gwasObj (gwas): a gwas S4 object
#' 
#' @return The input gwas object with added pruned data
#' @export
pruneData <- function(gwasObj, prune_dist, gene_annotate){
  inputData <- gwasObj@unPrunedData
  checkPruneDataParams(inputData, gene_annotate)
  
  if(nrow(inputData) == 0){
    prunedData <- getEmptyUnprunedData()
  } else if (nrow(inputData) > 0) { 
    ## Simple distance based pruning
    t0 <- as.double(Sys.time())
    res_sigs <- prune.distance(prunedData, surround = prune_dist, sorted = TRUE)
    ## sort by match(chrom, c(as.character(1:22), 'X'))
    # columns from original res to include
    rescols <- c('pos', 'ref', 'alt', 'pval', 'rsq', 'freq', 'beta', 'se') 
    prunedData <- data.table::data.table(cbind(res_sigs,
                                               prunedData[res_sigs$row, rescols,
                                                          with = FALSE]))
    setnames(res, rescols, paste0(rescols, '_index'))
    prunedData[ , row := NULL]
    t1 <- as.double(Sys.time())
    gtx_info(paste0('Pruned to {nrow(gwas@unPrunedData)} separate signals ',
                    'in {round(t1 - t0, 3)}s.'))
    
    if (gene_annotate) {
      t0 <- as.double(Sys.time())
      prunedData[ , gene_annotation := gene.annotate(chrom, pos_index)]
      t1 <- as.double(Sys.time())
      gtx_info('Gene annotation added in {round(t1 - t0, 3)}s.')
      if (is.null(getOption('gtx.dbConnection_cache_genes', NULL)) && (t1 - t0) > 1.) {
        gtx_warn('Use gtxcache() to speed up gene annotation')
      }
    }
    ## IN A FUTURE WORK we will introspect the existence of joint/conditional results
    ## and use if appropriate
    ## or support (via a Hail wrapper) on the fly LD-based pruning
    ## or a HUGE table of pairwise LD values (which could also be used for
    ## on-the-fly conditioning on 1 signal)
  }
  
  gwasObj@prunedData <- prunedData
  
}

getCommonGwasParams <- function(obj, dbc, analysis, maf_ge, rsq_ge, emac_ge,
                                case_emac_ge, plot_fastbox){
  ## Plot description
  ## no handling of entity is required in title
  main <- gtxanalysis_label(analysis = analysis, entity = NULL,
                            nlabel = TRUE, dbc = dbc)
  ## Filtering description
  ## in future we may need to pass maf_lt and rsq_lt as well  
  fdesc <- gtxfilter_label(maf_ge = maf_ge, rsq_ge = rsq_ge,
                           emac_ge = emac_ge, case_emac_ge = case_emac_ge, 
                           analysis = analysis)
  
  ## Make single query for either case
  t0 <- as.double(Sys.time())
  pvals <- getDataFromDB(connectionType = 'SQL',
                         connectionArguments = list(dbc,
                                                    sprintf('SELECT chrom, pos, pval
                                                            FROM %sgwas_results
                                                            WHERE %s AND %s AND pval IS NOT NULL;',
                                                            gtxanalysisdb(analysis),
                                                            gtxwhat(analysis1 = analysis),
                                                            gtxfilter(pval_le = 10^-plot_fastbox,
                                                                      maf_ge = maf_ge, rsq_ge = rsq_ge,
                                                                      emac_ge = emac_ge, case_emac_ge = case_emac_ge, 
                                                                      analysis = analysis,
                                                                      dbc = dbc)),
                                                    uniq = FALSE)
  )
  t1 <- as.double(Sys.time())
  gtx_info('Manhattan/QQplot results query returned {nrow(pvals)} rows in {round(t1 - t0, 3)}s.')
  ymax <- max(10, ceiling(-log10(min(pvals$pval))))
  if (ymax > plot_ymax) { # Hard coded threshold makes sense for control of visual display
    ymax <- plot_ymax
    warning('Truncating P-values at 1e-', ymax)
    truncp <- TRUE
  } else {
    truncp <- FALSE
  }
  
  gwasObj@commonPlotArgs <- list(ymax = ymax,
                                 truncp = truncp,
                                 pvals = pvals,
                                 main = main,
                                 fdesc = fdesc)

}

manhattanValues <- function(gwasObj, analysis, manhattan_interspace,
                            manhattan_col, manhattan_thresh, dbc){
  
  pvals <- gwasObj@commonPlotArgs$pvals
  truncp <- gwasObj@commonPlotArgs$truncp
  
  ## get min and max positions by chromosome, should this be a constant instead of lookup every time?
    t0 <- as.double(Sys.time())
    mmpos <- getDataFromDB(connectionType = 'SQL',
                           connectionArguments = list(dbc, 
                                                      sprintf('SELECT chrom, min(pos) AS minpos, max(pos) AS maxpos
                                                              FROM %sgwas_results
                                                              WHERE %s GROUP BY chrom;', 
                                                              gtxanalysisdb(analysis), 
                                                              gtxwhat(analysis1 = analysis)),
                                                      uniq = FALSE))
    mmpos <- mmpos[order_chrom(mmpos$chrom), ]
    mmpos$offset <- c(0, cumsum(as.double(mmpos$maxpos - mmpos$minpos + manhattan_interspace)))[1:nrow(mmpos)] - mmpos$minpos + manhattan_interspace
    mmpos$midpt <- 0.5*(mmpos$maxpos + mmpos$minpos) + mmpos$offset
    mmpos$col <- rep(manhattan_col, length.out = nrow(mmpos))       
    t1 <- as.double(Sys.time())
    gtx_info('Computed chromosome offsets in {round(t1 - t0, 3)}s.')
    
    t0 <- as.double(Sys.time())
    pvals$plotpos <- mmpos$offset[match(pvals$chrom, mmpos$chrom)] + pvals$pos
    pvals$plotcol <- mmpos$col[match(pvals$chrom, mmpos$chrom)]
    my_xlim <- c(min(mmpos$minpos + mmpos$offset), max(mmpos$maxpos + mmpos$offset)) + c(-1, 1)*manhattan_interspace
    my_ylim <- c(0, ymax)
    
    gwasObj@commonPlotArgs$pvals <- pvals
    gwasObj@manhattanPlotArgs <- list(my_xlim = my_xlim,
                                      my_ylim = my_ylim,
                                      mmpos = mmpos)
    
    return(gwasObj)
}

qqplotValues <- function(gwasObj, analysis, plot_fastbox, maf_ge, rsq_ge,
                         emac_ge, case_emac_ge, dbc){

    pvals <- gwasObj@commonPlotArgs$pvals
    truncp <- gwasObj@commonPlotArgs$truncp
    
    t0 <- as.double(Sys.time())
    nump <- as.integer(
      getDataFromDB(connectionType = 'SQL',
                    connectionArguments = list(getOption('gtx.dbConnection'),
                                               sprintf('SELECT count(1) AS nump
                                                       FROM %sgwas_results
                                                       WHERE %s AND %s AND pval IS NOT NULL;',
                                                       gtxanalysisdb(analysis),
                                                       gtxwhat(analysis1 = analysis),
                                                       gtxfilter(pval_gt = 10^-plot_fastbox,
                                                                 maf_ge = maf_ge, rsq_ge = rsq_ge,
                                                                 emac_ge = emac_ge, case_emac_ge = case_emac_ge,
                                                                 analysis = analysis,
                                                                 dbc = dbc)),
                                               uniq = TRUE))$nump) + nrow(pvals)
    pe <- (rank(pvals$pval) - 0.5)/nump # expected p-values
    t1 <- as.double(Sys.time())
    gtx_info('Counted truncated P-values and computed expected P-values in {round(t1 - t0, 3)}s.')
    
    gwasObj@qqplotPlotArgs <- list(pe = pe,
                                   nump = nump)
    
    return(gwasObj)
  
}

#' Prepare data from plotting manhattan and/or qqplot
#' 
#' @inheritParams gwas
#' @inheritParams pruneData
#' 
#' @return A gwas S4 class object with information that will be used by the plotting function
#' 
#' @export
prepareGwasPlot <- function(gwasObj, dbc, analysis, maf_ge, rsq_ge, emac_ge,
                            case_emac_ge, plot_fastbox, style){
  commonGwasParams <- getCommonGwasParams(gwasObj, dbc, analysis, maf_ge,
                                          rsq_ge, emac_ge, case_emac_ge,
                                          plot_fastbox)
  
  if ('manhattan' %in% style){
    gwasObj <- manhattanValues(gwasObj, analysis, manhattan_interspace,
                               manhattan_col, manhattan_thresh, dbc)
  }
  if ('qqplot' %in% style){
    gwasObj <- qqplotValues(gwasObj, analysis, plot_fastbox, maf_ge, rsq_ge,
                            emac_ge, case_emac_ge, dbc)
  }
  
  return(gwasObj)
         
  
}


plotManhattan <- function(gwasObj, manhattan_thresh, plot_fastbox){
  
  mmpos <- gwasObj@manhattanPlotArgs$mmpos
  my_xlim <- gwasObj@manhattanPlotArgs$my_xlim
  my_ylim <- gwasObj@manhattanPlotArgs$my_ylim
  pvals <- gwasObj@commonPlotArgs$pvals
  main <- gwasObj@commonPlotArgs$main
  truncp <- gwasObj@commonPlotArgs$truncp
  fdesc <- gwasObj@commonPlotArgs$fdesc
  
  #plotting
  plot.new()
  plot.window(my_xlim, my_ylim)
  points(pvals$plotpos, pmin(-log10(pvals$pval), ymax), 
         pch = 19, cex = 0.5, col = pvals$plotcol)
  if (plot_fastbox > 0) {
    for (idx in 1:nrow(mmpos)) {
      polygon(c(mmpos$minpos[idx], mmpos$maxpos[idx])[c(1, 2, 2, 1)] + mmpos$offset[idx],
              c(0, 0, plot_fastbox, plot_fastbox),
              density = NA, col = mmpos$col[idx], border = mmpos$col[idx], lwd = 9*.5) # value 9 here is a box fudge to make box line up with pch=19 points
    }
    polygon(my_xlim[c(1, 2, 2, 1)], c(0, 0, plot_fastbox, plot_fastbox), 
            density = 10, angle = 67.5, col = rgb(.75, .75, .75, .5), border = NA)
  }
  abline(h = -log10(manhattan_thresh), col = 'red', lty = 'dashed')
  #axis(1, at = mmpos$midpt, labels = rep(NA, nrow(mmpos)))
  lidx <- rep(1:2, length.out = nrow(mmpos))
  for (idx in 1:2) with(mmpos[lidx == idx, ], mtext(chrom, side = 1, line = c(0, 0.5)[idx], at = midpt, cex = 0.5))
  if (truncp) {
    ## would be nice to more cleanly overwrite y axis label
    axis(2, at = ymax, labels = substitute({}>=ymax, list(ymax = ymax)), las = 1)
    # could e.g. do setdiff(pretty(...), ymax)
  }
  axis(2, las = 1)
  xplt <- par("plt")[2] - par("plt")[1] # figure as fraction of plot, assumes no subsequent changes to par("mar")
  mtext(main, 3, 1,
        cex = min(1., xplt/strwidth(main, units = 'figure')))
  mtext(fdesc, 3, 0, cex = 0.5)
  mtext('Genomic position by chromosome', 1, 2)
  mtext(expression(-log[10](paste(italic(P), "-value"))), 2, 3)
  box()
  t1 <- as.double(Sys.time())
  gtx_info('Manhattan plot rendered in {round(t1 - t0, 3)}s.')
  
  ## aiming for style where non-genome-wide-significant signals are plotted faintly
  ## and perhaps with pseudo-points to avoid expensive overplotting
  ## then vertical lines linking to gene_annotation labels
  ##
  ## set ymax=-log10(min(pval))+0.5
  ## extend a line from each pruned signal vertically to ymax
  ## extend diagonal lines to equally spaced points along the x axis
  ## extend further vertical lines a short distance
  ## ??? hard to know the height on the plot that will be consumed by the labels
  ## unless we can query the aspect ratio and then convert strwidth()
  
  ## if points are drawn using pch = 19, cex = cex
  ## then relativly acceptable 'filling' blocks of solid colour may be drawn using
  ## polygon(x, y, density = NA, col = col, border = col, lwd = 7.5*cex)
  ## relatively acceptable thinning may be achieved with include flag w <- p >= runif(length(p), 0, p^2/.01)
}


plotQqplot <- function(gwasObj, qqplot_alpha, qqplot_col, plot_fastbox){
  
  fdesc <- gwasObj@commonPlotArgs$fdesc
  truncp <- gwasObj@commonPlotArgs$truncp
  ymax <- gwasObj@commonPlotArgs$ymax
  nump <- gwasObj@qqplotPlotArgs$nump
  pe <- gwasObj@qqplotPlotArgs$pe
  pvals <- gwasObj@commonPlotArgs$pvals
  main <- gwasObj@commonPlotArgs$main
  
  ##plotting
  t0 <- as.double(Sys.time())
  qq10.new(pmin = 10^-ymax) # annoying back conversion to p-value scale; FIXME qq10 code should directly support ymax
  qq10.envelope(nump, pmin = 10^-ymax, alpha = qqplot_alpha, col = "grey")
  points(-log10(pe), pmin(-log10(pvals$pval), ymax), 
         pch = 19, cex = 0.5, col = qqplot_col)
  
  if (plot_fastbox > 0) {
    box_x <- -log10(max(pe))
    box_y <- -log10(max(pvals$pval))
    lines(c(0, box_x),
          c(0, box_y),
          lwd = 9*.5, col = qqplot_col) # value 9 here is a box fudge to make box line up with pch=19 points
    polygon(c(0, box_x)[c(1, 2, 2, 1)], c(0, box_y)[c(1, 1, 2, 2)], 
            density = 10, angle = 67.5, col = rgb(.75, .75, .75, .5), border = NA)
  }
  if (truncp) {
    ## would be nice to more cleanly overwrite y axis label
    axis(2, at = ymax, labels = substitute({}>=ymax, list(ymax = ymax)), las = 1)
    # could e.g. do setdiff(pretty(...), ymax)
  }
  xplt <- par("plt")[2] - par("plt")[1] # figure as fraction of plot, assumes no subsequent changes to par("mar")
  mtext(main, 3, 1,
        cex = min(1., xplt/strwidth(main, units = 'figure')))
  mtext(fdesc, 3, 0, cex = 0.5)
  box()
  t1 <- as.double(Sys.time())
  gtx_info('QQ plot rendered in {round(t1 - t0, 3)}s.')
  
  
  ## Query total number of associations with p>thresh and passing other filters
  ##    in order to compute expected pvalues for the contents of res (need NOT TO PRUNE above...)
  ## Would we add a second set of points in grey for the filtered-out variants?
  ## -- and a nice legend using expressions for MAF>= & Rsq>=
}


#' Plotting function for gwas
#' @inheritParams gwas
#' @inheritParams pruneData
#' @export
setMethod("plot", signature(x = 'gwas', y = 'missing'),
          function(x, y, style, qqplot_alpha = NULL, qqplot_col = NULL,
                   plot_fastbox = NULL){
            if ('manhattan' %in% style){
              plotManhattan(x, manhattan_thresh, plot_fastbox)
            }
            if ('qqplot' %in% style){
              plotQqplot(x, qqplot_alpha, qqplot_col, plot_fastbox)
            }
          })

## pval_significance determines the threshold used to declare significance
## pval_plot determines the threshold for true plotting


#' Functions for reporting and visualizing a complete gwas
#' 
#' @param analysis: The key value for the GWAS analysis to summarize
#' @param style: Character vector specifying plot style(s). Can be 'manhattan'
#' or 'qqplot'. style = 'none' will produce no plot.
#' @param pval_thresh: P-value threshold for top hits table
#' @param maf_ge: MAF greater-or-equal threshold
#' @param rsq_ge: Imputation r-squared greater-or-equal threshold
#' @param emac_ge: UNDEFINED
#' @param case_emac_ge: UNDEFINED
#' @param gene_annotate (logical): whether to annotate top hits with gene names
#' @param plot_ymax: Y axis maximum, -log10(p) scale, for plots
#' @param prune_dist: UNDEFINED
#' @param manhattan_thresh: P-value for threshold line on Manhattan plot
#' @param manhattan_col: Chromosome colour cycle for Manhattan plot
#' @param manhattan_interspace: Chromosome interspacing for Manhattan plot
#' @param qqplot_col: Colour for QQ plot
#' @param qqplot_alpha: Null envelope alpha value for QQ plot
#' @param plot_fastbox: Box height for fast approximate plotting
#' @param zrok: UNDEFINED
#' @param dbc: Database connection
#' @param connectionType (character): The type of connection. Only SQL allowed for now
#' 
#' @description This high level function summarizes and visualizes results
#' (i.e. complete summary statistics) from a single GWAS analysis.
#' Currently, a distance-pruned table of top hits, and optional
#' Manhattan and QQ plots, are generated.
#' The distance pruning and gene annotation is currently very slow and
#' needs to be optimized.
#' The Manhattan and QQ plots use an approximation to speed up both the data
#' query and the plot rendering, where boxes or lines are drawn instead of
#' individual points, for y-axis values between 0 and
#' \code{plot_fastbox}.  Faint grey hatching is shown over this part of
#' the plot to indicate that the true data are not being plotted.
#' Y-axis points are truncated, but this needs to be made more clear in the
#' plot.
#' 
#' @return Data frame of pruned data and pplot of style 'manhattan' and/or 'qqplot'
#' 
#' @import data.table
#' @export
gwas <- function(analysis,
                 style = c('manhattan', 'qqplot'), 
                 pval_thresh = 5e-08, maf_ge, rsq_ge, emac_ge, case_emac_ge, 
                 gene_annotate = TRUE,
                 plot_ymax = 30,
                 prune_dist = 500000L,
                 manhattan_thresh = 5e-08,
                 manhattan_col = c('#064F7C', '#6D97BD'),
                 manhattan_interspace = 50e6,
                 qqplot_col = '#064F7C',
                 qqplot_alpha = 0.01,
                 plot_fastbox = 2, 
                 zrok = FALSE,
                 dbc = getOption("gtx.dbConnection", NULL),
                 connectionType = 'SQL') {
    gtxdbcheck(dbc)
    
    ## FIXME unclear how to handle analysis with entities
    ## FIXME should throw error if entity_type is not NULL
    
  
    #get unpruned data
    gwasObj <- getUnprunedData(connectionType = connectionType, analysis, pval_thresh,
                               maf_ge, rsq_ge, emac_ge, case_emac_ge, dbc)
    
    #prune the data
    gwasObj <- pruneData(gwasObj, prune_dist, gene_annotate)
    
    if(style == 'none'){
      gwasObj <- prepareGwasPlot(gwasObj, dbc, analysis, maf_ge, rsq_ge, emac_ge,
                                 case_emac_ge, plot_fastbox, style)
      
      plot(gwasObj, style, qqplot_alpha, qqplot_col, plot_fastbox)
    }

    return(as.data.frame(res)) # don't return as a data.table   
}

