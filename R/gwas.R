## functions for reporting and visualizing a complete gwas

## pval_significance determines the threshold used to declare significance
## pval_plot determines the threshold for true plotting

#' 
#' 
#' @param analysis
#' @param style
#' @param pval_thres
#' @param maf_ge
#' @param rsq_ge
#' @param emac_ge
#' @param case_emac_ge
#' @param gene_annotate
#' @param plot_ymax
#' @param prune_dist
#' @param manhattan_thresh
#' @param manhattan_interspace
#' @param qqplot_col
#' @param qqplot_alpha
#' @param plot_fastbox
#' @param zrocx
#' @param entity
#' @param dbc
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
		 entity = FALSE,
                 dbc = getOption("gtx.dbConnection", NULL)) {
    gtxdbcheck(dbc)
    
    ## FIXME unclear how to handle analysis with entities
    ## FIXME should throw error if entity_type is not NULL
    
    t0 <- as.double(Sys.time())
    res <- sqlWrapper(dbc,
                      sprintf('SELECT chrom, pos, ref, alt, pval, rsq, freq, beta, se
                               FROM %sgwas_results
                               WHERE %s AND %s AND pval IS NOT NULL ORDER BY chrom, pos;', 
                              gtxanalysisdb(analysis),
                              gtxwhat(analysis1 = analysis),
                              gtxfilter(pval_le = pval_thresh,
                                        maf_ge = maf_ge, rsq_ge = rsq_ge,
                                        emac_ge = emac_ge, case_emac_ge = case_emac_ge, 
                                        analysis = analysis,
                                        dbc = dbc)),
                      uniq = FALSE,
                      zrok = zrok)
    res <- data.table::data.table(res) # in future, sqlWrapper will return data.table objects always
    t1 <- as.double(Sys.time())
    gtx_info('Significant results query returned {nrow(res)} rows in {round(t1 - t0, 3)}s.')
    
    if(nrow(res) == 0){
      res <- data.table::data.table(signal     = NA, chrom        = NA, pos_start  = NA, 
                        pos_end    = NA, num_variants = NA, min_pval   = NA, 
                        pos_index  = NA, ref_index    = NA, alt_index  = NA, 
                        pval_index = NA, rsq_index    = NA, freq_index = NA, 
                        beta_index = NA, se_index     = NA)
    } else if (nrow(res) > 0) { 
        ## Simple distance based pruning
        t0 <- as.double(Sys.time())
        res_sigs <- prune.distance(res, surround = prune_dist, sorted = TRUE)
        ## sort by match(chrom, c(as.character(1:22), 'X'))
        rescols <- c('pos', 'ref', 'alt', 'pval', 'rsq', 'freq', 'beta', 'se') # columns from original res to include
        res <- data.table::data.table(cbind(res_sigs, res[res_sigs$row, rescols, with = FALSE]))
        setnames(res, rescols, paste0(rescols, '_index'))
        res[ , row := NULL]
        t1 <- as.double(Sys.time())
        gtx_info('Pruned to {nrow(res)} separate signals in {round(t1 - t0, 3)}s.')
        if (gene_annotate) {
            t0 <- as.double(Sys.time())
            res[ , gene_annotation := gene.annotate(chrom, pos_index)]
            t1 <- as.double(Sys.time())
            gtx_info('Gene annotation added in {round(t1 - t0, 3)}s.')
            if (is.null(getOption('gtx.dbConnection_cache_genes', NULL)) && (t1 - t0) > 1.) {
                gtx_warn('Use gtxcache() to speed up gene annotation')
            }
        }
        ## IN A FUTURE WORK we will introspect the existence of joint/conditional results
        ## and use if appropriate
        ## or support (via a Hail wrapper) on the fly LD-based pruning
        ## or a HUGE table of pairwise LD values (which could also be used for on-the-fly conditioning on 1 signal)
    }
    

    ## Plot description
    ## no handling of entity is required in title
    main <- gtxanalysis_label(analysis = analysis, entity = NULL, nlabel = TRUE, dbc = dbc)
    ## Filtering description
    ## in future we may need to pass maf_lt and rsq_lt as well  
    fdesc <- gtxfilter_label(maf_ge = maf_ge, rsq_ge = rsq_ge,
                             emac_ge = emac_ge, case_emac_ge = case_emac_ge, 
                             analysis = analysis)
    
    if ('manhattan' %in% style || 'qqplot' %in% style) {
        ## Make single query for either case
        t0 <- as.double(Sys.time())
        pvals <- sqlWrapper(dbc,
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
    }

    if ('manhattan' %in% style) { # nice semantics
        ## get min and max positions by chromosome, should this be a constant instead of lookup every time?
        t0 <- as.double(Sys.time())
        mmpos <- sqlWrapper(dbc, 
                            sprintf('SELECT chrom, min(pos) AS minpos, max(pos) AS maxpos
                               FROM %sgwas_results
                               WHERE %s GROUP BY chrom;', 
                               gtxanalysisdb(analysis), 
                               gtxwhat(analysis1 = analysis)),
                            uniq = FALSE)
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

    if ('qqplot' %in% style) {
        t0 <- as.double(Sys.time())
        nump <- as.integer(sqlWrapper(getOption('gtx.dbConnection'),
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
                               uniq = TRUE)$nump) + nrow(pvals)
        pe <- (rank(pvals$pval) - 0.5)/nump # expected p-values
        t1 <- as.double(Sys.time())
        gtx_info('Counted truncated P-values and computed expected P-values in {round(t1 - t0, 3)}s.')

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
    
    if ('regionplots' %in% style) {
        ## loop over each row of res and draw a regionplot
        ## OR, do this if regionplot() called with no positional specifier
    }

    return(as.data.frame(res)) # don't return as a data.table   
}

