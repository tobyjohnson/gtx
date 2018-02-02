## functions for reporting and visualizing a complete gwas

## pval_significance determines the threshold used to declare significance
## pval_plot determines the threshold for true plotting

gwas <- function(analysis,
                 style = 'manhattan',
                 pval_thresh = 5e-08, maf_ge, rsq_ge,
                 dbc = getOption("gtx.dbConnection", NULL)) {
    gtxdbcheck(dbc)
    
    ## FIXME unclear how to handle analysis with entities
    
    t0 <- as.double(Sys.time())
    res <- sqlWrapper(dbc,
                      sprintf('SELECT chrom, pos, ref, alt, pval, rsq, freq, beta, se
                               FROM %s.gwas_results
                               WHERE %s AND %s;', 
                              gtxanalysisdb(analysis),
                              gtxwhat(analysis1 = analysis),
                              gtxfilter(pval_le = pval_thresh, maf_ge = maf_ge, rsq_ge = rsq_ge,
                                        analysis = analysis,
                                        dbc = dbc)),
                      uniq = FALSE)
    res <- data.table(res) # in future, sqlWrapper will return data.table objects always
    t1 <- as.double(Sys.time())
    gtxlog('Results query returned ', nrow(res), ' rows in ', round(t1 - t0, 3), 's.')
    
    t0 <- as.double(Sys.time())
    oo <- res[ , order(suppressWarnings(as.integer(chrom)), pos)]
    res <- res[oo , ]
    t1 <- as.double(Sys.time())
    gtxlog('Ordered by chrom/pos in ', round(t1 - t0, 3), 's.')
    
    ## Simple distance based pruning
    if (nrow(res) > 0) {
        t0 <- as.double(Sys.time())
        ww <- res[ , prune.distance(chrom, pos, pval)]
        res <- res[ww] # ??? test this
        t1 <- as.double(Sys.time())
        gtxlog('Distance based pruning reduced to ', nrow(res), ' rows in ', round(t1 - t0, 3), 's. [PLEASE OPTIMIZE ME!]')
        ## important to prune before using the (currently) expensive gene.annotate() function
        t0 <- as.double(Sys.time())
        res[ , gene_annotation := gene.annotate(chrom, pos)]
        t1 <- as.double(Sys.time())
        gtxlog('Gene annotation added in ', round(t1 - t0, 3), 's. [PLEASE OPTIMIZE ME!]')
        ## IN A FUTURE WORK we will introspect the existence of joint/conditional results
        ## and use if appropriate
        ## or support (via a Hail wrapper) on the fly LD-based pruning
        ## or a HUGE table of pairwise LD values (which could also be used for on-the-fly conditioning on 1 signal)
        
    }
    return(res)

    ## works in progress
    if ('manhattan' %in% style) { # nice semantics
        ## get min and max positions by chromosome, should this be a constant instead of lookup every time?
        mmpos <- sqlWrapper(dbc, 
                            sprintf('SELECT chrom, min(pos) AS minpos, max(pos) AS maxpos
                               FROM %s.gwas_results
                               WHERE %s;'
                              gtxanalysisdb(analysis), 
                              gtxwhat(analysis1 = analysis)),
                            uniq = FALSE)
        ## copy the python code to compute per-chromosome offsets

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
        ## Query total number of associations with p>thresh and passing other filters
        ##    in order to compute expected pvalues for the contents of res (need NOT TO PRUNE above...)
        ## Would we add a second set of points in grey for the filtered-out variants?
        ## -- and a nice legend using expressions for MAF>= & Rsq>=
    }
    
    if ('regionplots' %in% style) {
        ## loop over each row of res and draw a regionplot
        ## OR, do this if regionplot() called with no positional specifier
    }
    
}

