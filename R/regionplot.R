## regionplot.new and associated functions
## assumes database connection is provided by getOption("gtx.dbConnection")

regionplot <- function(analysis, # what analysis (entity should be next)
                       chrom, pos_start, pos_end, pos,  
                       hgncid, ensemblid, rs, surround = 500000, # where
                       entity, 
                       maf_ge, rsq_ge, emac_ge, case_emac_ge, # which variants to include
                       priorsd = 1, priorc = 1e-5, cs_size = 0.95, # parameters for finemapping
                       plot_ymax = 300, # plot output control starts here... # 300 is too high FIXME
                       style = 'signals',
                       protein_coding_only = TRUE, # whether to only show protein coding genes in annotation below plot
                       highlight_style = 'circle', 
                       dbc = getOption("gtx.dbConnection", NULL)) {

  if (isTRUE(getOption('gtx.debug'))) {
	flog.threshold(DEBUG)
  }

  # check database connection
  gtxdbcheck(dbc)

  ## Obtain pvals by passing arguments through to regionplot.data()
  pvals <- regionplot.data(analysis = analysis,
                           chrom = chrom, pos_start = pos_start, pos_end = pos_end, pos = pos, 
                           hgncid = hgncid, ensemblid = ensemblid, rs = rs,
                           surround = surround, 
                           entity = entity,
                           style = style,
                           priorsd = priorsd, priorc = priorc, cs_size = cs_size, 
                           maf_ge = maf_ge, rsq_ge = rsq_ge,
                           emac_ge = emac_ge, case_emac_ge = case_emac_ge, 
                           dbc = dbc)

  ## Determine x-axis range from attr()ibutes of data
  ## was re-resolved from command line args
  xregion <- attr(pvals, 'region')
  #chrom <- xregion$chrom
  #pos_start = xregion$pos_start
  #pos_end = xregion$pos_end
  
  ## Determine entity from attr()ibutes of data
  xentity <- attr(pvals, 'entity')

  pvals <- within(pvals, impact[impact == ''] <- NA)
  pmin <- max(min(pvals$pval), 10^-plot_ymax)

  main <- gtxanalysis_label(analysis = analysis, entity = xentity, nlabel = TRUE, dbc = dbc)
  ## in future we may need to pass maf_lt and rsq_lt as well  
  fdesc <- gtxfilter_label(maf_ge = maf_ge, rsq_ge = rsq_ge, emac_ge = emac_ge, case_emac_ge = case_emac_ge, analysis = analysis)

  ## Highlight index SNP(s) if pos or rs argument was used
  gp <- data.frame(NULL)
  if (!missing(pos)) {
      ## funny syntax to avoid subset(pvals, pos=pos)
      gp <- rbind(gp, pvals[pvals$pos %in% pos, c('pos', 'pval')])
  }
  if (!missing(rs)) {
      # query this rs id
      qrs <- sqlWrapper(dbc,
                       sprintf('SELECT chrom, pos, ref, alt FROM sites WHERE %s;',
                               gtxwhere(chrom = xregion$chrom, rs = rs)),
                        uniq = FALSE, zrok = TRUE)
      gp <- rbind(gp, merge(pvals, qrs, all.x = FALSE, all.y = TRUE)[ , c('pos', 'pval')])
  }

  if ('classic' %in% style) {
    regionplot.new(chrom = xregion$chrom, pos_start = xregion$pos_start, pos_end = xregion$pos_end,
                 pmin = pmin, 
                 main = main, fdesc = fdesc, 
                 protein_coding_only = protein_coding_only, 
                 dbc = dbc)
    ## best order for plotting
    ## Plot all variants with VEP annotation as blue diamonds in top layer
    pvals <- pvals[order(!is.na(pvals$impact), -log10(pvals$pval)), ]
    with(pvals, regionplot.points(pos, pval,
                                  pch = ifelse(!is.na(impact), 23, 21),
                                  col = ifelse(!is.na(impact), rgb(0, 0, 1, .75), rgb(.33, .33, .33, .5)),
                                  bg = ifelse(!is.na(impact), rgb(.5, .5, 1, .75), rgb(.67, .67, .67, .5))))
    regionplot.highlight(gp, highlight_style = highlight_style)
  }
  if ('ld' %in% style) {
    regionplot.new(chrom = xregion$chrom, pos_start = xregion$pos_start, pos_end = xregion$pos_end,
                 pmin = pmin, 
                 main = main, fdesc = fdesc, 
                 protein_coding_only = protein_coding_only, 
                 dbc = dbc)
    ## best order for plotting
    ## Plot variants shaded by LD
    pvals <- pvals[order(!is.na(pvals$impact), -log10(pvals$pval)), ]
    ## current db format means r<0.01 is not included in table, hence substitute missing with zero
    pvals$r[is.na(pvals$r)] <- 0.
    pvals$alpha <- .5 + .25*pvals$r^2 # currently use a single alpha for col and bg for coding and noncoding, so precalculate
    with(pvals, regionplot.points(pos, pval,
                                  pch = ifelse(!is.na(impact), 23, 21),
                                  col = ifelse(!is.na(impact), rgb(0, 0, 0, alpha), rgb(.33, .33, .33, alpha)),
                                  bg = ifelse(!is.na(impact), 
                                              rgb((1 - r^2)*.67, (1 - r^2)*.67, .33*r^2 + .67, alpha),
                                              rgb(.33*r^2 + .67, (1 - r^2)*.67, (1 - r^2)*.67, alpha))))
    ## (re)draw the index variant for pairwise LD over the top, in larger size and with no transparency
    with(merge(attr(pvals, "index_ld"), pvals), 
         regionplot.points(pos, pval, 
                           cex = 1.25, 
                           pch = ifelse(!is.na(impact), 23, 21),
                           col = ifelse(!is.na(impact), rgb(0, 0, 0, alpha = 1), rgb(.33, .33, .33, alpha = 1)),
                           bg = ifelse(!is.na(impact), 
                                       rgb(0, 0, 1, alpha = 1),
                                       rgb(1, 0, 0, alpha = 1))))
    regionplot.highlight(gp, highlight_style = highlight_style)
  }
  if ('signal' %in% style) {
    regionplot.new(chrom = xregion$chrom, pos_start = xregion$pos_start, pos_end = xregion$pos_end,
                 pmin = pmin, 
                 main = main, fdesc = fdesc, 
                 protein_coding_only = protein_coding_only, 
                 dbc = dbc)
    ## best order for plotting is highest pp on top, with impact on top of that
    pvals <- pvals[order(!is.na(pvals$impact), pvals$pp_signal), ]
    pvals <- within(pvals, {
        pprel <- pp_signal/max(pp_signal, na.rm = TRUE) # relative posterior prob for colouring
        alpha <- .5 + .25*pprel # currently use a single alpha for col and bg for coding and noncoding, so precalculate
    })
    with(pvals, regionplot.points(pos, pval,
                                  pch = ifelse(!is.na(impact), 23, 21),
                                  cex = ifelse(cs_signal, 1.25, .75), 
                                  col = ifelse(!is.na(impact), rgb(0, 0, 0, alpha), rgb(.33, .33, .33, alpha)),
                                  bg = ifelse(cs_signal,
                                              ifelse(!is.na(impact),
                                                     rgb((1 - pprel)*.67, (1 - pprel)*.67, .33*pprel + .67, alpha),
                                                     rgb(.33*pprel + .67, (1 - pprel)*.67, (1 - pprel)*.67, alpha)),
                                              rgb(.67, .67, .67, .5))))
    regionplot.highlight(gp, highlight_style = highlight_style)
  }
  if ('signals' %in% style) {
    regionplot.new(chrom = xregion$chrom, pos_start = xregion$pos_start, pos_end = xregion$pos_end,
                 pmin = pmin, 
                 main = main, fdesc = fdesc, 
                 protein_coding_only = protein_coding_only, 
                 dbc = dbc)

    signals <- unique(na.omit(pvals$signal))
    if (length(signals) == 0L) {
        flog.debug('No CLEO results, using single signal results instead')
        ## no signals, so use single signal results as if they were the CLEO results
        pvals <- within(pvals, {
            cs_cleo <- cs_signal
            pp_cleo <- pp_signal
            signal <- ifelse(cs_signal, 'default', NA)
        })
        signals <- 'default'
    }
    ## best order for plotting
    ## Plot variants coloured/sized by credible set, with VEP annotation is diamond shape in top layer
    pvals <- pvals[order(!is.na(pvals$impact), pvals$pp_cleo), ]

    ## colvec is vector of identifying colour for each signal
    if (length(signals) == 1L) {
        colvec <- rgb(1, 0, 0) # red
    } else {
        colvec <- rainbow(length(signals), start = .6667, end = .3333) # blue-magenta-red-yellow-green avoiding cyan part of spectrum
    }

    bg_cleo <- rep(rgb(.67, .67, .67, alpha = .5), nrow(pvals))
    for (idx in 1:length(signals)) {
        ww <- which(pvals$signal == signals[idx])
        if (length(ww) > 0L) {
            pprel <- pvals$pp_cleo[ww]
            pprel <- pprel/max(pprel, na.rm = TRUE) # relative pp *within signal*
            bg_cleo[ww] <- rgb((col2rgb(colvec[idx])[1]/255 - .67)*pprel + .67, 
                               (col2rgb(colvec[idx])[2]/255 - .67)*pprel + .67,
                               (col2rgb(colvec[idx])[3]/255 - .67)*pprel + .67, 
                               alpha = .5)
        }
    }

    with(pvals, regionplot.points(pos, pval,
                                  pch = ifelse(!is.na(impact), 23, 21), 
                                  cex = ifelse(cs_cleo, 1.25, .75),
                                  bg = bg_cleo, 
                                  col = ifelse(!is.na(impact), rgb(0, 0, 0, .5), rgb(.33, .33, .33, .5))))

    ## legend indicating signals
    legend("bottomleft", pch = 21, col = rgb(.33, .33, .33, .5),
           pt.bg = colvec, legend=paste0('#', signals),
           horiz=T, bty="n", cex=.5)
    regionplot.highlight(gp, highlight_style = highlight_style)
  } 

  return(invisible(NULL))
}

regionplot.data <- function(analysis,
                            chrom, pos_start, pos_end, pos, 
                            hgncid, ensemblid, rs, surround = 500000,
                            entity, 
                            style = 'signals',
                            priorsd = 1, priorc = 1e-5, cs_size = 0.95, 
                            maf_ge, rsq_ge,
                            emac_ge, case_emac_ge, 
                            dbc = getOption("gtx.dbConnection", NULL)) {
    ## check database connection
    gtxdbcheck(dbc)

    style <- tolower(style)
    
    ## Determine x-axis range from arguments
    xregion <- gtxregion(chrom = chrom, pos_start = pos_start, pos_end = pos_end, pos = pos, 
                         hgncid = hgncid, ensemblid = ensemblid, rs = rs, surround = surround,
                         dbc = dbc)
    chrom = xregion$chrom # note, overwriting command line arguments
    pos_start = xregion$pos_start
    pos_end = xregion$pos_end
    
    ## If required, determine entity and associated info including entity_label
    xentity <- gtxentity(analysis, entity = entity, hgncid = hgncid, ensemblid = ensemblid)
    
    ## always query marginal p-values
    ## seems more flexible to query CLEO results separately and merge within R code
    flog.debug('Querying marginal p-values')
    t0 <- as.double(Sys.time())
    pvals <- sqlWrapper(dbc,
                        sprintf('SELECT gwas_results.chrom, gwas_results.pos, gwas_results.ref, gwas_results.alt, pval, impact %s FROM %s.gwas_results LEFT JOIN vep USING (chrom, pos, ref, alt) WHERE %s AND %s AND %s AND %s AND pval IS NOT NULL;',
                                if (any(c('signal', 'signals') %in% style)) ', beta, se, rsq, freq' else '', 
                                gtxanalysisdb(analysis), 
                                gtxwhat(analysis1 = analysis),
                                if (!is.null(xentity)) sprintf('feature=\'%s\'', xentity$entity) else '(True)',
                                gtxwhere(chrom, pos_ge = pos_start, pos_le = pos_end, tablename = 'gwas_results'),
                                gtxfilter(maf_ge = maf_ge, rsq_ge = rsq_ge, emac_ge = emac_ge, case_emac_ge = case_emac_ge, analysis = analysis)),
                        uniq = FALSE)
    t1 <- as.double(Sys.time())
    flog.debug(paste0('Query returned ', nrow(pvals), ' variants in query region ', xregion$label, ' in ', round(t1 - t0, 3), 's.'))

    if (any(c('signal', 'signals') %in% style)) {
        flog.debug('Finemapping under single signal assumption')
        ## cs_only = FALSE since we still want to plot/return variants not in the credible set
        pvals <- fm_signal(pvals, priorsd = priorsd, priorc = priorc, cs_size = cs_size, cs_only = FALSE)
    }
    
    if ('signals' %in% style) {
        ## get CLEO credible sets only, since we effectively left join with this, for plot colouring
        fmres <- fm_cleo(analysis = analysis, chrom = chrom, pos_start = pos_start, pos_end = pos_end,
                         priorsd = priorsd, priorc = priorc, cs_size = cs_size, cs_only = TRUE)
        ## note fm_cleo already prints logging messages about number of signals
        if (nrow(fmres) > 0) {
            ## sort by decreasing posterior probability, match each
            ## row or pvals with *first* match in fmres
            ## thus linking any variants in more than one
            ## credible set, with the one for the signal for which
            ## it has higher probability
            fmres <- fmres[order(fmres$pp_cleo, decreasing = TRUE), ]
            pvals <- cbind(pvals,
                           fmres[match(with(pvals, paste(chrom, pos, ref, alt, sep = '_')),
                                       with(fmres, paste(chrom, pos, ref, alt, sep = '_'))),
                                 c('signal', 'cs_cleo', 'pp_cleo')])
            ## FIXME in case one variant is in more than one credible
            ## set, it would be better to cbind like we do above with signal,
            ## but then cbind with the *marginal* pp's from aggregating over signals
            pvals$pp_cleo[is.na(pvals$pp_cleo)] <- 0. # make zero to be safe with sorting/plotting
            pvals$cs_cleo[is.na(pvals$cs_cleo)] <- FALSE # make FALSE to be safe with tables/plotting
        } else {
            # do nothing
        }
    }

    if ('ld' %in% style) {
        # merge with LD information in userspace code since we need to
        # pull the p-values first to find the top hit or otherwise chosen variant

        # note if style='ld', sort order is by ld with index variant,
        # otherwise sort by pval (see else block)

        # UB-133 Mark pvals with LD
        region <- list(
          chrom = pvals$chrom[1],
          range = range(pvals$pos, na.rm=TRUE, finite=TRUE)
        )

        ld_check <- tbl(getOption("gtx.dbConnection"), 'ld') %>%
          filter(chrom1==region$chrom & 
                   pos1 >= region$range[1] & 
                   pos1 <= region$range[2]) %>%
          select(chrom1,pos1,ref1,alt1) %>%
          rename(chrom=chrom1,pos=pos1,ref=ref1,alt=alt1) %>%
          group_by(chrom,pos,ref,alt) %>%
          summarize(has_ld = 1) %>%
          collect() %>%
          merge(pvals, all.x = FALSE, all.y = TRUE)
        
        missing.ld <- ld_check %>% filter(is.na(has_ld)) %>% nrow
        if (missing.ld > 0) {
          flog.warn(paste0(missing.ld,  ' snps did not have LD data and will not',
                           ' be considered for index snp selection.'))
        }
        
        # Determine Index Variant
        pval1 <- NULL # store chrom/pos/ref/alt of the index variant
        if (!missing(pos)) {
          ## funny syntax to avoid subset(pvals, pos %in% pos)
          pval1 <- pvals[pvals$pos %in% pos, c('chrom', 'pos', 'ref', 'alt'), drop = FALSE]
          if (identical(nrow(pval1), 1L)) {
            pval1$r <- 1
            flog.debug(paste0('Querying pairwise LD with index chr', pval1$chrom, ':', pval1$pos, ':', pval1$ref, ':', pval1$alt, 
                   ' (selected by pos argument)'))
          } else {
            flog.debug(paste0('Skipping pos argument [ ', paste(pos, collapse = ', '), 
                   ' ] for pairwise LD index selection because matches ',
                   nrow(pval1), ' variants'))
            pval1 <- NULL
          }
        }
        if (is.null(pval1) && !missing(rs)) {
          # query this rs id
          qrs <- sqlWrapper(dbc,
                       sprintf('SELECT chrom, pos, ref, alt FROM sites WHERE %s;',
                               gtxwhere(chrom = chrom, rs = rs)),
                        uniq = FALSE, zrok = TRUE)
          # use merge as a way to subset on match by all of chrom/pos/ref/alt
          # note default merge is all.x = FALSE, all.y = FALSE
          pval1 <- ld_check %>% filter(has_ld == 1) %>%
            merge(qrs) %>%
            select(chrom, pos, ref, alt)
          
          if (identical(nrow(pval1), 1L)) {
            pval1$r <- 1
            flog.debug(paste0('Querying pairwise LD with index chr', pval1$chrom, ':', pval1$pos, ':', pval1$ref, ':', pval1$alt, 
                   ' (selected by rs argument)'))
          } else {
            flog.debug(paste0('Skipping rs argument [ ', paste(rs, collapse = ', '), 
                   ' ] for pairwise LD index selection because matches ',
                   nrow(pval1), ' variants'))
            pval1 <- NULL
          }
        }
        if (is.null(pval1)) {
          pval1 <- ld_check %>% 
            filter(has_ld == 1) %>% 
            filter(rank(pval, ties.method='first')==1) %>% 
            select(chrom,pos,ref,alt) %>%
            mutate(r=1)
          
          flog.debug(paste0('Querying pairwise LD with index chr', pval1$chrom, ':', pval1$pos, ':', pval1$ref, ':', pval1$alt,
                 ' (selected by smallest pval)'))
        }
        stopifnot(identical(nrow(pval1), 1L))
        ld1 <- sqlWrapper(dbc, 
                          sprintf('SELECT chrom2 AS chrom, pos2 AS pos, ref2 AS ref, alt2 AS alt, r FROM ld
                                   WHERE chrom1=\'%s\' AND pos1=%s AND ref1=\'%s\' AND alt1=\'%s\';',
                                  pval1$chrom, pval1$pos, pval1$ref, pval1$alt), # should we sanitize
                          uniq = FALSE)
        pvals <- merge(pvals, rbind(pval1, ld1), all.x = TRUE, all.y = FALSE)
        # sort by decreasing r^2 (will be resorted for plotting, so makes 
        # a .data() call work as a useful proxy search
        pvals <- pvals[order(pvals$r^2, decreasing = TRUE), ]
        attr(pvals, 'index_ld') <- pval1
    } else {
        ## sort by increasing pval
        pvals <- pvals[order(pvals$pval), ]
    }
    
    attr(pvals, 'analysis') <- analysis
    attr(pvals, 'region') <- xregion
    attr(pvals, 'entity') <- xentity
    return(pvals)
}

regionplot.new <- function(chrom, pos_start, pos_end, pos, 
                           hgncid, ensemblid, rs, surround = 500000, 
                           pmin = 1e-10, main, fdesc, 
			   protein_coding_only = TRUE,   
                           dbc = getOption("gtx.dbConnection", NULL)) {
  gtxdbcheck(dbc)
    
  ## Determine x-axis range from arguments
  xregion <- gtxregion(chrom = chrom, pos_start = pos_start, pos_end = pos_end, pos = pos, 
                       hgncid = hgncid, ensemblid = ensemblid, rs = rs, surround = surround,
                       dbc = dbc)
  chrom = xregion$chrom
  pos_start = xregion$pos_start
  pos_end = xregion$pos_end

  ## Determine y-axis upper limit
  ymax <- ceiling(-log10(pmin) + 0.5)

  ## Determine amount of y-axis space needed for gene annotation
  gl <- regionplot.genelayout(chrom, pos_start, pos_end, ymax, protein_coding_only = protein_coding_only)
  
  ## Set up plotting area
  plot.new()
  par(mar = c(4, 4, 4, 4) + 0.1, xaxs = 'i') # xaxs='i' stops R expanding x axis
  plot.window(c(pos_start, pos_end), range(gl$yline))
  
  ## Draw axes, and axis labels
  abline(h = 0, col = "grey")
  with(list(xpretty = pretty(c(pos_start, pos_end)*1e-6)),
       axis(1, at = xpretty*1e6, labels = xpretty))
  axis(2, at = pretty(c(0, ymax)), las = 1)
  mtext(paste0("chr", chrom, " position (Mb)"), 1, 3)
  mtext(expression(paste("Association ", -log[10](paste(P, "-value")))), 2, 3)

  ## Add recombination rate and gene annotation
  ##  regionplot.recombination determines position range from par("usr")
  regionplot.recombination(chrom, yoff = mean(gl$yline[2:3]))
  ##  regionplot.genedraw uses previously determined layout
  regionplot.genedraw(gl)

  ## Add title, scaled to fit if necessary, CHANGE TO USE mtext.fit() FIXME
  if (!missing(main)) {
    xplt <- par("plt")[2] - par("plt")[1] # figure as fraction of plot, assumes no subsequent changes to par("mar")
    mtext(main, 3, 1,
          cex = min(1., xplt/strwidth(main, units = 'figure')))
  }
  if (!missing(fdesc)) {
      mtext(fdesc, 3, 0, cex = 0.5)
  }

  ## Draw box last to overdraw any edge marks
  box()

  return(invisible(NULL))
}

regionplot.genedraw <- function(gl) {
  arrows(gl$genelayout$pos_start, gl$genelayout$y, x1 = gl$genelayout$pos_end,
         length = 0, lwd = 2, col = "blue")
  text(gl$genelayout$pos_mid, gl$genelayout$y - .25*strheight("M", cex = gl$cex),
       gl$genelayout$label,
       adj = c(.5, 1), font = 3, cex = gl$cex, col = "blue")
  return(invisible(NULL))
}

regionplot.genelayout <- function (chrom, pos_start, pos_end, ymax, cex = 0.75, 
				   protein_coding_only = TRUE, 
                                   dbc = getOption("gtx.dbConnection", NULL)) {
  gtxdbcheck(dbc)
  
  yplt <- par("plt")[4] - par("plt")[3] # figure as fraction of plot, assumes no subsequent changes to par("mar")
  xplt <- par("plt")[2] - par("plt")[1]
  xusr <- pos_end - pos_start # par("usr")[2] - par("usr")[1] 
  return(with(sqlWrapper(getOption('gtx.dbConnection_cache_genes', dbc), 
                       ## use SQL query to aggregate over multiple rows with same name to get whole span
                       sprintf('SELECT min(pos_start) AS pos_start, max(pos_end) AS pos_end, hgncid, ensemblid FROM genes WHERE %s %s GROUP BY ensemblid, hgncid ORDER BY pos_start', 
                               gtxwhere(chrom = chrom, pos_end_ge = pos_start, pos_start_le = pos_end),
			       if (protein_coding_only) 'AND genetype=\'protein_coding\'' else ''),
			uniq = FALSE, zrok = TRUE),
              {
                label <- ifelse(hgncid != '', as.character(hgncid), as.character(ensemblid)) # could also check NULL or NA?
                ## compute start and end plot positions for each gene, using larger of transcript line and gene name annotation
                pos_mid <- .5*(pos_start + pos_end)
                xwidth <- strwidth(label, units = 'figure', cex = cex)*xusr/xplt
                xpad <- max(10000, strwidth('M', units = 'figure', cex = cex)*xusr/xplt)
                xstart <- pmin(pos_start, pos_mid - .5*xwidth) - xpad
                xend <- pmax(pos_end, pos_mid + .5*xwidth) + xpad
                ## layout engine based on putting each gene on the highest line possible
                iline <- rep(NA, length(label))
                if (length(label) > 0) {
                  iline[1] <- 0
                  if (length(label) > 1) {
                    for (idx in 2:length(label)) {
                      iline[idx] <- with(rbind(aggregate(xend[1:(idx - 1)], by = list(tline = iline[1:(idx - 1)]), FUN = max),
                                               data.frame(tline = max(iline[1:(idx - 1)]) + 1, x = -Inf)),
                                         min(tline[x < xstart[idx]]))
                    }
                  }
                }
                fy <- max(iline + 1)*2*strheight("M", units = "figure", cex = cex)/yplt # fraction of total figure region needed
                if (fy > 0.6) {
                  warning('Squashing gene annotation to fit within .6 of plot area')
                  fy <- 0.6
                }
                yd1 = ymax/(0.9-fy) * 0.1
                yd2 = ymax/(0.9-fy) * fy
                yline <- c(-(yd1 + yd2), -yd1, 0, ymax)
                y <- -(iline/max(iline + 1)*yd2 + yd1)
                list(cex = cex, yline = yline, genelayout = data.frame(label, pos_start, pos_end, pos_mid, iline, y))
              }))
}

regionplot.points <- function(pos, pval,
                              pch = 21, bg = rgb(.67, .67, .67, .5), col = rgb(.33, .33, .33, .5), cex = 1, 
                              suppressWarning = FALSE) {
  ymax <- floor(par("usr")[4])
  y <- -log10(pval)
  f <- y > ymax
  points(pos, ifelse(f, ymax, y), pch = pch, col = col, bg = bg, cex = cex)
  if (any(f) && !suppressWarning) {
    # would be nice to more cleanly overwrite y axis label
    axis(2, at = ymax, labels = substitute({}>=ymax, list(ymax = ymax)), las = 1)
    text(mean(range(pos[f])), ymax, 'Plot truncated', adj = c(0.5, 0), cex = 0.5)
    return(FALSE)
  }
  return(TRUE)
}

regionplot.highlight <- function(pvals, highlight_style) {
    stopifnot(is.data.frame(pvals))
    if (all(c('pos', 'pval') %in% names(pvals))) {
        pvals <- subset(pvals, !is.na(pos) & !is.na(pval))
        if (nrow(pvals) > 0) {
            if ('circle' %in% highlight_style) {
                regionplot.points(pvals$pos, pvals$pval, pch = 1, cex = 2, col = rgb(0, 0, 0, .75))
            }
            if ('pos' %in% highlight_style) {
                text(pvals$pos, -log10(pvals$pval), pvals$pos, pos = 3, cex = 0.5)
            }
            ## in future will add options if 'rs' %in% highlight_style etc.
        }
    } else {
        # silently do nothing, e.g. if pvals==data.frame(NULL)
    }
    return(invisible(NULL))
}

regionplot.recombination <- function(chrom, pos_start, pos_end, yoff = -.5, 
                                     dbc = getOption("gtx.dbConnection", NULL)) {
  gtxdbcheck(dbc)

  ## Recombination segments (r) that fall wholly or partly within [pos_start, pos_end]
  ## have r.pos_end >= pos_start and r.pos_start <= pos_end
  if (missing(pos_start)) pos_start = floor(par("usr")[1])
  if (missing(pos_end)) pos_end = ceiling(par("usr")[2])

  with(sqlWrapper(dbc, 
                sprintf('SELECT pos_start, pos_end, recombination_rate FROM genetic_map WHERE %s', 
                        gtxwhere(chrom = chrom, pos_end_ge = pos_start, pos_start_le = pos_end)),
		uniq = FALSE, zrok = TRUE),
       {
         abline(h = yoff, col = "grey")
         yscale <- (par("usr")[4] - yoff)*.9/max(recombination_rate)
         lines(c(pos_start, pos_end[length(pos_end)]),
               yoff + c(recombination_rate, recombination_rate[length(recombination_rate)])*yscale,
               type = "s", col = "cyan3")
         with(list(ypretty = pretty(c(0, max(recombination_rate)))),
              axis(4, at = yoff + ypretty*yscale, labels = ypretty, las = 1))
       })
  mtext("Recombination rate (cM/Mb)", 4, 3)
  return(invisible(NULL))
}

#library for region plot
#Author: Li Li, modified based on Toby's code



#draw recomination rate 
recomb.draw <- function(chr, from.bp, to.bp, ylo, yhi, recomb){
  if (is.integer(chr)) chr <- as.character(chr)
  if(chr == "23") chr <- "X"
  chr <- ifelse(substr(chr, 1, 3) == "chr", chr, paste("chr", chr, sep = ""))
  #recomb <- read.table(paste("genetic_map_chr", chr, ".txt", sep=""), header=T)
  #load recombination info
  recomb.chr <- recomb[[chr]][-1]
  idx.start<- which(recomb.chr[,1] < from.bp)
  idx.end<- which(recomb.chr[,1] > to.bp)
  keep.recomb <-recomb.chr[idx.start[length(idx.start)]: idx.end[1],]
  big.range <- yhi -ylo
  max.curr<- max(keep.recomb[,2], na.rm = T)
  #print (paste("Max recomb rate in the region is", max.curr))
  ticks.curr <- pretty(c(0, max.curr))
  max.curr<- ticks.curr[length(ticks.curr )]
  lines(keep.recomb[,1], ylo + ( ( keep.recomb[,2] / max.curr ) * big.range), type="l", 
        col="lightblue", lwd=1)
  axis(4, at=ylo + big.range /max.curr * ticks.curr, labels=ticks.curr, las=1, hadj = 0)
  mtext("Recomb rate (cM/Mb)", side=4, at=(ylo+big.range/2), line=2.5) 
}



################################################################################
#gwas: gwas results with columns for MARKER, chr, POS (bp), qcflag, genoflag, pvalues
#chr pos: hit snp
#flanking: flanking rang (kb) for plotting
#col.pvals: column names for p values to be plotted by different shapes (21:25, 3,4,7:14) 
#plabel: label for p value column
#col.r2: column name for r2 to index SNP. Background colored by LD to index SNP (white to red)
################################################################################
regionplot.multiP <- function(gwas, chr, pos, flanking = 250, col.r2=NULL,
                              col.pvals = c("P1","P2"), plabel = NULL, gencode = gencode, recomb=recomb,
                              main = "", main.sub = NULL, miny = 8, crity = -log10(5e-08)) {
  to.bp <- min(pos+flanking *1000, max(gwas$POS, na.rm=T))
  from.bp <- max(pos - flanking *1000, min(gwas$POS, na.rm = T))
  stopifnot(to.bp > from.bp)
  stopifnot(all(col.pvals %in% names(gwas)))
  sel <- gwas$CHROM == chr & gwas$POS >= from.bp & gwas$POS <= to.bp
  cat(sum(sel), "SNPs selected\n")
  pos.sel<- gwas$POS[sel]
  minP <- 1e-50 #minimum value for p value
  for( c in col.pvals) gwas[[c]][!is.na(gwas[[c]]) & gwas[[c]] < minP] <- minP
  lp <- -log10(as.matrix(subset(gwas, sel, col.pvals))) # lp = -log10(p)
  if (is.null(crity)) crity <- -log10(0.05/nrow(lp))
  if (is.null(plabel)) plabel <- col.pvals
  r2.sel <- rep(1, sum(sel))
  if(!is.null(col.r2)) {
    r2.sel<- gwas[[col.r2]][sel]
    r2.sel[is.na(r2.sel)]<- 0
    r2.sel <- pmin(1, pmax(0, 1-r2.sel))
  }
 
 
  plot.new()
  yhi<- max(c(miny, 1.05 * max(lp, na.rm = TRUE)))
  y.tick<- pretty(c(0, yhi))
  yhi <- y.tick[length(y.tick)]
  scale  <- yhi/8
  plot.window(c(from.bp, to.bp), c(-6 *scale,  yhi))
  
  abline(h = seq(0, 7, 1), col = "lightgrey", lty = "dotted")
  abline(h = crity, col = "red", lty = "dashed")
  ppch <- c(21:25, 3,4,7:14) #rep(21:25, 5)
  index <- pos.sel == pos 
  for (idx in 1:ncol(lp)) {
    points(pos.sel[-index], lp[-index , idx], pch = ppch[idx], col = "black",  bg = rgb(1, r2.sel[-index], r2.sel[-index]), cex = 1)
    if(any(index)) #index snp
      points(pos, lp[index, idx],  pch = ppch[idx], col = rgb(0, 0, 1),  bg = rgb(0, 0, 1), cex = 1.3)
  }   
  legend("topright",   ncol = 2,   bty = "n",  pch = ppch[1:ncol(lp)], 
           text.col =  "black", cex = 0.8,  legend = plabel)
  for(i in 0:10){
    polygon(from.bp + (c(0, 1, 1, 0)+i) *(to.bp-from.bp)/80, yhi- c(0.5, 0.5, 1, 1)*scale, 
            border = "grey", col = rgb(1, 1-i/10, 1-i/10)) 
  }
  text(from.bp, yhi - 0 *scale , expression(paste("LD Map Type:","r"^"2")), adj = 0, cex = 0.8)
  text(from.bp+seq(0.5, 10.5, 2) * (to.bp-from.bp)/80, rep(yhi, 6)- 1.5*scale, 
       labels =c(0, 0.2, 0.4, 0.6, 0.8,1),  cex = 0.7)
  
  xm <- pretty(c(from.bp, to.bp)*1e-6) # x marks
  axis(1, at = xm*1e6, labels = xm)
  axis(2, at = y.tick, las = 1)
  mtext(expression(-log[10](italic(P))), side=2, at=yhi/2, line=2.5) 
  recomb.draw(chr, from.bp, to.bp, -1*scale, yhi, recomb)                        
  gene.draw(paste("chr", chr, sep = ""), from.bp, to.bp, gencode, 
            yhi = -2*scale, ylo = -6*scale,exony = 0.05*scale, genecex = 0.7)
  mtext(paste("Chromosome ", chr, " genomic position (Mb)", sep = ""), side = 1, line = 2.5)
  title(main = main)  # ylab = expression(-log[10](italic(P))),
  if(!is.null(main.sub))
    title(main = main.sub, cex.main = 1, col.main = "blue", line = 0.5)
  box()
  return(invisible(NULL))
}


