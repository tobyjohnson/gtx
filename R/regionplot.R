## regionplot.new and associated functions
## assumes database connection is provided by getOption("gtx.dbConnection")

regionplot <- function(analysis,
                       chrom, pos_start, pos_end, pos, 
                       hgncid, ensemblid, rs, surround = 500000,
                       entity, 
                       style = 'signals', 
                       protein_coding_only = TRUE, # whether to only show protein coding genes in annotation below plot
                       maf_ge, rsq_ge, 
                       dbc = getOption("gtx.dbConnection", NULL)) {
  # check database connection
  gtxdbcheck(dbc)

  ## Determine x-axis range from arguments
  xregion <- gtxregion(chrom = chrom, pos_start = pos_start, pos_end = pos_end, pos = pos, 
                       hgncid = hgncid, ensemblid = ensemblid, rs = rs, surround = surround,
                       dbc = dbc)
  chrom = xregion$chrom
  pos_start = xregion$pos_start
  pos_end = xregion$pos_end

  ## Determine entity, if required, for each analysis
  xentity <- gtxentity(analysis, entity = entity, hgncid = hgncid, ensemblid = ensemblid)

  ## Obtain pvals
    
  if (identical(style, 'classic')) {
    ## basic query without finemapping
      pvals <- regionplot.data(analysis = analysis,
                               chrom = xregion$chrom, pos_start = xregion$pos_start, pos_end = xregion$pos_end,
                               entity = xentity$entity,
                               style = style,
                               maf_ge = maf_ge, rsq_ge = rsq_ge, 
                               dbc = dbc)$pvals
  } else if (identical(style, 'signals')) {
    ## plot with finemapping annotation

    ## need to think more carefully about how to make sure any conditional analyses are either fully in, or fully out, of the plot
    pvals <- regionplot.data(analysis = analysis,
                             chrom = xregion$chrom, pos_start = xregion$pos_start, pos_end = xregion$pos_end,
                             entity = xentity$entity,
                             style = style,
                             maf_ge = maf_ge, rsq_ge = rsq_ge, 
                             dbc = dbc)$pvals
    signals <- sort(unique(na.omit(pvals$signal)))
    nsignals <- length(signals)

    if (nsignals == 0) {
      message("No signals from conditional/joint analyses overlap this region")
      ## No signals from conditional stepwise analyses at the threshold those analyses were done to
      ## If at least one P<=1e-4, fine map based on assumption of one signal
      ## Otherwise, no fine mapping (the fixed prior, not dependent on allele freq,
      ##   gives unexpected results for very weak association signals)
      if (min(pvals$pval) <= 1e-4) {
        pvals <- within(pvals, {
          posterior <- norm1(abf.Wakefield(beta, se, priorsd = 1, log = TRUE), log = TRUE)
          posterior <- ifelse(!is.na(posterior), posterior, 0) # bad hack
          c95 <- credset(posterior)
          ## Colour with red->grey scale for the one signal
          colour <- rgb((255 - 171)*posterior/max(posterior)+171, 
                        (0 - 171)*posterior/max(posterior)+171, 
                        (0 - 171)*posterior/max(posterior)+171, 
                        alpha = 127, maxColorValue = 255)
        })
        nsignals <- 1 # used for legend below # SHOULD HAVE SOME KIND OF ANNOTATION WHEN THERE ARE NO REAL SIGNALS LIKE THIS
        colvec <- rgb(255, 0, 0, maxColorValue = 255) # used for legend below
      } else {
        pvals <- within(pvals, {
          posterior <- NA # set just so that ordering for plotting does not throw error
	  c95 <- FALSE
          colour <- rgb(171, 171, 171, alpha = 127, maxColorValue = 255)
        })
      }
    } else {
      message(nsignals, " signal(s) from conditional/joint analyses overlap this region") 
      ## one or more signals, same code works for one as >one except for setting up nice colour vector
      if (nsignals == 1) {
        colvec <- rgb(255, 0, 0, maxColorValue = 255)
      } else {
        colvec <- rainbow(nsignals, start = .6667, end = .3333) # blue-magenta-red-yellow-green avoiding cyan part of spectrum
      }
      pvals <- within(pvals, {
        posterior <- rep(NA, length(pval))
        colour <- rep(NA, length(pval))
        c95 <- rep(NA, length(pval))
        for (idx in 1:nsignals) {
          f <- !is.na(signal) & signal == signals[idx]
          if (any(f)) {
            posterior[f] <- norm1(abf.Wakefield(beta_cond[f], se_cond[f], priorsd = 1, log = TRUE), log = TRUE)
            posterior[f] <- ifelse(!is.na(posterior[f]), posterior[f], 0) # bad hack
            c95[f] <- credset(posterior[f])
            colour[f] <- rgb((col2rgb(colvec[idx])[1] - 171)*posterior[f]/max(posterior[f])+171, 
                             (col2rgb(colvec[idx])[2] - 171)*posterior[f]/max(posterior[f])+171, 
                             (col2rgb(colvec[idx])[3] - 171)*posterior[f]/max(posterior[f])+171, 
                             alpha = 127, maxColorValue = 255)
          }
          f <- NULL
      }})
      # ensure one point per variant, using the annotation for the signal with the highest posterior
      pvals <- pvals[order(pvals$posterior, decreasing = TRUE), ]
      pvals <- pvals[!duplicated(with(pvals, paste(chrom, pos, ref, alt, sep = "_"))), ]
    }
  } else {
    stop('style ', style, ' not recognised')
  }

  pvals <- within(pvals, impact[impact == ''] <- NA)
  pmin <- min(pvals$pval)

  pdesc <- sqlWrapper(dbc, 
                      sprintf('SELECT label, ncase, ncontrol, ncohort FROM analyses WHERE %s;',
                              gtxwhat(analysis1 = analysis)))
  main <- if (!is.na(pdesc$ncase[1]) && !is.na(pdesc$ncontrol[1])) {
              sprintf('%s, n=%i vs %i', pdesc$label[1], pdesc$ncase[1], pdesc$ncontrol[1])
  	  } else if (!is.na(pdesc$ncohort[1])) {
      	      sprintf('%s, n=%i', pdesc$label[1], pdesc$ncohort[1])
	  } else {
	      sprintf('%s, n=?', pdesc$label[1])
  	  }
  if (!is.null(xentity)) main <- paste(xentity$entity_label, main)
  fdesc <- gtxfilter_label(maf_ge = maf_ge, maf_lt = maf_lt, rsq_ge = rsq_ge, rsq_lt = rsq_lt, analysis = analysis)
  
  regionplot.new(chrom = chrom, pos_start = pos_start, pos_end = pos_end,
                 pmin = pmin, 
                 main = main,
                 protein_coding_only = protein_coding_only, 
                 dbc = dbc)
                                       
  if (identical(style, 'classic')) {
    ## best order for plotting
    pvals <- pvals[order(!is.na(pvals$impact), -log10(pvals$pval)), ]
    ## Plot all variants with VEP annotation as blue diamonds in top layer
    with(pvals, regionplot.points(pos, pval,
                                  pch = ifelse(!is.na(impact), 23, 21),
                                  col = ifelse(!is.na(impact), rgb(0, 0, 1, .75), rgb(.33, .33, .33, .5)),
                                  bg = ifelse(!is.na(impact), rgb(.5, .5, 1, .75), rgb(.67, .67, .67, .5))))
  } else if (identical(style, 'signals')) {
    ## best order for plotting
    pvals <- pvals[order(!is.na(pvals$impact), pvals$posterior), ]
    ## Plot variants coloured/sized by credible set, with VEP annotation is diamond shape in top layer
    with(pvals, regionplot.points(pos, pval,
                                  pch = ifelse(!is.na(impact), 23, 21), 
                                  cex = ifelse(c95, 1.25, .75), bg=colour, 
                                  col = ifelse(!is.na(impact), rgb(0, 0, 0, .5), rgb(.33, .33, .33, .5))))
    # legend indicating signals
    if (nsignals > 0) legend("bottomleft", pch = 21, col = rgb(.33, .33, .33, .5), pt.bg = colvec, legend=1:nsignals, horiz=T, bty="n", cex=.5)
  } 

  ## Highlight index SNP if pos or rs argument was used
  if (!missing(pos)) {
      ## funny workaround to avoid subset(pvals, pos=pos)
      ## means highlight won't work if the site is in gwas_results but not in sites
      gp <- sqlWrapper(dbc,
                       sprintf('SELECT chrom, pos FROM sites WHERE %s;',
                               gtxwhere(chrom = chrom, pos = pos)),
                       uniq = FALSE) # how to handle multiple sites at same pos?
      with(subset(pvals, pos == gp$pos),
           regionplot.points(pos, pval, pch = 1, cex = 2, col = "black"))
  }
  if (!missing(rs)) {
      gp <- sqlWrapper(dbc,
                       sprintf('SELECT chrom, pos, ref, alt FROM sites WHERE %s;',
                               gtxwhere(chrom = chrom, rs = rs))) # default uniq = TRUE
      # note gtxwhere() with combination of chrom and rs will cause sqlWrapper error if index site not on this chrom
      if (gp$chrom == chrom) {
          with(subset(pvals, pos == gp$pos & ref == as.character(gp$ref) & alt == as.character(gp$alt)),
               regionplot.points(pos, pval, pch = 1, cex = 2, col = "black"))
      }
  }

  return(invisible(NULL))
}

regionplot.data <- function(analysis,
                            chrom, pos_start, pos_end, pos, 
                            hgncid, ensemblid, rs, surround = 500000,
                            entity, 
                            style = 'signals',
                            maf_ge, rsq_ge, 
                            dbc = getOption("gtx.dbConnection", NULL)) {
    ## check database connection
    gtxdbcheck(dbc)
    
    ## Determine x-axis range from arguments
    xregion <- gtxregion(chrom = chrom, pos_start = pos_start, pos_end = pos_end, pos = pos, 
                         hgncid = hgncid, ensemblid = ensemblid, rs = rs, surround = surround,
                         dbc = dbc)
    chrom = xregion$chrom
    pos_start = xregion$pos_start
    pos_end = xregion$pos_end
    
    ## Determine entity, if required, for each analysis
    xentity <- gtxentity(analysis, entity = entity, hgncid = hgncid, ensemblid = ensemblid)
    
    if (identical(style, 'classic')) {
        ## basic query without finemapping
        pvals <- sqlWrapper(dbc,
                            sprintf('SELECT gwas_results.pos, gwas_results.ref, gwas_results.alt, pval, impact FROM %s.gwas_results LEFT JOIN vep USING (chrom, pos, ref, alt) WHERE %s AND %s AND %s AND %s AND pval IS NOT NULL;',
                                    gtxanalysisdb(analysis), 
                                    gtxwhat(analysis1 = analysis),
                                    if (!is.null(xentity)) sprintf('feature=\'%s\'', xentity$entity) else '(True)',
                                    gtxwhere(chrom, pos_ge = pos_start, pos_le = pos_end, tablename = 'gwas_results'),
                                    gtxfilter(maf_ge = maf_ge, rsq_ge = rsq_ge, analysis = analysis)),
                            uniq = FALSE)
    } else if (identical(style, 'signals')) {
        ## plot with finemapping annotation
        
        ## need to think more carefully about how to make sure any conditional analyses are either fully in, or fully out, of the plot
        pvals <- sqlWrapper(dbc,
                            sprintf('SELECT t1.chrom, t1.pos, t1.ref, t1.alt, beta, se, pval, signal, beta_cond, se_cond, pval_cond, impact FROM (SELECT gwas_results.chrom, gwas_results.pos, gwas_results.ref, gwas_results.alt, beta, se, pval, signal, beta_cond, se_cond, pval_cond FROM %s.gwas_results LEFT JOIN %s.gwas_results_cond USING (chrom, pos, ref, alt, analysis) WHERE %s AND %s AND %s AND %s AND pval IS NOT NULL) as t1 LEFT JOIN vep using (chrom, pos, ref, alt);',
                                    gtxanalysisdb(analysis), 
                                    gtxanalysisdb(analysis), 
                                    gtxwhat(analysis1 = analysis, tablename = 'gwas_results'),
                                    if (!is.null(xentity)) sprintf('feature=\'%s\'', xentity$entity) else '(True)',
                                    gtxwhere(chrom, pos_ge = pos_start, pos_le = pos_end, tablename = 'gwas_results'), 
                                    gtxfilter(maf_ge = maf_ge, rsq_ge = rsq_ge, analysis = analysis, tablename = 'gwas_results')),
                            uniq = FALSE)
    }
    pvals <- pvals[order(pvals$pval), ]
    return(list(region = xregion, analysis = analysis, entity = xentity, pvals = pvals))
}

regionplot.new <- function(chrom, pos_start, pos_end, pos, 
                           hgncid, ensemblid, rs, surround = 500000, 
                           pmin = 1e-10, main, 
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
  mtext(fdesc, 3, 0, cex = 0.5)

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
  return(with(sqlQuery(dbc, 
                       ## use SQL query to aggregate over multiple rows with same name to get whole span
                       sprintf('SELECT min(pos_start) AS pos_start, max(pos_end) AS pos_end, hgncid, ensemblid FROM genes WHERE %s %s GROUP BY ensemblid, hgncid ORDER BY pos_start', 
                               gtxwhere(chrom = chrom, pos_end_ge = pos_start, pos_start_le = pos_end),
			       if (protein_coding_only) 'AND genetype=\'protein_coding\'' else '')),
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

regionplot.recombination <- function(chrom, pos_start, pos_end, yoff = -.5, 
                                     dbc = getOption("gtx.dbConnection", NULL)) {
  gtxdbcheck(dbc)

  ## Recombination segments (r) that fall wholly or partly within [pos_start, pos_end]
  ## have r.pos_end >= pos_start and r.pos_start <= pos_end
  if (missing(pos_start)) pos_start = floor(par("usr")[1])
  if (missing(pos_end)) pos_end = ceiling(par("usr")[2])

  with(sqlQuery(dbc, 
                sprintf('SELECT pos_start, pos_end, recombination_rate FROM genetic_map WHERE %s', 
                        gtxwhere(chrom = chrom, pos_end_ge = pos_start, pos_start_le = pos_end))),
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


