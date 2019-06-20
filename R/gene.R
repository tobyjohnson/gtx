gene.label <- function(hgnc, ensg) {
    ## consider adding a db query option if hgnc is missing
    stopifnot(length(hgnc) == length(ensg))
    return(ifelse(!is.na(hgnc) & hgnc != '', hgnc, ensg)) # FIXME when ENSG ids change to integers, use gtxlabel()
}

#' @export
gene.annotate <- function(chrom, pos,
                          protein_coding_only = TRUE, 
                          dbc = getOption("gtx.dbConnection", NULL)) {
    stopifnot(length(chrom) == length(pos))
    gtxdbcheck(dbc)
    
    ## Current implementation involves a large number of SQL queries, which may be slow...
    if (identical(length(chrom), 0L)) return(NULL)
    
    ann <- sapply(1:length(chrom), function(idx) {
        res <- sqlWrapper(getOption('gtx.dbConnection_cache_genes', dbc), 
                          sprintf('SELECT hgncid, ensemblid FROM genes WHERE %s ORDER BY pos_start',
                                  ## FIXME filter on protein_coding_only
                                  gtxwhere(chrom = chrom[idx], pos_end_ge = pos[idx], pos_start_le = pos[idx])),
                          uniq = FALSE, zrok = TRUE)
        ## FIXME, need to support an alternative logic in gtxwhere to avoid serial queries
        ## If one more more genes overlap the query position, return like "[ABC123]" or "[ABC123;ABC124]"
        if (nrow(res) >= 1) {
            return(paste0('[', paste(unique(gene.label(res$hgncid, res$ensemblid)), collapse = ';'), ']'))
        }
        ## Otherwise, query nearest genes to left and right (arbitrary tie breaking, FIXME to be perfect)
        resl <- sqlWrapper(getOption('gtx.dbConnection_cache_genes', dbc),
                           sprintf('SELECT * FROM genes WHERE %s ORDER BY pos_end DESC LIMIT 1',
                                   ## FIXME filter on protein_coding_only
                                   gtxwhere(chrom = chrom[idx], pos_end_ge = pos[idx]-1000000, pos_end_le = pos[idx])),
                           zrok = TRUE)
        if (nrow(resl) == 0) {
            resl.s <- '---' # no genes within 1Mb
        } else {
            if (pos[idx] - resl$pos_end <= 10000) { # within 10kb
                resl.s <- paste0(gene.label(resl$hgnc, resl$ensemblid), '-')
            } else if (pos[idx] - resl$pos_end <= 100000) { # within 100kb
                resl.s <- paste0(gene.label(resl$hgnc, resl$ensemblid), '--')
            } else if (pos[idx] - resl$pos_end <= 1000000) { # within 1Mb
                resl.s <- paste0(gene.label(resl$hgnc, resl$ensemblid), '---')
            } else {
                stop('internal error in gene.annotate resl construction')
            }
        }
        resr <- sqlWrapper(getOption('gtx.dbConnection_cache_genes', dbc), 
                           sprintf('SELECT * FROM genes WHERE %s ORDER BY pos_start ASC LIMIT 1',
                                   ## FIXME filter on protein_coding_only
                                   gtxwhere(chrom = chrom[idx], pos_start_le = pos[idx]+1000000, pos_start_ge = pos[idx])),
                           zrok = TRUE)
        if (nrow(resr) == 0) {
            resr.s <- '---' # no genes within 1Mb
        } else {
            if (resr$pos_start - pos[idx] <= 10000) { # within 10kb
                resr.s <- paste0('-', gene.label(resr$hgnc, resr$ensemblid))
            } else if (resr$pos_start - pos[idx] <= 100000) { # within 100kb
                resr.s <- paste0('--', gene.label(resr$hgnc, resr$ensemblid))
            } else if (resr$pos_start - pos[idx] <= 1000000) { # within 1Mb
                resr.s <- paste0('---', gene.label(resr$hgnc, resr$ensemblid))
            } else {
                stop('internal error in gene.annotate resr construction')
            }
        }
        return(paste0(resl.s, '[]', resr.s))
    })
    return(ann)
}

#' @export
gene.nearest <- function(chr, pos, genetable) {
  stopifnot(length(chr) == length(pos))
  if (is.integer(chr)) chr <- as.character(chr)
  chr <- ifelse(substr(chr, 1, 3) == "chr", chr, paste("chr", chr, sep = ""))
  stopifnot(all(substr(chr, 1, 3) == "chr"))
  stopifnot(is.data.frame(genetable))
  stopifnot(all(c("name2", "chrom", "txStart", "txEnd") %in% names(genetable)))
  
  return (sapply(1:length(pos), function(idx) {
    gin <- unique(genetable$name2[which(genetable$chrom == chr[idx] & genetable$txStart+1 <= pos[idx] & genetable$txEnd+1 >= pos[idx])])
    if (length(gin) > 0) {
      return(paste(gin, collapse = ","))
    } else {
      ## only one gleft, gright using which.min
      gleft <- genetable$name2[which.min(ifelse(genetable$chrom == chr[idx] & genetable$txStart+1 <= pos[idx], pos[idx] - (genetable$txEnd+1), Inf))] # 
      gright <- genetable$name2[which.min(ifelse(genetable$chrom == chr[idx] & genetable$txEnd+1 >= pos[idx], genetable$txStart+1 - pos[idx], Inf))]
      return(paste(c(gleft, gright), collapse = "-"))
    }
  }))
}

#' @export
gene.draw <- function (chr, leftpos, rightpos, genetable, nodraw = NULL, genesep = 10000, 
                       hlines.min = NULL, yhi = -1, ylo = -5, exony = 0.05, genecex = 1) {
  if (is.integer(chr)) chr <- as.character(chr)
  if (substr(chr, 1, 3) != "chr") chr <- paste("chr", chr, sep = "")
  stopifnot(all(substr(chr, 1, 3) == "chr")) # should not be all as required to be length 1
  stopifnot(is.data.frame(genetable))
  stopifnot(all(c("name2", "chrom", "txStart", "txEnd", "exonStarts", "exonEnds") %in% names(genetable)))
  if (yhi < ylo) {
    ytmp <- yhi
    yhi <- ylo
    ylo <- ytmp
    rm(ytmp)
  }

  ## subset of genetable in the plotting region
  wg <- subset(genetable, chrom == chr & txStart + 1 <= rightpos & txEnd + 1 >=leftpos & !(name2 %in% nodraw))
  ## break into plotting groups same name and overlapping
  ## needs sort by name
  wg <- wg[with(wg, order(name2, txStart)), , drop = FALSE]
  if (nrow(wg) > 0) {
    if (nrow(wg) == 1) {
      wg$plotGroup <- 1
    } else {
      wg$plotGroup <- cumsum(c(TRUE, sapply(2:nrow(wg), function(idx) {
        return(all((wg$name2[1:(idx - 1)] != wg$name2[idx]) | (wg$txEnd[1:(idx - 1)] < wg$txStart[idx])))
      })))
    }

    ## reorder in plotting order
    wg$plotGroup <- rank(aggregate(txStart ~ plotGroup, data = wg, FUN = min)$txStart,
                         ties.method = "first")[wg$plotGroup]

    ## number of plotting groups
    gnum <- max(wg$plotGroup)
    ## autocalc hlines
    gstart <- aggregate(txStart ~ plotGroup, data = wg, FUN = min)$txStart
    gend <- aggregate(txEnd ~ plotGroup, data = wg, FUN = max)$txEnd
    if (gnum <= 1) {
      hlines <- 1
    } else {
      hlines <- min(which(c(sapply(1:(gnum - 1), function(hlines) {
        return(all(gend[1:(gnum - hlines)] + genesep < 
                   gstart[(1 + hlines):gnum]))
      }), TRUE)))
    }
    hlines <- max(c(hlines, hlines.min)) # hlines.min forces a minimum number of lines
    
    xmul <- 1
    line <- 0
    for (g1 in 1:gnum) {
      y <- yhi - (yhi - ylo) * line/hlines
      wwg <- which(wg$plotGroup == g1)
      txStart <- min(wg$txStart[wwg] + 1)
      txEnd <- max(wg$txEnd[wwg] + 1)
      txStrand <- paste(if ("-" %in% wg$strand[wwg]) 
                        "<"
      else NULL, if ("+" %in% wg$strand[wwg]) 
                        ">"
      else NULL, sep = "-")
      lines(c(txStart, txEnd) * xmul, rep(y, 2), lwd = 1)
      tmp <- unique(data.frame(exonStart = as.integer(unlist(strsplit(wg$exonStarts[wwg], 
                                 ","))) + 1, exonEnd = as.integer(unlist(strsplit(wg$exonEnds[wwg], 
                                               ","))) + 1))
      for (ex in 1:nrow(tmp)) {
        polygon(c(tmp$exonStart[ex], tmp$exonEnd[ex], tmp$exonEnd[ex], 
                  tmp$exonStart[ex]) * xmul, y + rep(c(-1, 1) * 
                                                     exony, each = 2), col = "yellow")
      }
      mp <- min(max((txStart + txEnd)/2, leftpos), rightpos)
      text(mp * xmul, y, paste(unique(wg$name2[wwg]), txStrand), pos = 1, offset = 0.3, 
           cex = genecex)
      line <- (line + 1)%%hlines
    }
  } else {
    # nothing to plot
  }
  return(invisible(NULL))
}

XXgene.draw <- function(chr, leftpos, rightpos, genetable, nodraw = NULL,
                      genesep = 10000, hlines.min = NULL,
                      yhi = -1, ylo = -5,
                      exony = 0.05, genecex = 1) {
  ## should make placement algorithm work better when the first line is entirely occupied by a very long transcript...
  if (is.integer(chr)) chr <- as.character(chr)
  if (substr(chr, 1, 3) != "chr") chr <- paste("chr", chr, sep = "")
  stopifnot(all(substr(chr, 1, 3) == "chr"))
  stopifnot(is.data.frame(genetable))
  stopifnot(all(c("name2", "chrom", "txStart", "txEnd", "exonStarts", "exonEnds") %in% names(genetable)))
  if (yhi < ylo) {ytmp <- yhi; yhi <- ylo; ylo <- ytmp; rm(ytmp)}
  
  genes <- setdiff(unique(genetable$name2[genetable$chrom == chr & 
                                genetable$txStart+1 <= rightpos &
                                genetable$txEnd+1 >= leftpos]), nodraw)
  gstart <- sapply(genes, function(gene) min(genetable$txStart[genetable$name2 == gene & genetable$chrom == chr]))
  gend <- sapply(genes, function(gene) max(genetable$txEnd[genetable$name2 == gene & genetable$chrom == chr]))
  genes <- genes[order(gstart)]
  gstart <- gstart[match(genes, names(gstart))]
  gend <- gend[match(genes, names(gend))]

  if (length(genes) <= 1) {
    hlines <- 1
  } else {
    hlines <- min(which(c(sapply(1:(length(genes)-1), function(hlines) {
      return(all(gend[1:(length(genes)-hlines)] + genesep < gstart[(1+hlines):length(genes)]))}), TRUE)))
  }
  hlines <- max(c(hlines, hlines.min))
  xmul <- 1 # could be set !=1 if e.g. plotting coords are Mb not bp
  line <- 0
  for (gene in genes) {
    y <- yhi - (yhi - ylo)*line/hlines
    wg <- which(genetable$name2 == gene & genetable$chrom == chr)
    txStart <- min(genetable$txStart[wg]+1)  # UCSC table coordinates are 0-based
    txEnd <- max(genetable$txEnd[wg]+1)  # UCSC table coordinates are 0-based
    txStrand <- paste(if ("-" %in% genetable$strand[wg]) "<" else NULL,
                      if ("+" %in% genetable$strand[wg]) ">" else NULL,
                      sep = "-")
    lines(c(txStart, txEnd) * xmul, rep(y, 2), lwd = 1)
    tmp <- unique(data.frame(exonStart = as.integer(unlist(strsplit(genetable$exonStarts[wg], ",")))+1,
                             exonEnd = as.integer(unlist(strsplit(genetable$exonEnds[wg], ",")))+1))
    for (ex in 1:nrow(tmp)) {
      polygon(c(tmp$exonStart[ex], tmp$exonEnd[ex], tmp$exonEnd[ex], tmp$exonStart[ex]) * xmul,
              y + rep(c(-1, 1) * exony, each = 2), col = "yellow")
    }
    mp <- min(max((txStart+txEnd)/2, leftpos), rightpos)
    text(mp * xmul, y, paste(gene, txStrand), pos = 1, offset = 0.3, cex = genecex)
    ## genecex==1 works well for CairoPNG(width=6*200,height=3*200,res=200)
    line <- (line+1) %% hlines
  }
}

#' @import data.table
#' @export
prune.distance <- function(data, surround = 500000L, sorted = FALSE) {
  if (nrow(data) == 0) return (NULL) # FIXME should return data table with zero rows
  data <- data.table::data.table(data)
  if (!sorted) {
      oo <- data[ , order(chrom, pos)]
      data <- data[oo, ]
  }
  cs <- c(1, which(data$chrom[-1] != data$chrom[-nrow(data)]) + 1) # row numbers of chrom starts (first row within each chrom block)
  dp <- c(NA, data$pos[-1] - data$pos[-nrow(data)]) # delta position
  dp[cs] <- NA # set dp to NA for all chrom start row
  ## if adjacent variants are greater than 2*surround bp apart,
  ## then when the signals are extended by +-surround on both sides,
  ## they will not overlap
  ssl <- is.na(dp) | dp > 2*surround # signal start at this row, logical
  ss <- which(ssl) # row indices where signals start
  data[ , signal := cumsum(ssl)] # index of signals, integers over whole genome, internal use within this function only
  signals <- data[ ,
                  list(chrom = unique(chrom), pos_start = min(pos), pos_end = max(pos),
                       num_variants = length(pos), min_pval = min(pval), row = which.min(pval)),
                  by = signal]
  if (sorted) {
      signals[ , row := row + ss - 1]
  } else {
      signals[ , row := NULL] # do not return row indices if input unsorted
      ## FIXME would be better to store original order and invert the permutation
  }
  signals[ , signal := NULL]
  return(signals)
}

#' @export
prune.genes <- function(chr, pos, genetable, genes, win.size = 10000) {
  stopifnot(length(chr) == length(pos))
  if (is.integer(chr)) chr <- paste("chr", chr, sep = "")
  stopifnot(all(substr(chr, 1, 3) == "chr"))
  prune <- rep(FALSE, length(chr))
  for (idx in which(genetable$name2 %in% genes)) {
    prune[chr == genetable$chrom[idx]
          & pos >= (genetable$txStart[idx]+1 - win.size)
          & pos <= (genetable$txEnd[idx]+1 + win.size)] <- TRUE
  }
  return(prune)
}
  


