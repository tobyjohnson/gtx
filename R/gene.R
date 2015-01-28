gene.annotate <- function(chr, pos, genetable, win.size = 10000) {
  stopifnot(length(chr) == length(pos))
  if (is.integer(chr)) chr <- as.character(chr)
  chr <- ifelse(substr(chr, 1, 3) == "chr", chr, paste("chr", chr, sep = ""))
  stopifnot(all(substr(chr, 1, 3) == "chr"))
  stopifnot(is.data.frame(genetable))
  stopifnot(all(c("name2", "chrom", "txStart", "txEnd") %in% names(genetable)))
  
  return (sapply(1:length(pos), function(idx) return(paste(unique(genetable$name2[which(genetable$chrom == chr[idx] & genetable$txStart+1 <= pos[idx]+win.size & genetable$txEnd+1 >= pos[idx]-win.size)]), collapse = ","))))
}

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

prune.distance <- function(chr, pos, pval, win.size = 500000) {
  stopifnot(length(pos) == length(chr))
  stopifnot(length(pval) == length(chr))
  if (length(chr) == 0) return (NULL)
  prune <- rep(NA, length(chr))
  for (idx in order(pval)) {
    prune[idx] <- all(c(Inf, abs(pos[idx] - pos[prune & chr == chr[idx]])) > win.size, na.rm = TRUE)
    ## note that prune == TRUE implies pval <= pval[idx]
  }
  return(prune)
}

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
  


