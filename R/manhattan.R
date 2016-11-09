contrasting.rainbow <- function (x, ...) {
  x <- as.integer(x)[1]
  stopifnot(x >= 1)
  ## make pretty colour cycle, by...
  ## ...finding non-unit divisors of ler...
  xnud <- setdiff(which(x %% 1:x == 0), 1)
  ## ...and thus largest coprime of x that is <=x/2...
  xcop <- max(which(sapply(1:(x/2), function(div) all(div %% xnud != 0))))
  ## ...and thus constrasting traversal of x
  xtrav <- (1:x*xcop) %% x + 1
  ## fallback; no coprime exists so use random shuffle
  if (!all(1:x %in% xtrav)) xtrav <- sample(x)
  return(rainbow(x, ...)[xtrav])
}

plotpos.by.chr <- function (chr, pos, gap = 5e7, chrset = c(1:22, "XY", "X", "Y", "M")) {
  stopifnot(length(pos) == length(chr))
  chr <- sub("^CHR", "", toupper(as.character(chr)))
  chrset <- toupper(as.character(chrset))
  plotpos <- rep(NA, length(pos))
  ox <- 0
  for (cx in chrset) {
    this <- which(chr == cx)
    if (length(this) > 0) {
      plotpos[this] <- pos[this] - min(pos[this], na.rm = TRUE) + ox + gap
      ox <- max(plotpos[this], na.rm = TRUE)
    }
  }
  return(plotpos)
}

plotcol.by.chr <- function (chr, cols = NULL, col.missing = "black", chrset = c(1:22, "XY", "X", "Y", "M")) {
  chrset <- toupper(as.character(chrset))
  ## We want to leave "gaps" for missing autosomes but not for missing sex/other
  chrnum <- match(sub("^CHR", "", toupper(as.character(chr))), chrset)
  cols <- if (is.null(cols)) contrasting.rainbow(length(chrset)) else rep(cols, length.out = length(chrset))
  chrnum[is.na(chrnum)] <- max(chrnum, na.rm = TRUE) + 1
  cols <- c(cols, rep(col.missing, max(0, max(chrnum) - length(cols))))
  return(cols[chrnum])
}

axis.by.chr <- function (chr, plotpos, side = 1, lines = 2) {
  stopifnot(length(plotpos) == length(chr))
  chr <- sub("^CHR", "", toupper(as.character(chr)))
  chrset <- unique(chr[!is.na(plotpos)])
  tickpos <- sapply(chrset, function (cx) return(0.5*(min(plotpos[chr == cx], na.rm = TRUE) + max(plotpos[chr == cx], na.rm = TRUE))))
  names(tickpos) <- chrset
  tickpos <- sort(tickpos)
  axis(side = side, at = tickpos, labels = FALSE)
  lidx <- rep(1:lines, length.out = length(tickpos))
  for (lx in unique(lidx)) mtext(names(tickpos[lidx == lx]), side = side, line = lx, at = tickpos[lidx == lx])
  return(tickpos)
}

manhattan <- function(p, SNP, chr, pos, pmin = NULL, ...) {
  if (!missing(SNP)) {
    chr <- vapply(strsplit(SNP, ":"), function(ss) return(ss[1]), character(1))
    pos <- as.integer(vapply(strsplit(SNP, "[:_]"), function(ss) return(ss[2]), character(1)))
  }
  plotpos <- plotpos.by.chr(chr, pos)
  plot(plotpos, -log10(p),
       col = plotcol.by.chr(chr, c("grey75", "cyan4")),
       ylim = c(0, max(-log10(c(p, pmin)), na.rm = TRUE)),
       xaxt = "n", yaxt = "n", ann = FALSE, 
       ...)
  axis(2, las = 1)
  axis.by.chr(chr, plotpos)
  title(xlab = "Genomic position by chromosome", 
        ylab = expression(-log[10](paste(italic(P), "-value"))))
  box()
  return(invisible(NULL))
}

  
