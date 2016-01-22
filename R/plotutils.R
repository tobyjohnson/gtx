textgrid <- function(tmat,
                     x0 = 0, y0 = 0, x1 = x0 + 1, y1 = y0 + 1,
                     xalign = "l",
                     xjust = 0.5, yjust = 0, 
                     colsep = 3, rowsep = 1,
                     include.rownames, include.colnames,
                     draw = TRUE) {

  if (x1 < x0) {tmp <- x1; x1 <- x0; x0 <- tmp; rm(tmp)}
  if (y1 < y0) {tmp <- y1; y1 <- x0; y0 <- tmp; rm(tmp)}
  ## stop if x0==x1 or y0==y1

  ## coerce to matrix
  if (missing(include.rownames)) include.rownames <- !is.null(rownames(tmat))
  if (missing(include.colnames)) include.colnames <- !is.null(colnames(tmat))
  if (include.rownames) {
    tmat <- cbind(rownames(tmat), apply(tmat, 2, format))
  } else {
    tmat <- apply(tmat, 2, format)
  }
  if (include.colnames) {
    tmat <- rbind(colnames(tmat), tmat)
  }

  ## make alignment flags all "l" or "r" and of required length
  xalign <- rep(tolower(substr(xalign, 1, 1)), length.out = ncol(tmat))
  stopifnot(all(xalign %in% c("l", "r")))
  
  ## Remember previous setting
  oldcex <- par("cex")
  newcex <- oldcex
  
  ## Adaptively scale since actual strwidth drops in steps at discrete cex values
  while(TRUE) {
    xpos <- apply(tmat, 2, function(x) max(strwidth(x, units = "user")) +
                  strwidth("M", units = "user")*colsep)
    ypos <- apply(tmat, 1, function(x) max(strheight(x, units = "user")) +
                  strheight("M", units = "user")*rowsep)
    if (sum(xpos) <= (x1 - x0) && sum(ypos) <= (y1 - y0)) break
    newcex <- par("cex")/max(sum(xpos)/(x1 - x0), sum(ypos)/(y1 - y0))
    ## Check for stuck in infinite loop?
    if (newcex < 1e-4) stop("No font is small enough")
    par(cex = newcex)
  }

  
  
  xpos <- unname(x0 + (x1 - x0 - sum(xpos))*xjust + c(0, cumsum(xpos)))
  ypos <- unname(y1 - (y1 - y0 - sum(ypos))*yjust - c(0, cumsum(ypos)))
  xoff <- strwidth("M", units = "user")*colsep/2.
  yoff <- strheight("M", units = "user")*rowsep/2.

  draw <- rep(draw, length.out = ncol(tmat))
  for (idx in which(draw)) {
    if (xalign[idx] == "l") {
      text(rep(xpos[idx] + xoff, nrow(tmat)), ypos[-(nrow(tmat) + 1)] - yoff, 
           tmat[ , idx], adj = c(0, 1))
    } else if (xalign[idx] == "r") {
      text(rep(xpos[idx + 1] - xoff, nrow(tmat)), ypos[-(nrow(tmat) + 1)] - yoff, 
           tmat[ , idx], adj = c(1, 1))
    }
  }

  par(cex = oldcex)
  return(list(cex = newcex, xpos = xpos, ypos = ypos, xoff = xoff, yoff = yoff))
}

alphaize <- function(col, alpha = 0.5) {
  return(apply(col2rgb(col, alpha = FALSE),
               2,
               function(rgbvals) return(do.call(rgb, as.list(c(rgbvals/255, alpha = alpha))))))
}

legendgrid <- function(x, y = NULL,
                       legend,
                       xalign = "l", 
                       xjust = 0, yjust = 0, 
                       colsep = 3, rowsep = 1,
                       include.colnames, 
                       lty, lwd, col, pch, bg,
                       debug = FALSE) {

  ## currently allows pch to be a vector or matrix (to draw multiple symbols for each legend line)
  ## may break if used in unexpected ways
  ## FIXME should have separate control of line colour, and point colour and fill
  
  ## Calculate bounding box (x0, y0) (x1, y1) and alignment xjust yjust
  if (is.character(x) && identical(length(x), 1L)) {
    x <- tolower(x)
    ## Character values of x, "top", "topleft", "center" etc are all satisfied
    ## by setting the bounding box equal to the whole user plotting area
    ## and using different justification within that bounding box.
    ## The matching here is a bit fast and loose, consider requiring
    ## strict matches to "top", "topleft" etc.
    usr <- par("usr")
    x0 <- usr[1]
    x1 <- usr[2]
    y0 <- usr[3]
    y1 <- usr[4]
    xjust <- if (grepl("left", x, fixed = TRUE)) 0 else if (grepl("right", x, fixed = TRUE)) 1 else 0.5
    yjust <- if (grepl("top", x, fixed = TRUE)) 0 else if (grepl("bottom", x, fixed = TRUE)) 1 else 0.5
  } else {
    xy <- xy.coords(x, y)
    if (identical(length(xy$x), 1L)) {
      x0 <- xy$x[1]
      x1 <- par("usr")[2] # use right edge of whole user plotting area
    } else if (identical(length(xy$x), 2L)) {
      x0 <- xy$x[1]
      x1 <- xy$x[2]
    } else {
      stop("invalid 'x' argument")
    }
    if (identical(length(xy$y), 1L)) {
      y0 <- xy$y[1]
      y1 <- par("usr")[4] # use bottom edge of whole user plotting area
    } else if (identical(length(xy$y), 2L)) {
      y0 <- xy$y[1]
      y1 <- xy$y[2]
    } else {
      stop("invalid 'y' argument")
    }
  }
      
  if (missing(include.colnames)) include.colnames <- !is.null(colnames(legend))

  if (!missing(pch) && is.matrix(pch)) {
    ## FIXME check col and bg 
    em <- paste(rep("M", ncol(pch)), collapse = "")
  } else {
    em <- "M"
  }

  legend <- cbind(em, apply(legend, 2, format)) # make dummy column for drawing legend symbols
  xalign <- rep(tolower(substr(xalign, 1, 1)), length.out = ncol(legend))
  stopifnot(all(xalign %in% c("l", "r")))
  xalign <- c("l", xalign) # fix alignment flags for dummy column
  
  if (missing(lwd)) lwd <- par("lwd")
  if (missing(col)) col <- par("col")
  if (missing(bg)) bg <- par("bg")

  tmp <- textgrid(legend, x0, y0, x1, y1,
                  xalign = xalign, xjust = xjust, yjust = yjust, colsep = colsep, rowsep = rowsep,
                  include.rownames = FALSE, include.colnames = include.colnames,
                  draw = c(FALSE, rep(TRUE, ncol(legend) - 1)))

  if (debug) {
    abline(v = tmp$xpos, lty = "dotted")
    abline(h = tmp$ypos, lty = "dotted")
  }
  
  ypos <- tmp$ypos[1:nrow(legend) + if(include.colnames) 1 else 0] - tmp$yoff - 0.5*strheight("M", units = "user", cex = tmp$cex)
  
  for (idx in 1:length(ypos)) {
    if (!missing(lty)) {
      lines(tmp$xpos[1:2], rep(ypos[idx], 2), 
            cex = tmp$cex,
            lty = rep(lty, length.out = length(ypos))[idx],
            lwd = rep(lwd, length.out = length(ypos))[idx],
            col = rep(col, length.out = length(ypos))[idx])
    }
    if (!missing(pch)) {
      if (is.matrix(pch)) {
        xs <- seq(from = 0, to = 1, length.out = ncol(pch))
        points(xs*(tmp$xpos[1] + colsep/2*strwidth("M", units = "user", cex = tmp$cex)) +
               (1-xs)*(tmp$xpos[2] - colsep/2*strwidth("M", units = "user", cex = tmp$cex)), 
               rep(ypos[idx], ncol(pch)),
               cex = tmp$cex,
               pch = pch[idx, ],
               col = col[idx, ],
               bg = bg[idx, ])
      } else {
        points(mean(tmp$xpos[1:2]), ypos[idx],
               cex = tmp$cex,
               pch = rep(pch, length.out = length(ypos))[idx],
               col = rep(col, length.out = length(ypos))[idx], 
               bg = rep(bg, length.out = length(ypos))[idx])
      }
    }
  }
  return(tmp)
}

## plot.new()
## plot.window(0:1, 0:1)
## box()
## legendgrid(matrix(round(rnorm(100), 2), ncol = 10),
##            lty = "solid", pch = 1:10, col = rainbow(10))

## plot.new()
## plot.window(0:1, 0:1)
## plot(rnorm(100))
## box()
## legendgrid("bottomleft",
##            legend = as.data.frame(matrix(round(rnorm(9), 2), ncol = 3)),
##            xalign = rep(c("r", "l"), 5),
##            lty = "solid", pch = 1:10, col = rainbow(10))

# misc utility function for plotting

prettye <- function(x) {
  x <- as.character(x)
  return(sapply(x, function(x1) {
    x1s <- unlist(strsplit(x1, "e"))
    if (length(x1s) != 2) return(x1)
    return(eval(substitute(expression(MMM %*% 10^EEE), list(MMM = as.numeric(x1s[1]), EEE = as.numeric(x1s[2])))))
  }))
}

