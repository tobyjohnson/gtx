otplot <- function(object, x, data, ylab, xlab, col, col.mix = c("white", "black"), style = "percent", ycrop = FALSE, yext = 1.25) UseMethod("otplot", object)

otplot.character <- function(object, x, data, ylab, xlab, col, col.mix = c("white", "black"), style = "percent", ycrop = FALSE, yext = 1.25) {
  stopifnot(is.data.frame(data))
  stopifnot(object %in% names(data))
  stopifnot(x %in% names(data))
  if (missing(ylab)) ylab <- object
  if (missing(xlab)) xlab <- x
  y <- as.factor(data[ , object])
  x <- as.factor(data[ , x])
  include <- !is.na(y) & !is.na(x)
  ## exclude missing data BEFORE passing to internal function
  internal.otplot(y[include], x[include],
                  ylab, xlab, col, col.mix, style, ycrop, yext)
}  

otplot.numeric <- function(object, x, data, ylab, xlab, col, col.mix = c("white", "black"), style = "percent", ycrop = FALSE, yext = 1.25) {
  stopifnot(identical(length(object), length(x)))
  ## assert that y and x are vectors
  if (missing(ylab)) ylab <- deparse(substitute(object))
  if (missing(xlab)) xlab <- deparse(substitute(x)) # broken!
  object <- as.factor(object)
  x <- as.factor(x)
  include <- !is.na(object) & !is.na(x)
  internal.otplot(object[include], x[include],
                  ylab, xlab, col, col.mix, style, ycrop, yext)
}

otplot.factor <- function(object, x, data, ylab, xlab, col, col.mix = c("white", "black"), style = "percent", ycrop = FALSE, yext = 1.25) {
  stopifnot(identical(length(object), length(x)))
  ## assert that y and x are vectors
  if (missing(ylab)) ylab <- deparse(substitute(object))
  if (missing(xlab)) xlab <- deparse(substitute(x)) # broken!
  x <- as.factor(x)
  include <- !is.na(object) & !is.na(x)
  internal.otplot(object[include], x[include],
                  ylab, xlab, col, col.mix, style, ycrop, yext)
}

internal.otplot <- function(y, x, ylab, xlab, col, col.mix, style, ycrop, yext) {

  ## Blend col for each level of x, with col.mix to span range of y, to make matrix of colours
  if (missing(col))
    col <- rainbow(length(levels(x)))
  col <- rep(col, length.out = length(levels(x)))
  if (!identical(length(col.mix), 2L)) stop ("invalid col.mix")
  colbars <- do.call(cbind, lapply(col, function(col1)
                                   sapply(seq(from = 0.125, to = 0.875, length.out = length(levels(y))),
                                          function(i) {
                                            if (i <= 0.5) {
                                              do.call(rgb, as.list(((1-2*i)*col2rgb(col.mix[1]) + 2*i*col2rgb(col1))/255))
                                            } else {
                                              do.call(rgb, as.list(((2-2*i)*col2rgb(col1) + (2*i-1)*col2rgb(col.mix[2]))/255))
                                            }
                                          })))
  
  yt <- table(y, x)
  style <- tolower(style)
  if (style == "percent") {
    yt <- apply(yt, 2, function(n) 100*n/sum(n))
    ylim <- c(0, 100)
    ## ytlab = label for entries of yt, as opposed to ylab = label for levels of y
    ytlab <- paste("Percentage within", xlab)
  } else if (style == "percentall") {
    yt <- 100*yt/sum(yt)
    ylim <- c(0, 100)
    ytlab <- "Percentage overall"
  } else if (style == "fraction") {
    yt <- apply(yt, 2, function(n) n/sum(n))
    ylim <- c(0, 1)
    ytlab <- paste("Fraction within", xlab)
  } else if (style == "fractionall") {
    yt <- yt/sum(yt)
    ylim <- c(0, 1)
    ytlab <- "Fraction overall"
  } else if (style == "count") {
    ## yt <- yt # no change
    ylim <- range(pretty(c(0, max(yt))))
    ytlab <- "Count"
  } else {
    stop("invalid style")
  }
  if (ycrop) ylim <- range(pretty(c(0, max(yt))))
  barplot(yt, beside = TRUE, col = colbars, 
          xaxt = "n",
          ylim = ylim*yext, yaxt = "n")
  xat <- (length(levels(y)) + 1) * (1:length(levels(x)) - 0.5) + 0.5
  mtext(levels(x), side = 1, line = 1, at = xat) # because axis(1) may not to include all labels
  n <- table(x)
  mtext(paste("N=", n, sep = ""), side = 1, line = 2, at = xat)
  axis(2, at = pretty(ylim), las = 1)
  box()
  title(xlab = xlab, ylab = ytlab)
          
  foo <- matrix(levels(y), ncol = 1); colnames(foo) <- ylab
  legendgrid("topright",
             legend = foo,
             pch = matrix(rep(22, nrow(colbars)*ncol(colbars)), ncol = ncol(colbars)),
             col = matrix(rep("black", nrow(colbars)*ncol(colbars)), ncol = ncol(colbars)),
             bg = colbars, debug = FALSE, colsep = 1)
  
  return(invisible(NULL))
}

