qtplot <- function(object, x, data, ylab, xlab, ylim, col, points.method = "beeswarm", boxplot.layer = "over", ...) UseMethod("qtplot", object)

qtplot.character <- function(object, x, data, ylab, xlab, ylim, col, points.method = "beeswarm", boxplot.layer = "over", ...) {
  stopifnot(is.data.frame(data))
  stopifnot(object %in% names(data))
  stopifnot(x %in% names(data))
  if (missing(ylab)) ylab <- object
  if (missing(xlab)) xlab <- x
  y <- data[ , object]
  x <- as.factor(data[ , x])
  include <- !is.na(y) & !is.na(x)
  ## exclude missing data BEFORE passing to internal function
  internal.qtplot(y[include], x[include],
                  ylab, xlab, ylim, col, points.method, boxplot.layer)
}  

qtplot.numeric <- function(object, x, data, ylab, xlab, ylim, col, points.method = "beeswarm", boxplot.layer = "over", ...) {
  stopifnot(identical(length(object), length(x)))
  ## assert that y and x are vectors
  if (missing(ylab)) ylab <- deparse(substitute(object))
  if (missing(xlab)) xlab <- deparse(substitute(x)) # broken!
  x <- as.factor(x)
  include <- !is.na(object) & !is.na(x)
  internal.qtplot(object[include], x[include],
                  ylab, xlab, ylim, col, points.method, boxplot.layer)
}

internal.qtplot <- function(y, x, ylab, xlab, ylim, col, points.method, boxplot.layer) {
  if (missing(ylim)) ylim <- range(y)
  if (missing(col)) col <- rainbow(length(levels(x)))
  col <- rep(col, length.out = length(levels(x)))
  plot.new()
  plot.window(c(0.5, length(levels(x)) + 0.5), ylim)
  axis(2, las = 1) # horizontal y-axis labels
  mtext(levels(x), side = 1, line = 1, at = 1:length(levels(x))) # because axis(1) may not to include all labels
  n <- table(x)
  mtext(paste("N=", n, sep = ""), side = 1, line = 2, at = 1:length(n))
  box()
  title(xlab = xlab, ylab = ylab)
  if (!identical(boxplot.layer, "under") && !identical(boxplot.layer, "over") && !identical(boxplot.layer, "none")) {
    stop("qtplot called with invalid boxplot.layer, must be either 'under', 'over', or 'none'")
  }
  if (identical(boxplot.layer, "under")) {
    boxplot(y ~ x,
            outpch = NA, # do not draw outliers as will be drawn by stripchart
            add = TRUE, xaxt = "n", yaxt = "n") # need to disable x and y axes even with add=TRUE
  }
  if (identical(points.method, "stripchart")) {
    stripchart(y ~ x,
               vertical = TRUE, method = "jitter", jitter = 0.45, # use 90% of available width in each bin 
               add = TRUE, 
               xaxt = "n", yaxt = "n",
               col = col, pch = 1)
  } else if (identical(points.method, "beeswarm")) {
    beeswarm(y ~ x,
             corral = "random", # ensure points stay inside their bins
             add = TRUE,
             xaxt = "n", yaxt = "n",
             col = col)
  } else {
    stop("qtplot called with invalid points.method, must be either 'stripchart' or 'beeswarm'")
  }
  if (identical(boxplot.layer, "over")) {
    boxplot(y ~ x,
            outpch = NA, # do not draw outliers as will be drawn by stripchart
            add = TRUE, xaxt = "n", yaxt = "n") # need to disable x and y axes even with add=TRUE
  }
  return(invisible(NULL))
}
