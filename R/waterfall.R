waterfall <- function(object, x, data, ylab, ylim, bylevel = FALSE, col, ...) UseMethod("waterfall", object)

waterfall.character <- function(object, x, data, ylab, ylim, bylevel = FALSE, col, ...) {
  stopifnot(is.data.frame(data))
  stopifnot(object %in% names(data))
  stopifnot(x %in% names(data))
  if (missing(ylab)) ylab <- object
  y <- data[ , object]
  x <- data[ , x]
  include <- !is.na(y) & !is.na(x)
  ## exclude missing data BEFORE passing to internal function
  internal.waterfall(y[include], x[include],
                  ylab, ylim, bylevel, col)
}  

waterfall.numeric <- function(object, x, data, ylab, ylim, bylevel = FALSE, col, ...) {
  stopifnot(identical(length(object), length(x)))
  ## assert that y and x are vectors
  x <- as.factor(x)
  if (missing(ylab)) ylab <- deparse(substitute(object))
  include <- !is.na(object) & !is.na(x)
  internal.waterfall(object[include], x[include],
                     ylab, ylim, bylevel, col)
}

internal.waterfall <- function(y, x, ylab, ylim, bylevel = FALSE, col) {
  if (missing(col)) col <- rainbow(length(levels(x)))
  col <- rep(col, length.out = length(levels(x)))
  if (missing(ylim)) ylim <- extrange(y)
  o <- if (bylevel) order(x, -y) else order(-y)
  barplot(y[o], col = col[as.integer(x[o])], border = NA,
          las = 1, # horizontal y-axis labels
          ylab = ylab, ylim = ylim,
          xpd = FALSE) # xpd=FALSE keeps bars in region
  box()
  n <- table(x)
  legend("topright", pch=15, col = col, legend = paste(levels(x), "N =", n))
}

extrange <- function(x, ext = 0.04) {
  ## extended range needed to get barplot to look nice
  r <- range(x, na.rm = TRUE)
  return(r + c(-1, 1)*(r[2] - r[1])*ext)
}
