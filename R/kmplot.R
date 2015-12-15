kmplot <- function(object, x, data, ylab, xlab, ylim, xlim, col, lty,
                   legend.location = "topright", legend.cex = 1,
                   legend.median = FALSE, digits = 1, atrisk = FALSE, ...) UseMethod("kmplot", object)

kmplot.character <- function(object, x, data, ylab, xlab, ylim, xlim, col, lty,
                             legend.location = "topright", legend.cex = 1,
                             legend.median = FALSE, digits = 1, atrisk = FALSE, ...) {
  stopifnot(is.data.frame(data))
  stopifnot(object %in% names(data))
  stopifnot(x %in% names(data))
  if (missing(ylab)) ylab <- object
  # it will be an error for xlab to be missing
  stopifnot(is.Surv(data[ , object]))
  y <- data[ , object]
  x <- as.factor(data[ , x])
  include <- !is.na(y) & !is.na(x)
  ## exclude missing data BEFORE passing to internal function
  internal.kmplot(y[include], x[include],
                  ylab, xlab, ylim, xlim, col, lty,
                  legend.location, legend.cex, legend.median, digits, atrisk, ...)
}  

kmplot.Surv <- function(object, x, data, ylab, xlab, ylim, xlim, col, lty,
                        legend.location = "topright", legend.cex = 1, 
                        legend.median = FALSE, digits = 1, atrisk = FALSE, ...) {
  stopifnot(identical(nrow(object), length(x)))
  ## assert that x is a vector
  x <- as.factor(x)
  if (missing(ylab)) ylab <- deparse(substitute(object))
  include <- !is.na(object) & !is.na(x)
  internal.kmplot(object[include], x[include],
                  ylab, xlab, ylim, xlim, col, lty,
                  legend.location, legend.cex, legend.median, digits, atrisk, ...)
}

internal.kmplot <- function(y, x, ylab, xlab, ylim, xlim, col, lty,
                            legend.location, legend.cex, legend.median, digits, atrisk, ...) {
  xcounts <- table(x)
  wlevels <- which(xcounts > 0) # which levels of x will be included in the fit
  nlevels <- length(levels(x)) # total levels including those with zero count
  sf <- survfit(y ~ x)
  if (length(wlevels) != 1) {
    stopifnot(identical(length(sf$strata), length(wlevels)))
  }
  
  ## Ensure space for axis labels including at-risk below x axis
  oldmar <- par("mar")
  if (atrisk) {
    par(mar = pmax(oldmar, c(nlevels + 4.1, 3.1, 0, 0)))
  } else {
    par(mar = pmax(oldmar, c(3.1, 3.1, 0, 0)))
  }

  if (missing(xlab)) xlab <- "Time (unspecified units)"
  if (missing(xlim)) xlim <- c(0, max(sf$time))
  if (missing(ylim)) ylim <- NULL # cannot pass missing argument into plot.survfit
  if (missing(col)) col <- rainbow(nlevels)
  if (missing(lty)) lty <- 1:6
  col <- rep(col, length.out = nlevels)
  lty <- rep(lty, length.out = length(col))

  plot(sf, conf.int = FALSE,
       lwd = 1, col = col[wlevels], lty = lty[wlevels], 
       las = 1, ylim = if (missing(ylim)) NULL else ylim, xlim = xlim,
       ann = FALSE, ...)
  
  n <- table(x)
  if (identical(tolower(legend.location), "none")) {
    ## Don't plot legend
  } else {
    if (legend.median) {
      sftf <- format(round(summary(sf)$table[ , c("median", "0.95LCL", "0.95UCL"), drop = FALSE], digits = digits))
      sfl <- cbind(levels(x), xcounts, sftf[ , "median"],
                   paste("(", sftf[ , "0.95LCL"], "-", sftf[ , "0.95UCL"], ")", sep = ""))
      colnames(sfl) <- c("", "n", "Median", "(95% CI)")
      oldcex <- par("cex")
      par(cex = legend.cex)
      legendgrid(legend.location, legend = sfl,
                 xalign = c("l", "r", "r", "r"), colsep = 1, 
                 lwd = 1, col = col, lty = lty)
      par(cex = oldcex)
    } else {
      legend(legend.location, lwd = 1, col = col, lty = lty,
             legend = paste(levels(x), "N =", xcounts),
             cex = legend.cex)
    }
  }
  
  mtext(ylab, side = 2, line = 3)
  
  if (atrisk) {
    mtext(xlab, side = 1, line = 2)
    mtext("Numbers at risk:", side = 1, line = 3, at = min(xlim), adj = 0)
    xat <- pretty(xlim)
    xat <- xat[xat >= xlim[1] & xat <= xlim[2]]
    
    if (length(wlevels) == 1) {
      sf$strata <- sf$n # make dummy stratum for computing numbers at risk
    }
    for (jdx in 1:length(wlevels)) {
      rr <- (cumsum(c(0, sf$strata))[jdx] + 1):(cumsum(sf$strata)[jdx])
      # yright=0 makes extrapolation on RHS zero, f=1 makes approximation right-continuous
      mtext(approx(sf$time[rr], sf$n.risk[rr], xat, method = "constant", yright = 0, rule = 2, f = 1)$y,
            side = 1, line = jdx + 3, at = xat, col = col[wlevels[jdx]])
    }
  } else {
    mtext(xlab, side = 1, line = 3)
  }
 
  par(mar = oldmar)
  return(invisible(NULL))
}
