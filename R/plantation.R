plantation <- function(ntext, beta, ci.lo, ci.hi, se, alpha = 0.05,
                       FUN = I, pvals = TRUE, meta = TRUE,
                       xzero = 0, xlim, xticks, xlab = "Effect",
                       digits = 2,
                       groups = list("Fixed effect meta-analysis" = 1:length(beta))) {

  ## calculate CI limits from SE or vice versa
  if (!missing(ci.lo) && !missing(ci.hi) && missing(se)) {
    stopifnot(identical(length(ci.lo), length(beta)))
    stopifnot(identical(length(ci.hi), length(beta)))
    se <- (ci.hi - ci.lo)/(qnorm(alpha/2, lower.tail = FALSE) -  qnorm(alpha/2))
  } else if (!missing(se)) {
    stopifnot(identical(length(se), length(beta)))
    ci.lo <- beta + qnorm(alpha/2)*se
    ci.hi <- beta + qnorm(alpha/2, lower.tail = FALSE)*se
  } else {
    stop("Either ci.lo and ci.hi, or se, must be specified")
  }
  
  w <- se^-2
  if (meta && length(groups) > 0) {  
    ## calculate meta-analyses for all specified groups
    metas <- sapply(groups, function(g) {
      stopifnot(all(g %in% 1:length(beta)))
      beta.m <- sum(w[g]*beta[g])/sum(w[g])
      se.m <- sum(w[g])^-0.5
      return(c(beta = beta.m,
               se = se.m,
               ci.lo = beta.m + qnorm(alpha/2)*se.m,
               ci.hi = beta.m + qnorm(1 - alpha/2)*se.m))})
  } else {
    meta <- FALSE
    groups <- list()
  }
  
  ## Note effect sizes and CIs assumed to be similar scale and hence
  ## use format(digits=digits), but P-values may be very different magnitudes
  ## and hence use elementwise signif(digits=digits)
  ## apply single format call to all numeric values

  nvals <- format(as.double(FUN(c(beta, if (meta) metas["beta", ] else NULL,
                                  ci.lo, if (meta) metas["ci.lo", ] else NULL,
                                  ci.hi, if (meta) metas["ci.hi", ] else NULL))),
                  digits = digits)
  nn <- length(beta) + length(groups)
  btext <- paste(nvals[1:nn],
                 " [", nvals[(nn + 1):(2*nn)],
                 "-", nvals[(2*nn + 1):(3*nn)], "]", sep = "")
  ptext <- sapply(pnorm(-abs(c(beta, if (meta) metas["beta", ] else NULL))/
                        c(se, if (meta) metas["se", ] else NULL))*2,
                  function(p1) prettye(format(signif(p1, digits = digits))))
  ## Force 1 dp of precision for weights, no possibility for user choice
  wtext <- paste(format(round(100*w/sum(w), digits = 1)), "%", sep = "")

  if (missing(xlim)) xlim <- range(c(ci.lo, ci.hi), na.rm = TRUE) + c(-1, 1)*0.1*max(se)
  if (!missing(xticks)) xlim <- range(c(xlim, xticks))
  ## symmetric, xlim <- c(-1, 1)*max(abs(xlim))
  ## xlim <- range(pretty(xlim))
  
  em <- strwidth("M", units = "inches")*3
  nwid <- max(strwidth(c(ntext, names(groups)), units = "inches")) + em
  bwid <- max(strwidth(c(paste(xlab, "[95% CI]"), btext), units = "inches")) + em
  pwid <- if (pvals) max(strwidth(c("P-value", ptext), units = "inches")) + em else 0
  wwid <- if (meta) max(strwidth(c("Weight", wtext), units = "inches")) + em else 0

  if (nwid + bwid + pwid + wwid > 0.75*par("fin")[1]) {
    oldcex <- par("cex")
    newcex <- signif(oldcex*(0.75*par("fin")[1])/(nwid + bwid + pwid + wwid), 2)
    par(cex = newcex)
    warning ("Text would occupy more than three quarters of plot width; changing cex from ", oldcex, " to ", newcex)
    nwid <- nwid*newcex/oldcex
    bwid <- bwid*newcex/oldcex
    pwid <- pwid*newcex/oldcex
    wwid <- wwid*newcex/oldcex
  }
  nwid <- nwid/par("fin")[1]
  bwid <- bwid/par("fin")[1]
  pwid <- pwid/par("fin")[1]
  wwid <- wwid/par("fin")[1]
  
  fwid <- 1 - nwid - bwid - pwid - wwid
   
  xv <- c(-nwid/fwid*xlim[2] + (nwid/fwid + 1)*xlim[1],
          xlim,
          (bwid/fwid + 1)*xlim[2] - bwid/fwid*xlim[1],
          ((bwid+pwid)/fwid + 1)*xlim[2] - (bwid+pwid)/fwid*xlim[1],
          ((bwid+pwid+wwid)/fwid + 1)*xlim[2] - (bwid+pwid+wwid)/fwid*xlim[1])
  plot.new()
  plot.window(range(xv), c(-0.5, nn))
  ## abline(v = xv, col = "red") ## for debugging layout

  text(xv[3], nn, paste(xlab, "[95% CI]"), pos = 4)
  if (pvals) text(xv[4], nn, "P-value", pos = 4)
  if (meta) text(xv[5], nn, "Weight", pos = 4)
  
  for (idx in 1:(length(beta) + if (meta) length(groups) else 0)) {
    text(xv[1], nn - idx, c(ntext, names(groups))[idx], pos = 4)
    text(xv[3], nn - idx, btext[idx], pos = 4)
    if (pvals) text(xv[4], nn - idx, ptext[idx], pos = 4)
    if (meta) text(xv[5], nn - idx, wtext[idx], pos = 4)
  }  
  for (idx in 1:length(beta)) {
    ## Grey square
    ## what to do if symbol is outside xlim?  Currently not plotted
    if (beta[idx] >= xlim[1] & beta[idx] <= xlim[2]) {
      points(beta[idx], nn - idx, pch = 15, cex = 3*sqrt(w[idx]/max(w)), col = "grey")
    } # else draw nothing
    ## Left (lower) whisker
    if (ci.lo[idx] >= xlim[1]) {
      if (ci.lo[idx] < xlim[2]) {
        lines(c(min(beta[idx], xlim[2]), ci.lo[idx]), rep(nn - idx, 2))
      } # else draw nothing
    } else {
      if (beta[idx] >= xlim[1]) {
        arrows(min(beta[idx], xlim[2]), nn - idx, x1 = xlim[1], length = .05)
      } # else draw nothing
    }
    ## Right (upper) whisker
    if (ci.hi[idx] <= xlim[2]) {
      if (ci.hi[idx] > xlim[1]) {
        lines(c(max(beta[idx], xlim[1]), ci.hi[idx]), rep(nn - idx, 2))
      } # else draw nothing
    } else {
      if (beta[idx] <= xlim[2]) {
        arrows(max(beta[idx], xlim[1]), nn - idx, x1 = xlim[2], length = .05)
      } # else draw nothing
    }
    ## Tick at point estimate
    if (beta[idx] >= xlim[1] & beta[idx] <= xlim[2]) {
      points(beta[idx], nn - idx, pch = 3)
    } else { # point estimate out of range, print "warning triangle"
      warning("point estimate out of range")
      if (beta[idx] < xlim[1]) {
        points(xlim[1], nn - idx, pch = 24, cex = 3, col = "red", bg = "white")
        text(xlim[1], nn - idx, "!", col ="red")
      } else if (beta[idx] > xlim[2]) {
        points(xlim[2], nn - idx, pch = 24, cex = 3, col = "red", bg = "white")
        text(xlim[2], nn - idx, "!", col ="red")
      }
    }
  }
  if (meta) {
    for (idx in 1:ncol(metas)) {
      ## need to handle out of range nicely
      polygon(c(metas["ci.lo", idx], metas["beta", idx], metas["ci.hi", idx], metas["beta", idx]), ncol(metas) - idx + c(0, 0.2, 0, -0.2), col = "grey")
    }
  }

  if (identical(FUN, I)) {
    if (missing(xticks)) xticks <- axisTicks(xlim, log = FALSE)
  } else if (identical(FUN, exp)) {
    if (missing(xticks)) xticks <- log(axisTicks(xlim, log = TRUE))
  } else {
    ## hmmm
  }
  xticks <- xticks[xticks >= xlim[1] & xticks <= xlim[2]]
  axis(1, at = xticks, labels = FUN(xticks))
  mtext(xlab, 1, 2,
        at = 0.5*(min(xlim) + max(xlim)), # may be broken
        cex = par("cex")) # apparently does not follow global cex setting  
  abline(v = xzero, lty = "dotted")
  return(invisible(NULL))
}
