chi2ncp <- Vectorize(function(alpha, beta, df = 1) {
  if (is.na(alpha) || alpha < 0 || alpha > 1 || is.na(beta) || beta < 0 || beta > 1 || alpha > beta || is.na(df) || df < 0) return(NA)
  cv <- qchisq(alpha, df = df, lower.tail = FALSE) # critical value for chisq statistic
  ncpmax <- 1 # find upper limit by doubling until bracket achieved
  while (pchisq(cv, df = df, lower.tail = FALSE, ncp = ncpmax) < beta) ncpmax <- ncpmax*2  
  return(uniroot(function(ncp)
                 return(pchisq(cv, df = df, lower.tail = FALSE, ncp = ncp) - beta),
                 c(0, ncpmax))$root) # line search
})

chi2power <- Vectorize(function(alpha, ncp, df = 1) {
  if (is.na(alpha) || alpha < 0 || alpha > 1 || is.na(ncp) || ncp < 0 || is.na(df) || df < 0) return (NA)
  cv <- qchisq(alpha, df = df, lower.tail = FALSE) # critical value for chisq statistic
  return(pchisq(cv, df = df, lower.tail = FALSE, ncp = ncp))
})

fupower <- function(p0, n0, n1, beta0, invlog = FALSE, alpha = 0.05,
                    plty = c("solid", "dashed", "dotdash"), elty = "dotted",
                    pcol = "red", ecol = "blue") {
  chi0 <- qchisq(p0, df = 1, lower.tail = FALSE) # discovery chi2 statistic

  ## wlog assume discovery effect direction is positive
  xvec <- seq(from = min(0, sqrt(chi0) - 3), to = sqrt(chi0) + 3, length.out = 100)
  zcrit <- qnorm(alpha, lower.tail = FALSE) # for one-sided test

  n1 <- sort(na.omit(n1), decreasing = TRUE)

  plty <- rep(plty, length.out = length(n1))
  names(plty) <- as.character(n1)
  pcol <- rep(pcol, length.out = length(n1))
  names(pcol) <- as.character(n1)
  
  plot.new()
  plot.window(range(xvec), c(0, 1))
  abline(h = pretty(0:1), col = "grey")
  box()
  axis(2, las = 1)
  title(ylab = "Power")
  if (!missing(beta0)) {
    if (invlog) {
      pp <- pretty(exp(xvec/sqrt(chi0)*beta0))
      axis(1, at=log(pp)/beta0*sqrt(chi0), labels=pp)
      title(xlab = "Effect size (log scale)")
    } else {
      pp <- pretty(xvec/sqrt(chi0)*beta0)
      axis(1, at=pp/beta0*sqrt(chi0), labels=pp)
      title(xlab = "Effect size")
    }
  }
  
  pe <- dnorm(xvec, mean = sqrt(chi0))*sqrt(2*pi)*0.3
  lines(xvec, pe, lty = elty, col = ecol)
  
  ip <- sapply(n1, function(n) {
    pow <- pnorm(zcrit, mean = xvec*sqrt(n/n0), lower.tail = FALSE)
    lines(xvec, pow, lty = plty[[as.character(n)]], col = pcol[[as.character(n)]]) # side effect
    return(sum(pow*pe)/sum(pe))
  })

  plot.new()
  legend("topleft", bty = "o", 
         lty = c(plty, elty),
         col = c(pcol, ecol),
         legend = c(paste("N =", n1, "integr. power =", round(ip, 2)),
           "estimated effect size"))

  return(invisible(NULL))
}
