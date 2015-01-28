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
