## functions for approximate Bayes factors

## wakefield ABF inverted to be BF for alternate model relative to null
## assumes prior mean is zero
abf.Wakefield <- function(beta, se, priorsd, log = FALSE) {
  if (log) {
    return(log(sqrt(se^2/(se^2 + priorsd^2))) + 
           (beta/se)^2/2 * priorsd^2/(se^2 + priorsd^2))
  } else {
    return(sqrt(se^2/(se^2 + priorsd^2)) * 
           exp((beta/se)^2/2 * priorsd^2/(se^2 + priorsd^2)))
  }
}

## code used for 2011 AJHG paper used fixed quadrature grid, xx <- seq(from = -0.5, to = 0.5, length.out = 1001)

## ABF for normal prior calculated numerically, should be equal to that calculated by abf.wakefield()
abf.normal <- function(beta, se, priorscale, gridrange = 3, griddensity = 20) {
  gridrange <- max(1, gridrange)
  griddensity <- max(1, griddensity)
  xlim <- range(c(beta - gridrange*se, beta + gridrange*se, -gridrange*priorscale, gridrange*priorscale), na.rm = TRUE)
  xx <- seq(from = xlim[1], to = xlim[2], length.out = max(3, ceiling(griddensity*(xlim[2] - xlim[1])/min(se, na.rm = TRUE))))
  dxx <- xx[2] - xx[1]
  return(sapply(1:length(beta), function(idx) {
    sum(dxx * dnorm(xx/priorscale)/priorscale * dnorm(beta[idx], xx, se[idx]))/dnorm(beta[idx], 0, se[idx])
  }))
}

## ABF for t distribution prior (a priori beta/priorscale ~ standard t distribution)
abf.t <- function(beta, se, priorscale, df = 1, gridrange = 3, griddensity = 20) {
  gridrange <- max(1, gridrange)
  griddensity <- max(1, griddensity)
  xlim <- range(c(beta - gridrange*se, beta + gridrange*se, -gridrange*priorscale, gridrange*priorscale), na.rm = TRUE)
  xx <- seq(from = xlim[1], to = xlim[2], length.out = max(3, ceiling(griddensity*(xlim[2] - xlim[1])/min(se, na.rm = TRUE))))
  dxx <- xx[2] - xx[1]
  return(sapply(1:length(beta), function(idx) {
    sum(dxx * dt(xx/priorscale, df = df)/priorscale * dnorm(beta[idx], xx, se[idx]))/dnorm(beta[idx], 0, se[idx])
  }))
}

norm1 <- function(x, log = FALSE) {
  if (all(is.na(x))) return(x)
  if (log) {
    x <- x - max(x, na.rm = TRUE)
    x <- exp(x)    
  } else {
    ## This does not work if x contains NaNs or +Infs
    stopifnot(all(x >= 0, na.rm = TRUE))
    x <- x / max(x, na.rm = TRUE)
  }
  return(x / sum(x, na.rm = TRUE))
}

credset <- function(bf, cred = 0.95) {
    if (sum(!is.na(bf)) == 0L) return(rep(NA, length(bf)))
    bf <- bf/max(bf, na.rm = TRUE)
    ps <- cumsum(sort(bf, decreasing = TRUE, na.last = NA))
    return(rank(-bf, na.last = TRUE) <= min(which(ps >= cred*ps[length(ps)])))
}

## Colocalisation analysis, implementing method of Giambartolomei et al. 2014
coloc.fast <- function(data, rounded = 6,
                       priorsd1 = 1, priorsd2 = 1, priorc1 = 1e-4, priorc2 = 1e-4, priorc12 = 1e-5, 
                       beta1 = 'beta_1', se1 = 'se_1', beta2 = 'beta_2', se2 = 'se_2') {
  abf1 <- abf.Wakefield(data[[beta1]], data[[se1]], priorsd1, log = TRUE)
  abf2 <- abf.Wakefield(data[[beta2]], data[[se2]], priorsd2, log = TRUE)
  w <- which(!is.na(abf1) & !is.na(abf2))
  nv <- length(w)
  abf1 <- norm1(c(0, abf1[w]), log = TRUE)
  abf2 <- norm1(c(0, abf2[w]), log = TRUE)
  res <- data.frame(hypothesis = paste0("H", 0:4),
                    label = c("No association",
                      "One variant associated with phenotype 1 only",
                      "One variant associated with phenotype 2 only",
                      "Two variants separately associated with phenotypes 1 and 2",
                      "One variant associated with phenotypes 1 and 2"),
                    prior = norm1(c(1, priorc1*nv, priorc2*nv, priorc1*priorc2*nv*(nv - 1), priorc12*nv)),
                    bf = c(abf1[1]*abf2[1], 
                      sum(abf1[-1])*abf2[1]/nv, 
                      abf1[1]*sum(abf2[-1])/nv, 
                      (sum(abf1[-1])*sum(abf2[-1]) - sum(abf1[-1]*abf2[-1]))/(nv*(nv - 1)), 
                      sum(abf1[-1]*abf2[-1])/nv))
  res$bf <- res$bf/res$bf[1]
  res$posterior <- norm1(res$prior*res$bf)
  if (is.finite(rounded)) {
    res$posterior = round(res$posterior, rounded)
  }
  return(list(results = res, nvariants = length(w)))
}
