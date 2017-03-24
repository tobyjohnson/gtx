## functions for approximate Bayes factors

## wakefield ABF inverted to be BF for alternate model relative to null
## assumes prior mean is zero
abf.Wakefield <- function(beta, se, priorsd) return(sqrt(se^2/(se^2 + priorsd^2)) * exp((beta/se)^2/2 * priorsd^2/(se^2+priorsd^2)))

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

## Colocalisation analysis, implementing method of Giambartolomei et al. 2014
coloc.fast <- function(data, rounded = 6,
                      priorsd1 = 1, priorsd2 = 1, priorc1 = 1e-4, priorc2 = 1e-4, priorc12 = 1e-5, 
                      beta1 = 'beta_1', se1 = 'se_1', beta2 = 'beta_2', se2 = 'se_2') {
    abf1 <- abf.Wakefield(data[[beta1]], data[[se1]], priorsd1)
    abf2 <- abf.Wakefield(data[[beta2]], data[[se2]], priorsd2)
    w <- which(!is.na(abf1) & !is.na(abf2))
    abf1 <- abf1[w] / (1 + sum(abf1[w]))
    abf2 <- abf2[w] / (1 + sum(abf2[w]))
    res <- data.frame(hypothesis = paste0("H", 0:4),
                      posterior = c((1 - sum(abf1))*(1 - sum(abf2)), 
                                    priorc1 * sum(abf1)*(1 - sum(abf2)), 
                                    priorc2 * (1 - sum(abf1))*sum(abf2), 
                                    priorc1*priorc2 * (sum(abf1)*sum(abf2) - sum(abf1*abf2)), 
                                    priorc12 * sum(abf1*abf2)),
                      label = c("No association",
                                "One variant associated with phenotype 1 only",
                                "One variant associated with phenotype 2 only",
                                "Two variants separately associated with phenotypes 1 and 2",
                                "One variant associated with phenotypes 1 and 2"))
    res$posterior <- res$posterior/sum(res$posterior)
    if (is.finite(rounded)) res$posterior = round(res$posterior, rounded)
    return(res)
}
