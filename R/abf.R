## functions for approximate Bayes factors

#'
#' Calculate approximate Bayes factor (ABF) using method of Wakefield (2009)
#'
#'   Calculates an approximation to the Bayes factor for an alternative
#'   model where the parameter beta is a priori normal, against a smaller
#'   model where beta is zero, by approximating
#'   the likelihood function with a normal distribution.
#'
#'   See \dQuote{Bayes factors for genome-wide association studies:
#'     comparison with P-values} by John Wakefield, 2009, Genetic
#'   Epidemiology 33(1):79-86 at \url{http://dx.doi.org/10.1002/gepi.20359}.
#'
#'   The original definition was for the Bayes factor for the model where
#'   beta is zero (no association), relative to the model with a normal
#'   prior for beta (association).  The definition used for this function
#'   inverts the original definition, so that higher values of the Bayes
#'   factor (or log Bayes factor) indicate stronger support for
#'   the model with association.
#'
#'   For strong associations, calculating Bayes factors may be numerically
#'   unstable, and it is recommended to work on a log scale and rescale
#'   appropriately before attempting to calculate Bayes factors or
#'   posterior probabilities.
#'
#'
#' @param beta 	Vector of effect size estimates.
#' @param se 		Vector of associated standard errors.
#' @param priorsd 	Scalar specifying the standard deviation of the prior on true effect sizes.
#' @param log 		Whether to return results on a natural log scale.
#'
#' @return
#'    A vector of approximate Bayes factors, on a log scale if
#'    \code{log=TRUE}.  Higher values indicate stronger support for
#'    association (which is inverted relative to the original definition).
#'
#'
#' @examples
#' data(agtstats)
#' agtstats$pval <- with(agtstats, pchisq((beta/se.GC)^2, df = 1, lower.tail = FALSE))
#' max1 <- function(bf) return(bf/max(bf, na.rm = TRUE))
#' agtstats$BF.normal <- with(agtstats, max1(abf.Wakefield(beta, se.GC, 0.05)))
#' agtstats$BF.t <- with(agtstats, max1(abf.t(beta, se.GC, 0.0208)))
#' with(agtstats, plot(-log10(pval), log(BF.normal)))
#' with(agtstats, plot(-log10(pval), log(BF.t)))
#'
#' @author   Toby Johnson \email{Toby.x.Johnson@gsk.com}
#'
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
