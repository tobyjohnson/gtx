blockstats <- function(m1, m0, coefname = "GENOTYPE") {
  stopifnot(all.equal(class(m1), class(m0)))
  UseMethod("blockstats", m1)
}

blockstats.lm <- function(m1, m0, coefname = "GENOTYPE") {
  beta <- coef(m1)[coefname]
  se <- sqrt(vcov(m1)[coefname, coefname])
  n <- length(na.omit(m1$residuals)) # best way to extract
  pval <- pt(-abs(beta)/se, df = m1$df.residual)*2
  return(c(beta = beta, se = se, n = n, lrt = NA, pval = pval))
}

blockstats.glm <- function(m1, m0, coefname = "GENOTYPE") {
  beta <- coef(m1)[coefname]
  se <- sqrt(vcov(m1)[coefname, coefname])
  n <- length(na.omit(residuals(m1))) # seemingly no direct way to extract#
  lrt <- max(m0$deviance - m1$deviance, 0)
  return(c(beta = beta, se = se, n = n, lrt = lrt, pval = NA))
}

blockstats.coxph <- function(m1, m0, coefname = "GENOTYPE") {
  beta <- coef(m1)[coefname]
  se <- sqrt(vcov(m1)[coefname, coefname])
  lrt <- max(2*(m1$loglik[2] - m0$loglik[2]), 0)
  n <- m1$n
  return(c(beta = beta, se = se, n = n, lrt = lrt, pval = NA))
}

blockstats.clm <- function(m1, m0, coefname = "GENOTYPE") {
  beta <- coef(m1)[coefname]
  se <- sqrt(vcov(m1)[coefname, coefname])
  lrt <- max(2*(m1$logLik - m0$logLik), 0)
  n <- m1$n
  return(c(beta = beta, se = se, n = n, lrt = lrt, pval = NA))
}

blockassoc <- function(qcall, data, minimac,
                       usubjid = "USUBJID",
                       threshold.MAF = 0, threshold.Rsq = 0, 
                       out.signif = 6, use.compiled = FALSE) {
  stopifnot(is.call(qcall))

  if (use.compiled) {
    stop("Compiled functions not yet implemented")
    ## test for special cases that can be handled by compiled C++ code
    ## PUT ALL THIS INSIDE A tryCatch()
    if (all.equal(qcall[1], quote(lm()))) { # qcall is a call to lm()
      mc <- match.call(definition = lm, call = qcall, expand.dots = TRUE)
      ## only proceed if all arguments are ones we know how to handle:
      if (all(names(mc)[-1] %in% c("formula", "subset"))) {
        mc[1] <- quote(model.frame())
        mf <- eval(mc, envir = data) #rownames contain indices after subset and NA omission
        ## call my C++ f(data[as.integer(rownames(mf)), usubjid], mf)
      }
      ## fall through
    }
    ## fall through
    ## otherwise fall through to interpreted code
  }

  m0 <- eval(qcall, envir = data)

  if (!missing(minimac)) {
    ## should loop over multiple files if supplied
    info <- read.table(gzfile(paste(minimac, "info.gz", sep = ".")),
                       comment.char = "", quote = "", header = TRUE, as.is = TRUE)
    dose <- read.table(gzfile(paste(minimac, "dose.gz", sep = ".")),
                       comment.char = "", quote = "", header = FALSE,
                       colClasses = c(rep("character", 2), rep("numeric", nrow(info))))
    stopifnot(all(dose[ , 2] == "DOSE"))
  } else {
    stop("Missing minimac argument; other genotype formats not supported")
  }
  
  stopifnot(length(unique(data[[usubjid]])) == nrow(data))
  data.machid <- paste(data[[usubjid]], data[[usubjid]], sep = "->")

  data.include <- apply(!is.na(data), 1, all) & data.machid %in% dose[ , 1]
  cat("Analysing ", sum(data.include), "/", nrow(data),
      " subjects with non-missing phenotypes, covariates, and genotypes\n", sep = "")
  data <- data[data.include, , drop = FALSE]
  data.machid <- paste(data[[usubjid]], data[[usubjid]], sep = "->") # update 

  cat("Using ", length(data.machid), "/", nrow(dose),
      " genotyped subjects\n", sep = "")
  dose <- dose[match(data.machid, dose[ , 1]), 3:ncol(dose)]
  stopifnot(nrow(dose) == nrow(data))
  stopifnot(ncol(dose) == nrow(info))
  stopifnot(!any(is.na(dose)))

  ## this will not match the Rsq calculated by minimac, because minimac calculates the *ALLELIC* r2hat
  xbar <- colMeans(dose)
  x2bar <- colMeans(dose^2)
  info2 <- data.frame(analysed.Freq1 = signif(0.5*xbar, out.signif),
                      analysed.Rsq = signif((x2bar - xbar^2)/(xbar*(1-0.5*xbar)), out.signif))
  ## which to analyse?
  ww <- which(pmin(info2$analysed.Freq1, 1 - info2$analysed.Freq1) >= threshold.MAF & info2$analysed.Rsq >= threshold.Rsq)
  cat("Analysing ", length(ww), "/", nrow(info2), " variants with MAF >= ", threshold.MAF, " and Rsq >= ", threshold.Rsq, "\n", sep = "")
  
  assoc <- as.data.frame(t(vapply(ww, function(idx) {
    return(tryCatch({
      data$GENOTYPE <- dose[ , idx]
      m1 <- suppressWarnings(update(m0, formula = . ~ . + GENOTYPE, data = data))
      blockstats(m1, m0, coefname = "GENOTYPE")},
                    error = function(e) return(NA, NA, NA, NA, NA)))},
                                  c(beta.allele = 0., SE.allele = 0., N.analysed= 0, LRT.allele = 0., pvalue.allele = 0.))))
  if (all(is.na(assoc$pvalue.allele))) {
    assoc$pvalue.allele <- signif(pchisq(assoc$LRT.allele, df = 1, lower.tail = FALSE), out.signif)
  }
  assoc <- within(assoc, {
    pwald.allele <- signif(pnorm(-abs(beta.allele)/SE.allele)*2, out.signif)
    beta.allele <- signif(beta.allele, out.signif)
    SE.allele <- signif(SE.allele, out.signif)
    N.analysed <- signif(N.analysed, out.signif)
    LRT.allele <- signif(LRT.allele, out.signif)
    pvalue.allele <- signif(pvalue.allele, out.signif)
  })
  return(cbind(info[ww, ], info2[ww, ], assoc))
}
