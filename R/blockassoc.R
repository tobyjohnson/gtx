blockstats <- function(m1, m0, coefname = "GENOTYPE") {
  stopifnot(all.equal(class(m1), class(m0)[!class(m0) %in% c("coxph.null")]))
  ## should handle special case where class(m0)=="coxph.null"
  UseMethod("blockstats", m1)
}

blockstats.lm <- function(m1, m0, coefname = "GENOTYPE") {
  beta <- unname(coef(m1)[coefname])
  se <- unname(sqrt(vcov(m1)[coefname, coefname]))
  n <- length(na.omit(residuals(m1))) # seemingly no direct way to extract
  if(length(na.omit(residuals(m0))) != n) stop("Unequal sample size for m0 and m1!")
  pval <- pt(-abs(beta)/se, df = m1$df.residual)*2
  return(c(n = n, beta = beta, se = se, lrt = NA, pval = pval))
}

blockstats.glm <- function(m1, m0, coefname = "GENOTYPE") {
  beta <- unname(coef(m1)[coefname])
  se <- unname(sqrt(vcov(m1)[coefname, coefname]))
  n <- length(na.omit(residuals(m1))) # seemingly no direct way to extract
  if(length(na.omit(residuals(m0))) != n) stop("Unequal sample size for m0 and m1!")
  lrt <- max(m0$deviance - m1$deviance, 0)
  return(c(n = n, beta = beta, se = se, lrt = lrt, pval = NA))
}

blockstats.negbin <- function(m1, m0, coefname = "GENOTYPE") {
  beta <- unname(coef(m1)[coefname])
  se <- unname(sqrt(vcov(m1)[coefname, coefname]))
  n <- length(na.omit(residuals(m1))) # seemingly no direct way to extract
  if(length(na.omit(residuals(m0))) != n) stop("Unequal sample size for m0 and m1!")
  #lrt <- max(m0$deviance - m1$deviance, 0)
  lrt<- anova(m0, m1)[2,"LR stat."]
  return(c(n = n, beta = beta, se = se, lrt = lrt, pval = NA))
}

blockstats.coxph <- function(m1, m0, coefname = "GENOTYPE") {
  beta <- unname(coef(m1)[coefname])
  se <- unname(sqrt(vcov(m1)[coefname, coefname]))
  lrt <- max(2*(m1$loglik[length(m1$loglik)] - m0$loglik[length(m0$loglik)]), 0)
  n <- m1$n
  if(m0$n != n) stop("Unequal sample size for m0 and m1!")
  return(c(n = n, beta = beta, se = se, lrt = lrt, pval = NA))
}

blockstats.clm <- function(m1, m0, coefname = "GENOTYPE") {
  beta <- unname(coef(m1)[coefname])
  se <- unname(sqrt(vcov(m1)[coefname, coefname]))
  lrt <- max(2*(m1$logLik - m0$logLik), 0)
  n <- m1$n
  if(m0$n != n) stop("Unequal sample size for m0 and m1!")
  return(c(n = n, beta = beta, se = se, lrt = lrt, pval = NA))
}

blockassoc <- function(qcall, data, minimac,
                       usubjid = getOption("gtx.usubjid", "USUBJID"), 
                       threshold.MAF = 0, threshold.Rsq = 0,
                       threshold.pass = NULL, 
                       message.begin = "blockassoc", 
                       out.signif = 6, use.compiled = FALSE) {

  stopifnot(is.expression(qcall) || is.call(qcall))
  
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
  cat(message.begin, ": Analysing ", sum(data.include), "/", nrow(data),
      " subjects with non-missing phenotypes, covariates, and genotypes\n", sep = "")
  data <- data[data.include, , drop = FALSE]
  data.machid <- paste(data[[usubjid]], data[[usubjid]], sep = "->") # update 
  m0 <- eval(qcall, envir = data)
  cat(message.begin, ": Using ", length(data.machid), "/", nrow(dose),
      " genotyped subjects\n", sep = "")
  dose <- dose[match(data.machid, dose[ , 1]), 3:ncol(dose)]
  stopifnot(nrow(dose) == nrow(data))
  stopifnot(ncol(dose) == nrow(info))
  stopifnot(!any(is.na(dose)))

  ## this will not match the Rsq calculated by minimac, because minimac calculates the *ALLELIC* r2hat
  xbar <- colMeans(dose)
  x2bar <- colMeans(dose^2)
  info2 <- data.frame(analysed.Freq1 = signif(0.5*xbar, out.signif),
                      analysed.Rsq = signif((x2bar - xbar^2)/(xbar*(1-0.5*xbar)), out.signif),
                      SNP = info$SNP,
                      stringsAsFactors = FALSE)
  ## which to analyse?
  ## should there be an option to choose whether filtering on calculated or minimac-reported statistics?
  #cat(basename(minimac), ": Analysing n candidate variants\n", sep = "")
  ww <- which((pmin(info2$analysed.Freq1, 1 - info2$analysed.Freq1) >= threshold.MAF &
               info2$analysed.Rsq >= threshold.Rsq) |
              info2$SNP %in% threshold.pass)
  cat(message.begin, ": Analysing ", length(ww), "/", nrow(info2), " variants\n", sep = "")

  
  assoc <- as.data.frame(t(vapply(ww, function(idx) {
    return(tryCatch({
      data$GENOTYPE <- dose[ , idx]
      m1 <- suppressWarnings(update(m0, formula = . ~ . + GENOTYPE, data = data))
      blockstats(m1, m0, coefname = "GENOTYPE")},
                    error = function(e) return(c(NA, NA, NA, NA, NA))))},
                                  c(analysed.N = 0L, beta = 0., SE = 0., LRT = 0., pvalue = 0.))))
  if (all(is.na(assoc$pvalue))) {
    assoc$pvalue <- signif(pchisq(assoc$LRT, df = 1, lower.tail = FALSE), out.signif)
  }
  assoc <- within(assoc, {
    pwald <- signif(pnorm(-abs(beta)/SE)*2, out.signif)
    beta <- signif(beta, out.signif)
    SE <- signif(SE, out.signif)
    analysed.N <- as.integer(analysed.N)
    LRT <- signif(LRT, out.signif)
    pvalue <- signif(pvalue, out.signif)
  })
  return(cbind(info[ww, ], info2[ww, ], assoc))
}
