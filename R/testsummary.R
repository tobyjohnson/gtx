## this is a hack and could be made more object oriented!
coeff.extract <- function (object) UseMethod("coeff.extract", object)

coeff.extract.default <- function (object) {
  ## works for lm, glm, coxph, clm
  tmp <- data.frame(beta = coef(object), se = sqrt(diag(vcov(object))))
  colnames(tmp) <- c("Estimate", "Std. Error")
  return(tmp)
}

n.extract <- function (object) UseMethod("n.extract", object)

n.extract.default <- function (object) {
  ## works for coxph, clm
  return(n = tryCatch(object[["n"]], error = function(e) return(NA)))
}

n.extract.lm <- function (object) {
  ## works for lm, glm; coerces logicals to 0/1 for sum
  return(n = tryCatch(sum(!is.na(residuals(object))), error = function(e) return(NA)))
}

test.extract <- function(object, partial = NULL) UseMethod("test.extract", object)

test.extract.htest <- function(object, partial = NULL) {
  ## object with class "htest" is returned by chisq.test() and fisher.test()
  stopifnot(is.null(partial))
  ts <- if (object$method == "Fisher's Exact Test for Count Data") {
    if (!is.null(object$or)) {
      c(pval = object$p.value, effect = log(object$or), df = 1)
    } else {
      c(pval = object$p.value, effect = NA, df = NA) # !!!
    }
  } else if (object$method == "Pearson's Chi-squared test") {
    c(pval = object$p.value, effect = NA, df = object$parameter)
  } else if (object$method == "Welch Two Sample t-test") {
    c(pval = object$p.value, effect = object$estimate[2] - object$estimate[1], df = 1)
  } else if (object$method == "Pearson's product-moment correlation") {
    c(pval = object$p.value, effect = object$estimate, df = 1)
  } else {
    c(pval = object$p.value, effect = NA, df = NA) # !!!
  }
  names(ts) <- c("pval", "effect", "df"); return(ts)
}

test.extract.survdiff <- function(object, partial = NULL) {
  stopifnot(is.null(partial))
  df = length(object$n) - 1  # there must be a better way...
  ts <- c(pval = pchisq(object$chisq, df = df, lower.tail = FALSE), effect = NA, df = df)
  names(ts) <- c("pval", "effect", "df"); return(ts)
}

test.extract.lm <- function(object, partial = NULL) {
  objcoef <- coef(object)
  if (is.null(partial) && length(objcoef) == 2 && names(objcoef)[1] == "(Intercept)") partial <- names(objcoef)[2]
  ts <- if (!is.null(partial)) {
    tryCatch({
      objeffect <- objcoef[partial]
      stopifnot(!is.na(objeffect))
      objdrop <- drop1(object, partial, test = "Chisq")
      c(pval = objdrop[["Pr(>Chi)"]][2],
        effect = objeffect,
        df = objdrop[["Df"]][2])},
             error = function(e) return(c(pval = 1, effect = NA, df = 0)))
  } else {
    tryCatch({
      objf <- summary(object)$fstatistic
      c(pval = pf(objf["value"], objf["numdf"], objf["dendf"], lower.tail = FALSE),
        effect = NA, df = objf["numdf"])},
             error = function(e) return(c(pval = 1, effect = NA, df = 0)))
  }
  names(ts) <- c("pval", "effect", "df"); return(ts)
}                 
  
test.extract.glm <- function(object, partial = NULL) {
  objcoef <- coef(object)
  if (is.null(partial) && length(objcoef) == 2 && names(objcoef)[1] == "(Intercept)") partial <- names(objcoef)[2]
  ts <- if (!is.null(partial)) {
    ## same as for lm
    tryCatch({
      objeffect <- objcoef[partial]
      stopifnot(!is.na(objeffect))
      objdrop <- drop1(object, partial, test = "Chisq")
      c(pval = objdrop[["Pr(>Chi)"]][2],
        effect = objeffect,
        df = objdrop[["Df"]][2])},
             error = function(e) return(c(pval = 1, effect = NA, df = 0)))
  } else {
    tryCatch({
      df <- object$df.null - object$df.residual
      c(pval = pchisq(object$null.deviance - object$deviance, df = df, lower.tail = FALSE),
        effect = NA, df = df)},
             error = function(e) return(c(pval = 1, effect = NA, df = 0)))
  }
  names(ts) <- c("pval", "effect", "df"); return(ts)
}

test.extract.coxph <- function(object, partial = NULL) {
  objcoef <- coef(object)
  if (is.null(partial) && length(objcoef) == 1) partial <- names(objcoef)[1]
  ts <- if (!is.null(partial)) {
    ## same as for lm
    objeffect <- objcoef[partial]
    tryCatch({
      stopifnot(!is.na(objeffect))
      ## suppressWarnings needed when object has only one term, which is dropped
      objdrop <- suppressWarnings(drop1(object, partial, test = "Chisq"))
      c(pval = objdrop[["Pr(>Chi)"]][2],
        effect = objeffect,
        df = objdrop[["Df"]][2])},
             error = function(e) return(c(pval = 1, effect = NA, df = 0)))
  } else {
    ## compute LRT for all terms
    tryCatch({
      df <- sum(!is.na(object$coefficients)) # there must be a better way...
      c(pval = pchisq(2*(object$loglik[2] - object$loglik[1]), df = df, lower.tail = FALSE),
        effect = NA, df = df)},
             error = function(e) return(c(pval = 1, effect = NA, df = 0)))
  }
  names(ts) <- c("pval", "effect", "df"); return(ts)
}

