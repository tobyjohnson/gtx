## mps = "Modular Power Simulations|System"
##
## idea of a "driver" with interchangeable "bits"
##   driver combines "bits" of simulation code to produce power estimates
## powersmoother
##   combine power estimates across design points to estimate smooth power function
##   [add option to force power=alpha at null, spoof with power=alpha, w=1000000 at null; must be option since test may miscalibrate]
## powerblender
##   integrate power estimates across design points
##   [make multivariate...]

## should tFun "know" any of the parameters?  Or force it not to "know" by design?

mps.driver1 <- function(design1, xFun, yFun, tFun, nrep = 1000) {
  return(with(as.list(design1), {
    environment(xFun) <- environment()
    environment(yFun) <- environment()
    ## note tFun does not "know" design1
    x1 <- xFun()
    y1 <- yFun(x1)
    ## check all tFun return values are equal length
    z1 <- lapply(tFun, function(tFun1) return(tFun1(x1, y1)))
    check.pval <- names(z1)[sapply(z1, function(tt) !("pval" %in% names(tt)))]
    if (length(check.pval) > 0) stop("some tFun do not return pval : ", paste(check.pval, collapse = ","))
    check.length <- sapply(z1, length)
    numstat <- unique(check.length)
    if (length(numstat) != 1) stop("tFun do not all return equal length : ", paste(names(z1), "[", check.length, "]", sep = "", collapse = ","))
    ## generate zz using replicate
    replicate(nrep, {
      y1 <- yFun(x1)
      vapply(tFun, function(tFun1) return(tFun1(x1, y1)), double(numstat))})}))
  ## Should check str and dim of return value, but numstat not in scope
}

mps.summary <- function(zz, alpha = 0.05, sided = "two-sided") {
  stopifnot(length(dim(zz)) == 3)
  return(do.call(rbind,
                 lapply(1:dim(zz)[2], function(tFunIdx) {
                   thisP <- zz["pval", tFunIdx, , drop = TRUE]
                   if ("effect" %in% dimnames(zz)[[1]]) {
                     thisEffect <- zz["effect", tFunIdx, , drop = TRUE]
                   } else {
                     thisEffect <- rep(NA, length(thisP))
                   }
                   zok <- !is.na(thisP)
                   res <- cbind(#design[rep(designIdx, length(alpha)), , drop = FALSE],
                                data.frame(test = rep(dimnames(zz)[[2]][tFunIdx], length(alpha)),
                                           nsim = rep(sum(zok), length(alpha)),
                                           alpha = alpha,
                                           sided = sided,
                                           empirical.power = rep(NA, length(alpha))),
                                if (dim(zz)[1] > 1) {
                                  do.call("cbind", lapply(2:(dim(zz)[1]), function(statIdx) {
                                    thisStat <- zz[statIdx, tFunIdx, zok]
                                    tmp2 <- data.frame(rep(mean(thisStat), length(alpha)),
                                                       rep(var(thisStat), length(alpha)))
                                    names(tmp2) <- paste(c("mean", "var"), dimnames(zz)[[1]][statIdx], sep = ".")
                                    tmp2
                                  }))
                                } else NULL) # end cbind
                   ## in subset of results with zok and df == max(df)
                   ## NCPhat = max(mean(qchisq(thisP, lower=FALSE, df = df)) - df, 0)
                   ## See Saxena & Alam 1982 (Annals of Statistics 10(3):1012-1016)
                   ## concerning theoretical inadmissibility of MLE
                   ## ## In practical terms, using optim[ize] to obtain MLE
                   ## often falls over because dchisq(0, ... ) is often Inf, occurs for P==1
                   res$empirical.power <- vapply(1:length(alpha), function(adx) {
                     if (sided[adx] == "two-sided") {
                       workP <- thisP
                     } else if (sided[adx] == "greater") {
                       workP <- ifelse(thisEffect >= 0, thisP/2, 1 - thisP/2)
                     } else if (sided[adx] == "less") {
                       workP <- ifelse(thisEffect <= 0, thisP/2, 1 - thisP/2)
                     }
                     return(mean(workP[zok] <= alpha[adx]))
                   }, double(1))
                   return(res)
                 })))
}

mps.driver <- function(design, xFun, yFun, tFun,
                       nrep = 1000,
                       alpha = 0.05, sided = "two-sided",
                       verbose = FALSE) {
  ## 
  ## check arguments and coerce to correct types
  ## 
  if (is.null(xFun)) xFun <- function() return(NULL)
  stopifnot(is.function(xFun))
  if (is.null(yFun)) yFun <- function() return(NULL)
  stopifnot(is.function(yFun))
  if (is.function(tFun)) tFun <- list(test = tFun)
  tFun <- as.list(tFun)
  stopifnot(length(tFun) >= 1)
  stopifnot(all(sapply(tFun, function(tFun1) is.function(tFun1))))
  if (is.null(names(tFun))) names(tFun) <- if (length(tFun) == 1) "test" else paste("test", 1:length(tFun), sep = "")
  design <- as.data.frame(design)
  stopifnot(nrow(design) >= 1)
  nrep <- as.integer(nrep)
  stopifnot(nrep >= 1)
  alpha <- as.vector(alpha)
  stopifnot(length(alpha) >= 1)
  sided <- c("two-sided", "greater", "less")[match(rep(tolower(substr(as.vector(sided), 1, 1)), length.out = length(alpha)), c("t", "g", "l"))]
  if (any(is.na(sided))) stop ("sided must be two-sided, greater, or less (or abbrevations thereof)")

  ##
  ## loop over design points, running simulations and summarising results
  ##
  tmp <- do.call(rbind,
                 ## Simulate independently for each row of design
                 lapply(1:nrow(design), function(designIdx) { 
                   ## For each row of the design, generate zz
                   ## zz is a numstat*length(tFun)*nrep array of test statistics
                   ## typically numstat==3 for pval, effect, df
                   zz <- mps.driver1(design[designIdx, , drop = FALSE], xFun, yFun, tFun, nrep = nrep)
                   if (verbose) cat("Done simulations for", designIdx, "/", nrow(design), "design points\n")
                   return(cbind(# replicate this row of the design
                                design[rep(designIdx, dim(zz)[[2]]*length(alpha)), , drop = FALSE], 
                                # summary of test summaries
                                mps.summary(zz, alpha = alpha, sided = sided)))
                 }))
  rownames(tmp) <- 1:nrow(tmp)
  return(tmp)
}


#                                    maxDf <- max(thisDf[zok])
#                                    dfmax <- zok & thisDf == maxDf
#                                    ncp.ML <- if (any(dfmax)) optimize(function(qq) return(-sum(dchisq(thisStat[dfmax], 
#                                                                                                       df = thisDf[dfmax], 
#                                                                                                       ncp = qq, log = TRUE))),
#                                                                       c(0, max(thisStat[dfmax]) + 1))$minimum else NA


## use expand.grid for making design matrix
## strongly suggest use verbose names to avoid clashes
## provided functions use:  sampleSize, effectSize, alleleFrequency, populationAffected, sampleAffected
##
## you can compare two or more tests on the same set of simulated datasets;
## each test MUST produce a vector<double> of the same length with the same names
## element pval is required; effect is required for one-sided tests
## the built in examples all use (pval, effect, df)

## chi2df ----> testsummary
## for simulation based infrastructure, should build in support for:
##   1. exact tests (test returns p, not chi2) without spoofing equivalent chi2
##   2. test fails (no variation in explanatory or response variable, no events in coxph)
##   3. alpha=c(.05,.005,.0005), alt=c("two-sided", "greater", "greater")
##   4. expectation and variance of effect size estimate (check simulation algorithm, estimate bias and efficiency)
##   5. df>1 tests with no "direction"
##   6. "clever" estimation of chi2 NCP [do we NEED this?]
##
## This is a powerdriver so pval is most important
## ----> test.extract must always return(c(pval, effectsize, df))
## a generic function that works as an "adapter" to make different tests report necessary results in a consistent form
##
## lots of choices about most "useful" default behaviour.  Basic guidance is
## if partial == NULL do the (practical yet most) exact multi-df test against intercept-only model (not a Wald test unless absolutely necessary)
## if partial != NULL do a 1df test in most practical way (use drop1(test="Chisq", not a Wald test unless neccesary)
##
## if user wants other partial tests they have to write their own code



                                        #                                                            vapply(adx, function(adx) return(mean(pchisq(thisStat[zok], df = thisDf[zok], lower = FALSE) <=# a)), double(1)),
                                        #                                                            power.ncp = pchisq(qchisq(alpha, df = maxDf, lower.tail = FALSE), df = maxDf, ncp = ncp.ML, low#er.tail = FALSE))
                                        #                                    

# can spoof analytical calcs by using rows of X to enumerate possibilities, with columnhaving element 

## user supplies functions xFun(), yFun(X), zFun(X, Y) to be called sequentially
## typically, xFun() generates independent variables or does deterministic calcs
## typically  yFun(X) generates dependent variables or simulates
## typically  tFun(X, Y) generates test
## these functions need additional arguments specifying the CURRENT point in design space, e.g. sample size, effect size, allele freq
## how to pass?  s <- design[current, , drop = FALSE]
## Options considered:
## 1.    design vars are *free variables* in Xfun(), ...
##       1a.     with(s, {X <- Xfun(); Y <- Yfun(X); ...)
##       1b.     with(s, {environment(Xfun) <- environment(); X <- Xfun(); ...}
## 2.    design vars are passed explicitly
##       2a.     s$X <- do.call(Xfun, s); s$Y <- do.call(Yfun, s); z$Z <- do.call(Zfun, s)
##       2b.     s$X <- Xfun(s); s$Y <- Yfun(s); Z <- Zfun(s)
##
## 1a. ; Does not work because R uses lexical scoping, 
## 1b. Makes programming Xfun() etc easy; but could be error prone
## 2   s stores all information; but requires reinitialisation of s each time
## 2a. Xfun etc need to be defined with Xfun(maf, n, ...) to absorb extra args
## 2b. Xfun etc need to be like "Xfun <- function(p) with(p, FN)"


## Do we need to accomodate situation where df of test done depends on simulation
## e.g. low MAF and genotypic model will be mix of 1 and 2 df tests

## alpha not part of design because does not require indep simulations

## option to allow two-sided 



## powersmoother <- function(power, effectsize, w) {
##   if (length(power) <= 2) return(power)
##   stopifnot(length(effectsize) == length(power))
##   w <- rep(w, length.out = length(power))
##   return(pnorm(predict(glm(power~effectsize, family = binomial(link = probit), w = w))))
## }

## powerblender <- function(power, effectsize, effectsizehat, effectsizese) {
##   if (length(power) <= 1) return(power)
##   stopifnot(length(effectsize) == length(power))
##   o <- order(effectsize)
##   effectsize <- effectsize[o]
##   power <- power[o]
##   effectsized <- effectsize[-1] - effectsize[-length(effectsize)]
##   qw <- dnorm(effectsize, effectsizehat, effectsizese)*(c(0, effectsized) + c(effectsized, 0))/2 # quadrature weights for generalized trapezium rule
##   if (sum(qw) < .95) warning("design points (effectsize) do not well cover parameter uncertainty")
##   return(list(power = sum(power*qw)/sum(qw), design.cover = sum(qw)))
## }

mps.asymptotic <- function(design, alpha) {
  tmp <- do.call(rbind,
                 lapply(1:nrow(design), function(designIdx) { 
                   with(as.list(design[designIdx, , drop = FALSE]), {
                     ncp <- tryCatch({
                       freqVec <- c((1 - alleleFrequency)^2,
                                    2*(1 - alleleFrequency)*alleleFrequency,
                                    alleleFrequency^2)
                       genoVec <- (0:2 + c(0, dominanceCoeff, 0)) * effectSize
                       sampleSize*(sum(freqVec*genoVec^2)-sum(freqVec*genoVec)^2)},
                                     error = function(e) return(NA))
                     do.call(rbind, lapply(alpha, function(alpha1) {
                       return(c(powerContinuous = tryCatch(chi2power(alpha1, ncp = ncp), error = function(e) return(NA)),
                                powerBinary = tryCatch(chi2power(alpha1, ncp = ncp*sampleAffected*(1 - sampleAffected)), error = function(e) return(NA)),
                                powerSurvival = tryCatch(chi2power(alpha1, ncp = ncp*(1 - surviveEnd)), error = function(e) return(NA))))}))
                   })}))}
## needs cbinding up with alpha values and row of design...


