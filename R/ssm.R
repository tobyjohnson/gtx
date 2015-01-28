ssm.null <- function() return(NULL)

ssm.QT <- function(x) {
  ## The following free variables are used:
  if (FALSE) { # stop R CMD check complaining about no visible binding...
    sampleSize <- NA
    alleleFrequency <- NA
    dominanceCoeff <- NA
    effectSize <- NA
  }
  genotypeA <- rbinom(sampleSize, 2, alleleFrequency)
  return(data.frame(genotypeA = genotypeA,
                    genotypeR = ifelse(genotypeA == 2, 1, 0),
                    genotypeD = ifelse(genotypeA == 0, 0, 1),
                    genotypeH = ifelse(genotypeA == 1, 1, 0),
                    phenotype = scale((0:2 + c(0, dominanceCoeff, 0))[genotypeA + 1]*effectSize,
                                       center = TRUE, scale = FALSE) + rnorm(sampleSize)
                    ))
}
  
ssm.LM <- function(x, y) return(test.extract(lm(phenotype ~ genotypeA, data = y)))

ssm.CC <- function() {
  ## The following free variables are used:
  if (FALSE) { # stop R CMD check complaining about no visible binding...
    alleleFrequency <- NA
    dominanceCoeff <- NA
    effectSize <- NA
    populationAffected <- NA
    sampleAffected <- NA
  }
  ## populationAffected is proportion affected in a HW population with alleleFrequency
  freqVec <- c((1 - alleleFrequency)^2,
          2*(1 - alleleFrequency)*alleleFrequency,
          alleleFrequency^2)
  genoVec <- (0:2 + c(0, dominanceCoeff, 0)) * effectSize
  pgy1 <- function(alpha) return(freqVec / (1 + exp(-alpha - genoVec)))
  ## prob of G=AA,AB,BB given case
  pgy0 <- function(alpha) return(freqVec / (1 + exp(alpha + genoVec)))
  ## prob of G=AA,AB,BB given control
  ahat <- if (effectSize==0) {
    log(populationAffected/(1 - populationAffected))
  } else {
    uniroot(function(alpha) return(sum(pgy1(alpha)) - populationAffected),
            sort(log(populationAffected/(1 - populationAffected)*c(exp(-range(genoVec))))))$root
  }
  x <- data.frame(phenotype = rep(0:1, each = 3),
                  genotypeA = rep(0:2, times = 2),
                  genotypeR = rep(c(0, 0, 1), times = 2),
                  genotypeD = rep(c(0, 1, 1), times = 2),
                  genotypeH = rep(c(0, 1, 0), times = 2),
                  populationProb = c(pgy0(ahat), pgy1(ahat)),
                  sampleProb = c(pgy0(ahat)*(1 - sampleAffected)/(1 - populationAffected), pgy1(ahat)*sampleAffected/populationAffected))
  ## NCP1 is Pearson correlation between phenotype and genotypeA under sampleProb
  attr(x, "NCP1") <- ((sum(x$phenotype*x$genotypeA*x$sampleProb) - sum(x$phenotype*x$sampleProb)*sum(x$genotypeA*x$sampleProb))^2 /
                      ((sum(x$phenotype^2*x$sampleProb) - sum(x$phenotype*x$sampleProb)^2) *
                       (sum(x$genotypeA^2*x$sampleProb) - sum(x$genotypeA*x$sampleProb)^2)))
  return(x)
}

ssm.CCpopulation <- function(x) {
  ## The following free variables are used:
  if (FALSE) { # stop R CMD check complaining about no visible binding...
    sampleSize <- NA
  }
  ## unconditional simulation
  x$sample <- rmultinom(1, sampleSize, x$populationProb)
  return(x)
}

ssm.CCsample <- function(x) {
  ## The following free variables are used:
  if (FALSE) { # stop R CMD check complaining about no visible binding...
    sampleSize <- NA
  }
  ## conditional on affected marginally having expected sample counts
  x$sample <- 0
  for (aa in unique(x$phenotype)) {
    ii <- which(x$phenotype == aa)
    pp <- x$sampleProb[ii]
    x$sample[ii] <- rmultinom(1, sampleSize*sum(pp), pp) # automatically rounds N down
  }
  return(x)
}

ssm.GLM <- function(x, y) {
  return(test.extract(glm(phenotype ~ genotypeA, weights = sample, data = y, family = binomial)))
}

ssm.Surv <- function(x) {
  ## The following free variables are used:
  if (FALSE) { # stop R CMD check complaining about no visible binding...
    sampleSize <- NA
    alleleFrequency <- NA
    timeRecruit <- NA
    timeTotal <- NA
    dominanceCoeff <- NA
    effectSize <- NA
    hazardRate <- NA
  }
  stopifnot(timeTotal >= timeRecruit)
  genotypeA <- rbinom(sampleSize, 2, alleleFrequency)
  timeStart <- runif(sampleSize, 0, timeRecruit)
  timeFollowup <- timeTotal - timeStart
  timeEvent <- rexp(sampleSize,
                    rate = exp(scale((0:2 + c(0, dominanceCoeff, 0))[genotypeA + 1]*effectSize,
                      center = TRUE, scale = FALSE))*hazardRate)
  return(data.frame(genotypeA = genotypeA,
                    genotypeR = ifelse(genotypeA == 2, 1, 0),
                    genotypeD = ifelse(genotypeA == 0, 0, 1),
                    genotypeH = ifelse(genotypeA == 1, 1, 0),
                    phenotype = Surv(pmin(timeEvent, timeFollowup),
                                      ifelse(timeEvent <= timeFollowup, 1, 0))
                    ))
}

ssm.CoxPH <- function(x, y) {
  return(test.extract(coxph(phenotype ~ genotypeA, data = y)))
}
  
