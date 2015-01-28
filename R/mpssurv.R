trial.params <- Vectorize(function(hazardRate = NA, surviveEnd = NA, timeTotal = NA, timeRecruit = NA) {
  ## check valid parameter values
  if (!is.na(hazardRate) && hazardRate <= 0) stop("hazardRate must be positive")
  if (!is.na(surviveEnd) && (surviveEnd <= 0 || surviveEnd >= 1)) stop("surviveEnd must be between 0 and 1")
  if (!is.na(timeTotal) && timeTotal <= 0) stop("timeRecruit must be positive")
  if (!is.na(timeRecruit) && timeRecruit < 0) stop("timeRecruit must be non-negative")
  if (!is.na(timeTotal) && !is.na(timeRecruit) && timeRecruit > timeTotal) stop("timeRecruit must not be greater than timeTotal")
  if (is.na(hazardRate) + is.na(surviveEnd) + is.na(timeTotal) + is.na(timeRecruit) != 1) stop("One parameter must be NA")

  if (is.na(hazardRate) && !is.na(surviveEnd) & !is.na(timeTotal) & !is.na(timeRecruit)) {
    ## calculate hazardRate
    ## needs to handle cases where numerics break down, e.g. timeRecruit==timeTotal and surviveEnd very close to zero
    hazardRate <- try({
      ## lower bound occurs with instantaneous recruitment timeRecruit==0
      hrMin <- -log(surviveEnd)/timeTotal
      ## find upper bound numerically
      hrMax <- hrMin
      while (exp(-timeTotal*hrMax)*(if (timeRecruit*hrMax>0) expm1(timeRecruit*hrMax)/(timeRecruit*hrMax) else 1) > surviveEnd) hrMax <- hrMax*2
      ## find value between lower and upper bounds, /2 and *2 ensures stable numerics when timeRecruit is close to zero
      uniroot(function(hR) return(exp(-timeTotal*hR)*(if (timeRecruit*hR>0) expm1(timeRecruit*hR)/(timeRecruit*hR) else 1) - surviveEnd), c(hrMin/2, hrMax*2))$root}, silent = TRUE)
    if (class(hazardRate) == "try-error") {
      return(c(hazardRate = NA,
               surviveEnd = surviveEnd,
               timeTotal = timeTotal,
               timeRecruit = timeRecruit))
    } else {
      return(c(hazardRate = hazardRate,
               surviveEnd = surviveEnd,
               timeTotal = timeTotal,
               timeRecruit = timeRecruit))
    }
  } else if (!is.na(hazardRate) && is.na(surviveEnd) & !is.na(timeTotal) & !is.na(timeRecruit)) {
    ## calculate surviveEnd, closed form
    return(c(hazardRate = hazardRate,
             surviveEnd = exp(-timeTotal*hazardRate)*(if (timeRecruit*hazardRate>0) expm1(timeRecruit*hazardRate)/(timeRecruit*hazardRate) else 1),
             timeTotal = timeTotal,
             timeRecruit = timeRecruit))
  } else if (!is.na(hazardRate) && !is.na(surviveEnd) & is.na(timeTotal) & !is.na(timeRecruit)) {
    ## calculate timeTotal, solve closed form
    timeTotal <- log((if (timeRecruit*hazardRate>0) expm1(timeRecruit*hazardRate)/(timeRecruit*hazardRate) else 1)/surviveEnd)/hazardRate
    if (timeTotal >= timeRecruit) {
      return(c(hazardRate = hazardRate,
               surviveEnd = surviveEnd,
               timeTotal = timeTotal,
               timeRecruit = timeRecruit))
    } else {
      return(c(hazardRate = hazardRate,
               surviveEnd = surviveEnd,
               timeTotal = NA,
               timeRecruit = timeRecruit))
    }
  } else if (!is.na(hazardRate) && !is.na(surviveEnd) & !is.na(timeTotal) & is.na(timeRecruit)) {
    ## calculate timeRecruit, numerics
    ## lower bound on surviveEnd when timeRecruit == 0
    sMin <- exp(-timeTotal*hazardRate)
    ## upper bound on surviveEnd when timeRecruit == timeTotal
    sMax <- exp(-timeTotal*hazardRate)*expm1(timeTotal*hazardRate)/(timeTotal*hazardRate)
    if (surviveEnd >= sMin && surviveEnd <= sMax) {
      return(c(hazardRate = hazardRate,
                  surviveEnd = surviveEnd,
                  timeTotal = timeTotal,
                  timeRecruit = uniroot(function(tR) return(exp(-timeTotal*hazardRate)*(if (tR*hazardRate>0) expm1(tR*hazardRate)/(tR*hazardRate) else 1) - surviveEnd), c(0, timeTotal))$root))
    } else {
      ## not possible for any permissible value for timeRecruit
      return(c(hazardRate = hazardRate,
                  surviveEnd = surviveEnd,
                  timeTotal = timeTotal,
                  timeRecruit = NA))
    }
  } else {
    ## should never get here, but default in case
    return(c(hazardRate = NA,
             surviveEnd = NA,
             timeTotal = NA,
             timeRecruit = NA))
  }
})

## needs better bounds and range checking!
## trial.params(NA, .01, 1, 1) fails

## prob to survive >T is exp(-hT)
## T==t2 where t2=timeTotal
##    => S = exp(-h t2)
##
## T uniform random betweed t1 and t2; t2>t1
## where t2=timeTotal, t1=timeTotal-timeRecruit, and t2-t1=timeRecruit
## S = int(exp(-hT), dT, from=t1, to=t2) / (t2-t1)
## change variable u=-hT, careful with signs
## == int(-1/h * exp(u), d(u), from = -h t1, to = -h t2) / (t2-t1)
## == int( 1/h * exp(u), d(u), from = -h t2, to = -h t1) / (t2-t1)
## == {exp(-h t1) - exp(-h t2)}/{h(t2-t1)}
## == exp(-h t2) * {exp(-h [t1-t2]) - 1}/{h(t2-t1)}
## == exp(-h t2) * {exp(+h [t2-t1]) - 1}/{h(t2-t1)}
##
## special case t1=0, t2=timeRecruit==timeTotal
## S = {1 - exp(-h t2)}/{h t2}

