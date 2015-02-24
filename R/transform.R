## convenience function, convert +-Inf to NA
finitise <- function(x) return(ifelse(is.finite(x), x, NA))

## convenience function, convert <[=]0 to NA
positivise <- function(x, strict = FALSE) {
  if (strict) {
    return(ifelse(x > 0, x, NA))
  } else {
    return(ifelse(x >= 0, x, NA))
  }
}

## convenience function 
twopq <- function(p) return(2*p*(1 - p))

##
## zise() performs normal quantile transformation
##
zise <- function (x, only = NULL, by = NULL) {
  if (length(x) < 1) return(x)
  if (is.null(only)) {
    only <- which(!is.na(x))
  } else {
    stopifnot(length(only) == length(x))
    only <- which(only & !is.na(x))
  }
  zx <- rep(NA, length(x))
  if (length(only) >= 1) {
    if (is.null(by)) {
      zx[only] <- qnorm((rank(x[only]) - 0.5)/length(only))
    } else {
      stopifnot(length(by) == length(x))
      by <- factor(by)
      for (by1 in levels(by)) {
        onlyby <- intersect(only, which(by == by1))
        if (length(onlyby) >= 1) {
          zx[onlyby] <- qnorm((rank(x[onlyby]) - 0.5)/length(onlyby))
        }
      }
    }
  }
  return(zx)
}

zise.old <- function (x) {
  return(qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x))))
}

ile <- function (x, levels.out) {
  x <- as.double(x)
  ## if no levels.out specified, use three levels named 1,2,3
  if (missing(levels.out)) levels.out <- 1:3
  ## if levels.out is an integer of length 1, assume numbered levels
  if (identical(length(levels.out), 1L)) levels.out <- 1:as.integer(levels.out)
  n <- length(levels.out)
  xd <- ceiling(n * (rank(x, na.last = "keep") - 0.5)/sum(!is.na(x)))
  return(factor(levels.out[xd], levels.out))
}



##
## usage: safe(min, x), etc.
##
safe <- function(FUN, x, ...) {
  stopifnot(is.function(FUN))
  x <- as.vector(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  return(FUN(x, ...))
}

valuesof <- function(x, sep = ",", na.convert = NA) {
  x <- as.vector(x)
  x[is.na(x)] <- na.convert
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  paste(as.character(unique(x)), collapse = sep)
}
