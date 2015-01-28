## for a series of intervals that may overlap or contain each other,
## compute the minimum size of a cover.
## E.g. applications time intervals, physical bp intervals
mincover <- function(x.begin, x.end, trunc.begin = NULL, trunc.end = NULL) {
  stopifnot(length(x.begin) == length(x.end))
  if (!is.null(trunc.begin)) stopifnot(length(trunc.begin) == 1)
  if (!is.null(trunc.end)) stopifnot(length(trunc.end) == 1)
  if (length(x.begin) == 0) return(0)
  oo <- order(x.begin)
  x.begin <- as.integer(x.begin[oo])
  x.end <- as.integer(x.end[oo])
  if (length(x.begin) > 1) {
    x.begin[2:length(x.begin)] <- sapply(2:length(x.begin), function(ii) max(c(x.begin[ii], x.end[1:(ii-1)] + 1)))
    ## ii-th interval chopped to guarantee not to overlap intervals 1..(ii-1) => no overlap
    ## chopped parts guaranteed to be covered by other intervals => no loss of cover
    ## chopping may make some intervals have begin > end, excluded from sum by pmax below
  }
  if (!is.null(trunc.begin)) x.begin <- pmax(x.begin, as.integer(trunc.begin))
  if (!is.null(trunc.end)) x.end <- pmin(x.end, as.integer(trunc.end))
  return(sum(pmax(0, x.end - x.begin + 1)))
}
