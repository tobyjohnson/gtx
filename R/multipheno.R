multipheno <- function(z, cor.use = NULL, cor.method = "spearman") {
  z <- as.matrix(z)
  if (is.null(cor.use)) {
    rmat <- cor(z, method = cor.method, use = "complete")
  } else {
    stopifnot(length(cor.use) == nrow(z))
    rmat <- cor(z[which(cor.use), ], method = cor.method, use = "complete")
  }
  ir <- solve(rmat) # rmat usually small so okay to do this inefficiently
  zir <- tcrossprod(z, ir) # compute z *%* t(ir) == z %*% ir  because ir is symmetric
  ob <- apply(zir, 1, sum)/sqrt(sum(ir))
  pval.ob <- pchisq(ob^2, df = 1, lower.tail = FALSE)
  t2 <- sapply(1:nrow(z), function(idx) return(sum(zir[idx, , drop = TRUE] * z[idx, , drop = TRUE])))
  pval.t2 <- pchisq(t2, df = ncol(z), lower.tail = FALSE)
  return(list(rmatrix = rmat, OB = ob, OB.pval = pval.ob, T2 = t2, T2.pval = pval.t2))
}

multipheno.T2 <- function(z, cor.use = NULL, cor.method = "spearman") {
  z <- as.matrix(z)
  if (is.null(cor.use)) {
    rmat <- cor(z, method = cor.method, use = "complete")
  } else {
    stopifnot(length(cor.use) == nrow(z))
    rmat <- cor(z[which(cor.use), ], method = cor.method, use = "complete")
  }
  ir <- solve(rmat) # rmat usually small so okay to do this inefficiently
  zir <- tcrossprod(z, ir) # compute z *%* t(ir) == z %*% ir  because ir is symmetric
  t2 <- sapply(1:nrow(z), function(idx) return(sum(zir[idx, , drop = TRUE] * z[idx, , drop = TRUE])))
  pval <- pchisq(t2, df = ncol(z), lower.tail = FALSE)
  return(list(rmatrix = rmat, T2 = t2, pval = pval))
}

multipheno.OBrien <- function(z, cor.use = NULL, cor.method = "spearman") {
  z <- as.matrix(z)
  if (is.null(cor.use)) {
    rmat <- cor(z, method = cor.method, use = "complete")
  } else {
    stopifnot(length(cor.use) == nrow(z))
    rmat <- cor(z[which(cor.use), ], method = cor.method, use = "complete")
  }
  ## original implemention, was not exactly same weighting as O'Brien for >2 phenotypes
  ## ob <- apply(z, 1, sum)/sqrt(sum(rmat))
  ## pval <- pchisq(ob^2, df = 1, lower.tail = FALSE)
  ## new method, most powerful against alternative with E(Z1)==E(Z2)==...
  ir <- solve(rmat) # rmat usually small so okay to do this inefficiently
  zir <- tcrossprod(z, ir) # compute z *%* t(ir) == z %*% ir  because ir is symmetric
  ob <- apply(zir, 1, sum)/sqrt(sum(ir))
  pval <- pchisq(ob^2, df = 1, lower.tail = FALSE)
  return(list(rmatrix = rmat, OB = ob, pval = pval))
}
  
