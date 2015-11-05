textgrid <- function(tmat, x0 = 0, y0 = 0, x1 = x0 + 1, y1 = y0 + 1, colsep = 3, rowsep = 1) {
  if (x1 < x0) {tmp <- x1; x1 <- x0; x0 <- tmp; rm(tmp)}
  if (y1 < y0) {tmp <- y1; y1 <- x0; y0 <- tmp; rm(tmp)}
  ## stop if x0==x1 or y0==y1

  ## coerce to matrix
  tmat <- cbind(rownames(tmat), apply(as.matrix(tmat), 2, as.character))
  tmat <- rbind(colnames(tmat), tmat)
  
  ## Remember previous setting
  oldcex <- par("cex")
  newcex <- oldcex
  
  ## Adaptively scale since actual strwidth drops in steps at discrete cex values
  while(TRUE) {
    xpos <- apply(tmat, 2, function(x) max(strwidth(x, units = "user")) +
                  strwidth("M", units = "user")*colsep)
    ypos <- apply(tmat, 1, function(x) max(strheight(x, units = "user")) +
                  strheight("M", units = "user")*rowsep)
    if (sum(xpos) <= (x1 - x0) && sum(ypos) <= (y1 - y0)) break
    newcex <- par("cex")/max(sum(xpos)/(x1 - x0), sum(ypos)/(y1 - y0))
    ## Check for stuck in infinite loop?
    if (newcex < 1e-4) stop("No font is small enough")
    par(cex = newcex)
  }

  xpos <- unname(x0 + c(0, cumsum(xpos)[-ncol(tmat)]))
  ypos <- unname(y1 - c(0, cumsum(ypos)[-nrow(tmat)]))
  for (idx in 1:nrow(tmat)) text(xpos, rep(ypos[idx], ncol(tmat)),
                                 tmat[idx, ], adj = c(0, 1))

  par(cex = oldcex)
  return(list(cex = newcex, xpos = xpos, ypos = ypos))
}

alphaize <- function(col, alpha = 0.5) {
  return(apply(col2rgb(col, alpha = FALSE),
               2,
               function(rgbvals) return(do.call(rgb, as.list(c(rgbvals/255, alpha = alpha))))))
}


