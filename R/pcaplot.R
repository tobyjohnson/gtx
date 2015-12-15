pcaplot <- function(pcs, f) {
  ## Coerce arguments to correct types
  pcs <- as.matrix(pcs)
  f <- as.factor(f)
  stopifnot(all(c("PC1", "PC2", "PC3") %in% colnames(pcs)))
  stopifnot(all.equal(length(f), nrow(pcs)))
  ## Exclude any subjects with missing PCs or descriptive factor
  inc1 <- apply(!is.na(pcs), 1, all) & !is.na(f)
  pcs <- pcs[inc1, , drop = FALSE]
  f <- f[inc1]
  ## Sort levels of f by count
  f <- factor(f, names(rev(sort(table(f)))))
  ## Find a (deterministic) order so that smallest groups tend to be plotted on top
  otmp <- rep(NA, length(f))
  for (f1 in levels(f)) otmp[f == f1] <- max(table(f)) - 1:sum(f == f1)
  o <- order(otmp, as.integer(f)) 
  nlev <- length(levels(f))
  colset <- rainbow(nlev, start = 0, end = 0.7, v = 0.7, alpha = 0.3)
  ## if rgb does not take alpha argument, use: colset <- alphaize(rainbow(nlev, start = 0, end = 0.7, v = 0.7), alpha = 0.5)
  pchset <- rep(21:25, length.out = nlev) # transparency filled symbols
  ## if open symbols are preferred, use: pchset <- rep(c(1, 0, 5, 2, 6), length.out = nlev) # open symbols

  scs <- split.screen(c(2, 2)) # Must use split.screen, not par(mfrow) or other mechanisms inside pipeplot

  screen(scs[1]); par(mar = c(2, 2, 0, 0) + 0.1)
  plot(pcs[o, "PC2"], pcs[o, "PC1"],
       pch = pchset[f][o], col = colset[f][o], bg = colset[f][o], 
       xaxt = "n", yaxt = "n", ann = FALSE)
  mtext("PC2", 1, 0); mtext("PC1", 2, 0)
  screen(scs[2]); par(mar = c(2, 2, 0, 0) + 0.1)
  plot(pcs[o, "PC3"], pcs[o, "PC1"],
       pch = pchset[f][o], col = colset[f][o], bg = colset[f][o], 
       xaxt = "n", yaxt = "n", ann = FALSE)
  mtext("PC3", 1, 0); mtext("PC1", 2, 0)
  screen(scs[3]); par(mar = c(2, 2, 0, 0) + 0.1)
  plot.new()
  plot.window(c(0, 1), c(0, 1))
  ## FIXME replace with:
  legendgrid("top",
             legend = matrix(paste(levels(f), " (n=", table(f), ")", sep = ""), ncol = 1),
             include.colnames = FALSE, 
             pch = pchset, col = colset, bg = colset)
#  tmp <- textgrid(matrix(paste(levels(f), " (n=", table(f), ")", sep = ""), ncol = 1),
#                  x0 = strwidth("M", units = "user", cex = 1),
#                  y0 = 0, x1 = 1, y1 = 1)
#  points(rep(0, nlev), tmp$ypos - 0.5*strheight("M", units = "user", cex = tmp$cex),
#         cex = tmp$cex, pch = pchset, col = colset, bg = colset)
  mtext("Legend", 3, 0)
  screen(scs[4]); par(mar = c(2, 2, 0, 0) + 0.1)
  plot(pcs[o, "PC3"], pcs[o, "PC2"],
       pch = pchset[f][o], col = colset[f][o], bg = colset[f][o], 
       xaxt = "n", yaxt = "n", ann = FALSE)
  mtext("PC3", 1, 0); mtext("PC2", 2, 0)
  close.screen(scs)            

  return(invisible(NULL))
}
