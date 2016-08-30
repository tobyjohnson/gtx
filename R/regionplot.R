#library for region plot
#Author: Li Li, modified based on Toby's code



#draw recomination rate 
recomb.draw <- function(chr, from.bp, to.bp, ylo, yhi, recomb){
  if (is.integer(chr)) chr <- as.character(chr)
  if(chr == "23") chr <- "X"
  chr <- ifelse(substr(chr, 1, 3) == "chr", chr, paste("chr", chr, sep = ""))
  #recomb <- read.table(paste("genetic_map_chr", chr, ".txt", sep=""), header=T)
  #load recombination info
  recomb.chr <- recomb[[chr]][-1]
  idx.start<- which(recomb.chr[,1] < from.bp)
  idx.end<- which(recomb.chr[,1] > to.bp)
  keep.recomb <-recomb.chr[idx.start[length(idx.start)]: idx.end[1],]
  big.range <- yhi -ylo
  max.curr<- max(keep.recomb[,2], na.rm = T)
  #print (paste("Max recomb rate in the region is", max.curr))
  ticks.curr <- pretty(c(0, max.curr))
  max.curr<- ticks.curr[length(ticks.curr )]
  lines(keep.recomb[,1], ylo + ( ( keep.recomb[,2] / max.curr ) * big.range), type="l", 
        col="lightblue", lwd=1)
  axis(4, at=ylo + big.range /max.curr * ticks.curr, labels=ticks.curr, las=1, hadj = 0)
  mtext("Recomb rate (cM/Mb)", side=4, at=(ylo+big.range/2), line=2.5) 
}



################################################################################
#gwas: gwas results with columns for MARKER, chr, POS (bp), qcflag, genoflag, pvalues
#chr pos: hit snp
#flanking: flanking rang (kb) for plotting
#col.pvals: column names for p values to be plotted by different shapes (21:25, 3,4,7:14) 
#plabel: label for p value column
#col.r2: column name for r2 to index SNP. Background colored by LD to index SNP (white to red)
################################################################################
regionplot.multiP <- function(gwas, chr, pos, flanking = 250, col.r2=NULL,
                              col.pvals = c("P1","P2"), plabel = NULL, gencode = gencode, recomb=recomb,
                              main = "", main.sub = NULL, miny = 8, crity = -log10(5e-08)) {
  to.bp <- min(pos+flanking *1000, max(gwas$POS, na.rm=T))
  from.bp <- max(pos - flanking *1000, min(gwas$POS, na.rm = T))
  stopifnot(to.bp > from.bp)
  stopifnot(all(col.pvals %in% names(gwas)))
  sel <- gwas$CHROM == chr & gwas$POS >= from.bp & gwas$POS <= to.bp
  cat(sum(sel), "SNPs selected\n")
  pos.sel<- gwas$POS[sel]
  minP <- 1e-50 #minimum value for p value
  for( c in col.pvals) gwas[[c]][!is.na(gwas[[c]]) & gwas[[c]] < minP] <- minP
  lp <- -log10(as.matrix(subset(gwas, sel, col.pvals))) # lp = -log10(p)
  if (is.null(crity)) crity <- -log10(0.05/nrow(lp))
  if (is.null(plabel)) plabel <- col.pvals
  r2.sel <- rep(1, sum(sel))
  if(!is.null(col.r2)) {
    r2.sel<- gwas[[col.r2]][sel]
    r2.sel[is.na(r2.sel)]<- 0
    r2.sel <- pmin(1, pmax(0, 1-r2.sel))
  }
 
 
  plot.new()
  yhi<- max(c(miny, 1.05 * max(lp, na.rm = TRUE)))
  y.tick<- pretty(c(0, yhi))
  yhi <- y.tick[length(y.tick)]
  scale  <- yhi/8
  plot.window(c(from.bp, to.bp), c(-6 *scale,  yhi))
  
  abline(h = seq(0, 7, 1), col = "lightgrey", lty = "dotted")
  abline(h = crity, col = "red", lty = "dashed")
  ppch <- c(21:25, 3,4,7:14) #rep(21:25, 5)
  index <- pos.sel == pos 
  for (idx in 1:ncol(lp)) {
    points(pos.sel[-index], lp[-index , idx], pch = ppch[idx], col = "black",  bg = rgb(1, r2.sel[-index], r2.sel[-index]), cex = 1)
    if(any(index)) #index snp
      points(pos, lp[index, idx],  pch = ppch[idx], col = rgb(0, 0, 1),  bg = rgb(0, 0, 1), cex = 1.3)
  }   
  legend("topright",   ncol = 2,   bty = "n",  pch = ppch[1:ncol(lp)], 
           text.col =  "black", cex = 0.8,  legend = plabel)
  for(i in 0:10){
    polygon(from.bp + (c(0, 1, 1, 0)+i) *(to.bp-from.bp)/80, yhi- c(0.5, 0.5, 1, 1)*scale, 
            border = "grey", col = rgb(1, 1-i/10, 1-i/10)) 
  }
  text(from.bp, yhi - 0 *scale , expression(paste("LD Map Type:","r"^"2")), adj = 0, cex = 0.8)
  text(from.bp+seq(0.5, 10.5, 2) * (to.bp-from.bp)/80, rep(yhi, 6)- 1.5*scale, 
       labels =c(0, 0.2, 0.4, 0.6, 0.8,1),  cex = 0.7)
  
  xm <- pretty(c(from.bp, to.bp)*1e-6) # x marks
  axis(1, at = xm*1e6, labels = xm)
  axis(2, at = y.tick, las = 1)
  mtext(expression(-log[10](italic(P))), side=2, at=yhi/2, line=2.5) 
  recomb.draw(chr, from.bp, to.bp, -1*scale, yhi, recomb)                        
  gene.draw(paste("chr", chr, sep = ""), from.bp, to.bp, gencode, 
            yhi = -2*scale, ylo = -6*scale,exony = 0.05*scale, genecex = 0.7)
  mtext(paste("Chromosome ", chr, " genomic position (Mb)", sep = ""), side = 1, line = 2.5)
  title(main = main)  # ylab = expression(-log[10](italic(P))),
  if(!is.null(main.sub))
    title(main = main.sub, cex.main = 1, col.main = "blue", line = 0.5)
  box()
  return(invisible(NULL))
}


