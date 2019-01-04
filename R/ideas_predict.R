#' Ideas make
#' 
#' Convert GWAS summary stats into posterior probabilities that can be used by IDEAS predict. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @author Karl Guo \email{karl.x.guo@@gsk.com}
#' @export
#' @param analysis GWAS analysis id. 
#' @param chrom Chromosome to query
#' @param pos Position
#' @param pval_le [Default = 5e-8] Max pval for GWAS 
#' @param case_emac_ge [Default = 25]
#' @return  Table of GWAS loci with posteriors for \code{ideas_predict}
#' @import dplyr
ideas_make = function(analysis, pval_le = 5e-08, chrom = NULL, pos = NULL, case_emac_ge = 25) {
  print("Generating input ...")
  my_gwas <- gwas(analysis, case_emac_ge = 25, pval_thresh = pval_le, 
                  style = FALSE, gene_annotate = FALSE)
  # TODO - Fix chrom from numeric to string - remove NULL, check missing
  if (is.numeric(chrom) == TRUE & is.numeric(pos) == TRUE) {
    my_gwas = my_gwas %>% filter(chrom == chrom, pos_index == pos)
  }
  all_significant_regions <- 
    my_gwas %>% 
    group_by(chrom, pos_index, pval_index) %>% 
    do(regionplot.data(analysis, chrom = .$chrom, pos = .$pos_index, 
                       case_emac_ge = 25, style = "signals", maf_ge > 0.001)) %>% 
    ungroup() %>% 
    select(pos_index, cs_signal, freq, rsq, beta, pp_signal, pval, chrom, pos, ref, alt) %>% 
    filter(cs_signal == TRUE) %>% 
    filter(pp_signal > 0) %>% 
    distinct() %>% 
    # TODO - incorportate ref & alt
    select(-ref, -alt)
  
  all_significant_regions$pval = -log10(all_significant_regions$pval)
  # all_significant_regions$chrom = paste("chr", all_significant_regions$chrom, sep = "") 
  
  return(all_significant_regions)
}

#' ideas_predict
#' 
#' Predict IDEAS states from input GWAS data supplied by \code{ideas_make}
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @author Karl Guo \email{karl.x.guo@@gsk.com}
#' @export
#' @param .data
#' @param statepref
#' @param states
#' @param permute
#' @param vlp
#' @param states
#' @param inter [Default = 5]
#' @param alpha [Default = 0.99]
#' @param layer [Default = 10]
#' @param permute [Default = FALSE]
#' @param binsz [Default = 200]
#' @param pcut [Default = 6]
#' @param CI [Default = 0.95]
#' @param seed [Default = NULL]
#' @return TBD
ideas_predict <- function(.data, statepref = "/home/kg746906/IDEASstate/bin98.", 
                          states = NULL, permute = FALSE, inter = 5, binsz = 200, 
                          pcut = 6, CI = 0.95, seed = NULL) {
  BB = inter
  set.seed(seed)
  x = as.data.frame(.data)
  t = NULL
  for (i in unique(x[, 1])) {
    tt = which(x[, 1] == i & x[, 9] > 0 & x[, 6] > 0 & x[, 7] > pcut)
    if (length(tt) > 0) {
      o = tt[order(x[tt, 6], decreasing = T)]
      p = cumsum(x[o, 6])
      t = c(t, o[1:min(which(p > CI * sum(x[o, 6])))])
    }
  }
  futile.logger::flog.info(glue::glue("Number of cred sets: {length(unique(x[t, 1]))}"))
  futile.logger::flog.info(glue::glue("Number of snps across cred sets: {length(t)}"))
  if (length(t) < 10) {
    return(NULL)
  }
  rt = ideas_lasso(chrom  = as.matrix(x[t, 8]), 
                   pos    = cbind(as.numeric(x[t, 9]) - binsz, as.numeric(x[t, 9]) + binsz), 
                   vindex = x[t, 1], 
                   vprior = x[t, 6], 
                   vlp    = x[t, 7], 
                   statepref, 
                   states, 
                   permute = permute, 
                   inter = inter)
  rt$sel = t
  return(rt)
}

#' ideas_plot
#' 
#' Plot \code{ideas_predict}
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @author Karl Guo \email{karl.x.guo@@gsk.com} 
#' @export
#' @param .data \code{ideas_predict} output
#' @param fout [Default = NULL]
#' @param cutoff [Default = 0.001]
#' @param group [Default = FALSE]
#' @param relaxed [Default = NA]
#' @param legend [Default = TRUE]
#' @return TBD
ideas_plot <- function(.data, fout = NULL, cutoff = 0.001, group = FALSE, relaxed = NA, legend = TRUE) {
  if (length(fout) > 0) {
    pdf(paste(fout, ".pdf", sep = ""), width = 16, height = 8)
  }
  
  cell_info = suppressMessages(read_tsv("/home/kg746906/IDEASstate/epigenomeID_mod_key.txt"))
  a = t(.data$n)
  colnames(a) = c("first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eigth", "ninth", "tenth")
  results = cbind(cell_info, a)
  
  
  if (is.numeric(relaxed) == TRUE & relaxed != 1) {
    results$first = rowSums(results[, c("first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eigth", "ninth", 
                                        "tenth")[1:relaxed]])
    results$first = results$first/sum(results$first)
  }
  results = results[order(results$first), ]
  
  if (group == TRUE) {
    groups = unique(dplyr::select(results, Group, Color))
    first = c()
    for (i in 1:dim(groups)[1]) {
      sum = filter(results, Group == groups[i, 1])$first %>% sum()
      first = c(first, sum)
    }
    z = matrix(0, nrow = dim(groups)[1], ncol = 1)
    results = cbind.data.frame(z, z, z, z, groups, first)
    colnames(results) = c("a", "b", "c", "d", "Group", "Coloc", "first")
    results = results[order(results$first), ]
  }
  
  if (is.numeric(cutoff) == TRUE) {
    results = dplyr::filter(results, first > cutoff)
  }
  
  if (group == FALSE) {
    labels = results$EpigenomeMnemonic
    labels2 = unique(dplyr::select(results, Group, Color))[, 1]
    cc2 = unique(dplyr::select(results, Group, Color))[, 2]
  }
  if (group == TRUE) {
    labels = results$Group
  }
  
  cclr = col2rgb(as.matrix(results[, 6]))
  hclr = apply(cclr, 2, function(x) {
    rgb2hsv(x[1], x[2], x[3])
  })
  cc = c(apply(hclr, 2, function(x) {
    hsv(x[1], x[2], x[3])
  }))
  par(mar = c(7, 4, 3, 3))
  x = barplot(results$first, col = cc, cex.names = 0.6, ylab = "Proportion of genetic signal")
  if (legend == TRUE & group == FALSE) {
    legend("topleft", inset = 0.02, cex = 0.7, fill = cc2, labels2, ncol = 2)
  }
  text(x, labels = labels, srt = 45, par("usr")[3], adj = c(1, 1.1), xpd = TRUE, cex = 0.7)
  results_tbl = cbind.data.frame(labels, results$first)
  colnames(results_tbl) = c("cell", "proportion")
  return(results_tbl)
  if (length(fout) > 0) {
    dev.off()
  }
}

#' ideas_loci
#' 
#' IDEAS predict data for a locus. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @author Karl Guo \email{karl.x.guo@@gsk.com}
#' @export
#' @param ideas_make_obj
#' @param ideas_predict_obj
#' @return TBD
ideas_loci = function(ideas_make_obj, ideas_predict_obj) {
  tissues = suppressMessages(read_tsv("/home/kg746906/IDEASstate/epigenomeID_mod_key.txt"))
  ind = sort(ideas_predict_obj$sel)
  input_subset = ideas_make_obj[ind, ]
  scores = t(t(ideas_predict_obj$score) * ideas_predict_obj$cellscore)
  colnames(scores) = tissues$EpigenomeMnemonic
  scores = (input_subset$pp_signal * scores)
  SNP_scores = cbind.data.frame(input_subset, scores)
  locus_scores = c()
  for (i in 1:dim(unique(SNP_scores[, c(1, 8)]))[1]) {
    index = unique(SNP_scores[, c(1, 8)])[i, 1]
    chr = unique(SNP_scores[, c(1, 8)])[i, 2]
    slice = dplyr::filter(SNP_scores, pos_index == index)
    ordered = names(sort(colSums(slice[, 10:136]), decreasing = TRUE))
    rank_1 = ordered[1]
    rank_2 = ordered[2]
    rank_3 = ordered[3]
    rank_4 = ordered[4]
    rank_5 = ordered[5]
    rank_6 = ordered[6]
    rank_7 = ordered[7]
    rank_8 = ordered[8]
    rank_9 = ordered[9]
    rank_10 = ordered[10]
    x = cbind.data.frame(chr, index, rank_1, rank_2, rank_3, rank_4, rank_5, rank_6, rank_7, rank_8, rank_9, rank_10)
    locus_scores = rbind.data.frame(locus_scores, x)
  }
  colnames(locus_scores)[1:2] = c("chr", "pos")

  return(locus_scores)
}


#' ideas_nearest_state
#' 
#' Find the closest IDEAS state(s) to a chrom:pos. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @author Karl Guo \email{karl.x.guo@@gsk.com}
#' @param chrom Chromosome to query
#' @param pos Position
#' @param statepref
#' @param staten
#' @param states [Default = NULL]
#' @return TBD
ideas_nearest_state <- function(chrom, pos, statepref, staten, states = NULL) {
  rt = NULL
  rt$X = rt$P0 = NULL
  for (i in unique(chrom)) {
    print(i)
    t = which(chrom == i)
    if (length(t) == 0) {
      next;
    }
    
    vpst = as.integer((pos[t, 1] - 1)/200)
    vped = as.integer((pos[t, 2] - 1)/200)
    # TODO - update to SQL
    # I think I can delete this if statement and move right to query 
    if (length(states) > 0) {
      g <- states[which(states[, 2] == i), ]
    } else {
      # TODO - update to SQL
      g <- fread(paste(statepref, i, ".state", sep = ""))
    }
    # Select only cell type data by removing cols
    # select(-id, -chrom, -pos_start, -pos_end, -posclass) 
    gp <- as.integer(as.matrix(g[, 3])/200)
    
    g <- as.matrix(g[, 5:131])
    celln = ncol(g)
    # Tabulate for each col
    p0 = apply(g, 2, function(x) {
      tabulate(x + 1, nbins = staten)
    })
    
    mst = match(vpst, gp)
    med = match(vped, gp)
    tt = which(is.na(mst) == T)
    if (length(tt) > 0) {
      for (j in tt) {
        tttt = which(gp > vpst[j])
        if (length(tttt) == 0) {
          mst[j] = -1
        } else {
          mst[j] = min(tttt)
        }
      }
    }
    tt = which(is.na(med) == T)
    if (length(tt) > 0) {
      for (j in tt) {
        tttt = which(gp < vped[j])
        if (length(tttt) == 0) {
          med[j] = -1
        } else {
          med[j] = max(tttt)
        }
      }
    }
    tt = which(is.na(mst) == F & is.na(med) == F & mst <= med & mst >= 0 & med >= 0)
    A = array(0, dim = c(length(vpst), celln * staten))
    if (length(tt) > 0) {
      for (l in 0:max(med[tt] - mst[tt])) {
        ttt = tt[which(med[tt] >= mst[tt] + l)]
        apply(cbind(ttt, g[mst[ttt] + l, ]), 1, function(x) {
          a = x[-1] + 1 + (1:celln - 1) * staten
          A[x[1], a] <<- A[x[1], a] + 1
          return(NULL)
        })
      }
    }
    rt$X = rbind(rt$X, cbind(t, A))
    rt$P0 = rbind(rt$P0, c(p0))
  }
  return(rt)
}

#' ideas_sum_pp_cross_cell_types
#' 
#' Sum posterior probabilities across cell types
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @author Karl Guo \email{karl.x.guo@@gsk.com}
#' @param prior
#' @param prob
#' @param cellweight
#' @param snpindex
#' @param confidence [Default = 0.95]
#' @param cut [Default = 0]
#' @param layer [Default = 10]
#' @return TBD
ideas_sum_pp_cross_cell_types <- function(prior, prob, cellweight, 
                                          snpindex, confidence = 0.95, cut = 0, layer = 10) {
  celln = dim(cellweight)[2]
  snpn = dim(cellweight)[1]
  
  p = rep(0, snpn)
  p0 = p
  bb = NULL
  for (j in unique(snpindex)) {
    tt = which(snpindex == j)
    p0[tt] = prior[tt]/sum(prior[tt] + 1e-100)
    p[tt] = prob[tt]/sum(prob[tt] + 1e-100)
    zz0 = sort(p0[tt], decreasing = TRUE)
    zz = sort(p[tt], decreasing = T)
    ttt = which(cumsum(zz) >= confidence)[1]
    p[tt[p[tt] < zz[ttt]]] = 0
    bb = rbind(bb, c(length(tt), mean(cumsum(zz0)), mean(cumsum(zz))))
  }
  p[p < cut] = 0
  n = NULL
  rc = t(apply(cellweight, 1, function(x) {
    rank(x)
  }))
  for (i in 1:layer) {
    l = apply(rc, 1, function(x) {
      which(x > celln - i)[1]
    })
    tn = rep(0, celln)
    for (j in 1:dim(rc)[1]) {
      tn[l[j]] = tn[l[j]] + p[j]
      rc[j, l[j]] = 0
    }
    n = rbind(n, tn)
  }
  n = n/sum(p)
  
  rt = NULL
  rt$n = n
  rt$prob = p
  rt$prior = p0
  rt$auc = bb
  
  return(rt)
}

#' ideas_lasso
#' 
#' Using gglasso, identify the most likely cell type per locus. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @author Karl Guo \email{karl.x.guo@@gsk.com}
#' @param chrom
#' @param pos
#' @param vindex
#' @param vprior
#' @param vlp
#' @param statepref
#' @param states
#' @param usegenomebackground [Default = TRUE]
#' @param inter [Default = 5]
#' @param alpha [Default = 0.99]
#' @param layer [Default = 10]
#' @param permute [Default = FALSE]
#' @return TBD
ideas_lasso <- function(chrom, pos, vindex, vprior, vlp, statepref, states = NULL, 
                        usegenomebackground = TRUE, inter = 5, alpha = 0.99, 
                        layer = 10, permute = FALSE){
  BB = inter
  # temp hard code
  staten = 34 
  # Fix to:
  # para_tbl <- tbl(impala, "gene_gwas.ideas_para")
  # staten <- para_tbl %>% tally()
  # staten = length(readLines(paste(statepref, "para", sep = ""))) - 1
  L = nrow(pos)
  
  oindex = vindex
  chrom = as.matrix(chrom)
  if (ncol(pos) == 1) {
    pos = as.matrix(pos[, 1])
    pos = cbind(pos, pos)
  } else {
    pos = as.matrix(pos)
  }
  if (usegenomebackground) {
    oi = vindex
    ov = chrom
    op = pos
    ol = vlp
    for (B in 1:BB) {
      tpos = op
      for (i in unique(ov)) {
        tt = which(ov == i)
        rg = range(op[tt, ])
        rg[1] = rg[1] - 1e+06
        if (rg[1] < 1) 
          rg[1] = 1
        rg[2] = rg[2] + 1e+06
        tsz = tpos[tt, 2] - tpos[tt, 1]
        tpos[tt, 1] = as.integer(runif(length(tt), rg[1], rg[2]))
        tpos[tt, 2] = tpos[tt, 1] + tsz
      }
      chrom = c(chrom, ov)
      pos = rbind(pos, tpos)
      vindex = c(vindex, paste("null", B, "_", oi, sep = ""))
      vprior = c(vprior, rep(0, dim(tpos)[1]))
      vlp = c(vlp, ol)
    }
  }
  
  message("Preprocessing...")
  trt = ideas_nearest_state(chrom, pos, statepref, staten, states)
  x = trt$X
  p0 = trt$P0
  trt$X = trt$P0 = NULL
  
  sel = sort(x[, 1])
  sel1 = sel[which(substr(vindex[sel], 1, 4) != "null")]
  vindex = vindex[sel]
  vprior = vprior[sel]
  vlp = vlp[sel]
  x = x[order(x[, 1]), -1]
  n = apply(x[, 1:staten], 1, sum)
  weight = rep(1, length(vindex))
  for (i in unique(vindex)) {
    t = which(vindex == i)
    weight[t] = 1/length(t)
  }
  t = which(substr(vindex, 1, 4) == "null")
  ovindex = vindex
  vindex[t] = substr(vindex[t], 1, as.integer(regexpr("_", vindex[t])) - 1)
  
  x = x/(n + 1e-100)
  xs = apply(x, 2, function(x) {
    sum(x * vprior)
  })/mean(vprior)
  celln = dim(x)[2]/staten
  t = which(vprior <= 0)
  
  print("Predicting ...")
  if (usegenomebackground == FALSE && length(t) > 10) {
    p0 = apply(x[t, ], 2, sum)
  } else {
    p0 = apply(array(p0, dim = c(length(p0)/dim(x)[2], dim(x)[2])), 2, sum)
  }
  p0 = p0/sum(p0[1:staten])
  
  statefold = log2((xs + 1)/(p0 * sum(xs[1:staten]) + 1))
  score = t(apply(x, 1, function(x) {
    apply(array(x * statefold, dim = c(staten, celln)), 2, sum)
  }))
  if (permute) 
    score = score[sample(dim(score)[1]), ]
  
  gid = as.matrix(vindex)
  y = log((vprior + 0.001)/(1 - vprior + 0.001))
  xm = apply(score, 1, median)
  e = NULL
  for (B in 1:BB) {
    mm = which(substr(vindex, 1, 4) != "null")
    mm = c(mm, which(vindex == paste("null", B, sep = "")))
    t = 1:dim(score)[2]
    rt = cv.gglasso(cbind(xm[mm], score[mm, t]), y[mm], group = c(1, rep(1 + 1:celln, each = dim(score)[2]/celln)), pred.loss = "L2", 
                    loss = "ls", nfolds = 10)
    te = coef(rt, s = "lambda.min")[-1, 1]
    ee = rep(0, 1 + dim(score)[2])
    ee[1 + t] = te[-1]
    ee[1] = te[1]
    t = which(ee[-1] != 0)
    print(c(B, length(t)))
    e = rbind(e, ee)
    if (B == 1) {
      f = c(cbind(xm[mm], score[mm, ]) %*% ee)
    }
  }
  e = apply(e, 2, mean)
  
  p0 = p = rep(0, L)
  tf = c(cbind(xm, score) %*% e)
  e = e[-1]
  p[sel1] = (vprior * exp(tf)/(vprior * exp(tf) + 1 - vprior))[substr(vindex, 1, 4) != "null"]
  mm = which(substr(vindex, 1, 4) != "null")
  mm = c(mm, which(vindex == paste("null", 1, sep = "")))
  print(cor(f, y[mm]))
  p0[sel1] = vprior[which(substr(vindex, 1, 4) != "null")]
  
  rt = NULL
  rt$zscore = cor(f, y[mm]) * sqrt(length(mm))
  rt$statescore = statefold
  rt$p = p
  rt$p0 = p0
  rt$score = array(0, dim = c(L, dim(score)[2]))
  rt$score[sel1, ] = score[which(substr(vindex, 1, 4) != "null"), ]
  rt$score0 = score[which(substr(vindex, 1, 4) == "null"), ]
  rt$x = x
  
  e[e < 0] = 0
  rt$cellscore = e
  w = t(t(rt$score) * e)
  
  w = t(apply(w, 1, function(x) {
    apply(array(x, dim = c(length(x)/127, 127)), 2, sum)
  }))
  rt$cellweight = array(0, dim = c(L, dim(w)[2]))
  rt$cellweight[sel1, ] = w
  
  e = apply(array(e, dim = c(length(e)/127, 127)), 2, sum)
  ss = t(apply(x, 1, function(x) {
    return(array(x * statefold, dim = c(staten, celln)) %*% e)
  }))
  ss = ss[which(substr(vindex, 1, 4) != "null"), ]
  rt$stateweight = array(0, dim = c(L, staten))
  rt$stateweight[sel1, ] = ss
  
  tt = ideas_sum_pp_cross_cell_types(rt$p0, rt$p, rt$cellweight, oindex, confidence = 0.95, layer = layer)
  rt$n = tt$n
  rt$adjprob = tt$prob
  rt$adjprior = tt$prior
  rt$auc = tt$auc
  rt$statepref = statepref
  
  return(rt)
}



# TODO
# ideas_wrapper (make, predict, plot)