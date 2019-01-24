#' Ideas make
#' 
#' Convert GWAS summary stats into posterior probabilities that can be used by IDEAS predict. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @author Karl Guo \email{karl.x.guo@@gsk.com}
#' @export
#' @family ideas_predict
#' @param analysis GWAS analysis id. 
#' @param chrom Chromosome to query
#' @param pos Position
#' @param pval_le [5e-8] Max pval for GWAS 
#' @param case_emac_ge [25] 
#' @return Table of GWAS loci with posteriors for \code{\link{ideas_predict}}
#' @import dplyr
#' @import purrr
ideas_make <- function(analysis, pval_le = 5e-08, chrom, pos, case_emac_ge = 25,
                      dbc = getOption("gtx.dbConnection", NULL)) {
  gtxdbcheck(dbc)
  gtx_info("ideas_make | Calculating GWAS top hits for: {analysis}.")
  # 
  if(!is_character(analysis)){
    gtx_error("ideas_make | analysis must be of type 'character'.");
    stop();
  }
  safely_gwas <- purrr::safely(gtx::gwas)
  exec <- safely_gwas(analysis, case_emac_ge = 25, pval_thresh = pval_le, 
                      style = FALSE, gene_annotate = FALSE)
  
  if(is_null(exec$result)){
    gtx_warn("ideas_make | gwas() returned zero results for analysis: {analysis} at pval_thresh: {pval_le}.")
    return(NULL);
  } else {
    my_gwas <- exec$result
  }
  
  if (!missing(chrom) & !missing(pos)) {
    gtx_debug("ideas_make | Filtering top hits for input chrom & pos.")
    if(!is_character(chrom)){
      gtx_error("ideas_make | param 'chrom' is not a character.")
    }
    if(!is.numeric(pos)){
      gtx_error("ideas_make | param 'pos' is not numeric.")
    }
    my_gwas = my_gwas %>% dplyr::filter(chrom == chrom, pos_index == pos)
    if(nrow(my_gwas) == 0){
      gtx_error("ideas_make | No GWAS results after filtering for chrom:{chrom} & pos:{pos}.")
    }
  }
  
  gtx_info("ideas_make | Calculating cred sets for each GWAS top hit.")
  all_significant_regions <- 
    my_gwas %>% 
    dplyr::group_by(chrom, pos_index, pval_index) %>% 
    dplyr::do(regionplot.data(analysis, chrom = .$chrom, pos = .$pos_index, 
                              case_emac_ge = 25, style = "signal", maf_ge = 0.001)) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(pos_index, cs_signal, freq, rsq, beta, pp_signal, pval, chrom, pos, ref, alt) %>% 
    dplyr::filter(cs_signal == TRUE) %>% 
    dplyr::filter(pp_signal > 0) %>% 
    dplyr::distinct() %>% 
    # TODO - incorportate ref & alt
    dplyr::select(-ref, -alt)
  
  if(nrow(all_significant_regions) == 0){
    gtx_error("ideas_make | No cred sets returned for analysis:{analysis}.")
  }
  
  all_significant_regions$pval = -log10(all_significant_regions$pval)
  
  return(all_significant_regions)
}

#' ideas_predict
#' 
#' Predict IDEAS states from input GWAS data supplied by \code{\link{ideas_make}}
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @author Karl Guo \email{karl.x.guo@@gsk.com}
#' @export
#' @family ideas_predict
#' @param .data \code{\link{ideas_make}} input data
#' @param db = [`gtx::config_db()`] Database contained IDEAS states data.
#' @param permute [FALSE] TODO
#' @param states [NULL] TODO
#' @param inter [5] TODO
#' @param alpha [0.99] TODO
#' @param layer [10] TODO
#' @param permute [FALSE] TODO
#' @param binsz [200] Padding around each cred set SNP.
#' @param pcut [6] TODO
#' @param CI [0.95] Cred set confidence interval.
#' @param seed [NULL] Set seed for RNG number generator.
#' @param states Select subset of states. Not currently used as of 2019-01-11.
ideas_predict <- function(.data,  states_data = NULL,
                          states = NULL, permute = FALSE, inter = 5, binsz = 200, 
                          pcut = 6, CI = 0.95, seed = NULL, db = config_db()) {
  if(missing(.data)){
    gtx_error("ideas_predict | no input .data.");
    stop();
  } else if(is_null(.data)){
    gtx_warn("ideas_predict | .data is NULL.");
    return(NULL);
  } else {
    x = as.data.frame(.data);
  }
  
  BB = inter
  set.seed(seed)
  t = NULL
  for (i in unique(x[, 1])) {
    tt = which(x[, 1] == i & x[, 9] > 0 & x[, 6] > 0 & x[, 7] > pcut)
    if (length(tt) > 0) {
      o = tt[order(x[tt, 6], decreasing = T)]
      p = cumsum(x[o, 6])
      t = c(t, o[1:min(which(p > CI * sum(x[o, 6])))])
    }
  }
  gtx_info("ideas_predict | Number of cred sets: {length(unique(x[t, 1]))}")
  gtx_info("ideas_predict | Number of snps across cred sets: {length(t)}")
  
  gtx_debug("ideas_predict | Starting ideas_lasso . . .");
  rt = ideas_lasso(chrom  = as.matrix(x[t, 8]), 
                   pos    = cbind(as.numeric(x[t, 9]) - binsz, as.numeric(x[t, 9]) + binsz), 
                   vindex = x[t, 1], 
                   vprior = x[t, 6], 
                   vlp    = x[t, 7], 
                   states_data = states_data,
                   permute = permute, 
                   inter = inter,
                   states = states, 
                   db = db)
  rt$sel = t
  gtx_info("ideas_predict | Prediction complete.")
  return(rt)
}

#' ideas_loci_summarize
#' 
#' IDEAS predict data for a locus. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @author Karl Guo \email{karl.x.guo@@gsk.com}
#' @export
#' @family ideas_predict
#' @param ideas_make_dat \code{\link{ideas_make}} object
#' @param ideas_predict_dat \code{\link{ideas_predict}} object
#' @param impala [getOption("gtx.impala", NULL)] Implyr impala connection
#' @param db [config_db()] Database to use for queries.
ideas_loci_summarize = function(ideas_make_dat, ideas_predict_dat) {
  gtx_debug("ideas_loci_summarize | Staring ideas_loci_summarize.")
  # Verify input
  if(missing(ideas_make_dat)){
    gtx_error("ideas_loci_summarize | param 'ideas_make_dat' must be supplied.")
    stop();
  }
  if(missing(ideas_predict_dat)){
    gtx_error("ideas_loci_summarize | param 'ideas_predict_dat' must be supplied.")
    stop();
  }
  if(is_null(ideas_make_dat) | is_null(ideas_predict_dat)){
    if(is_null(ideas_make_dat)){
      gtx_warn("ideas_loci_summarize | ideas_make_dat is NULL. Skipping.");
      return(NULL);
    } else {
      gtx_warn("ideas_loci_summarize | ideas_predict_dat is NULL. Skipping.");
      return(NULL);
    }
  }
  
  tissues <- ideas_predict_dat$ideas_metadata
  ind = sort(ideas_predict_dat$sel)
  input_subset = ideas_make_dat[ind, ]
  scores = t(t(ideas_predict_dat$score) * ideas_predict_dat$cellscore)
  colnames(scores) = tissues$epigenome_mnemonic
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
  colnames(locus_scores)[1:2] = c("chrom", "pos")

  gtx_debug("ideas_loci_summarize | complete.")
  
  return(locus_scores)
}


#' ideas_get_states
#' 
#' Find the closest IDEAS state(s) to a chrom:pos. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @author Karl Guo \email{karl.x.guo@@gsk.com}
#' @param chrom Chromosome to query
#' @param pos Position
#' @param staten Number of states. 
#' @param states [NULL]
#' @param impala [getOption("gtx.impala", NULL)] Implyr impala connection
#' @param db [config_db()] Database to use for queries.
ideas_get_states <- function(chrom, pos, staten, states = NULL, states_data,
                             impala = getOption("gtx.impala", NULL), db = config_db()) {
  impala <- validate_impala(impala = impala);
  # If we didn't pass states_data, load all the data based on the input chrom
  if(missing(states_data)){
    gtx_debug("ideas_get_states | Querying states from: {db}.ideas_states")
    sql_statement <- 
      glue::glue_collapse(
        c(glue::glue("SELECT * FROM {db}.ideas_states WHERE"), 
          glue::glue_collapse(glue::glue("chrom = \"{unique(chrom)}\""), sep = " OR ")), 
        sep = " ")
    
    safely_get_query <- purrr::safely(implyr::dbGetQuery)
    exec <- safely_get_query(impala, sql_statement)
    if (!is_null(exec$error)){
      gtx_error("ideas_get_states | unable to query states (1), error:\n{exec$error}");
      gtx_error("Failed SQL query:\n {sql_statement}");
      stop()
    } else {
      gtx_debug("ideas_get_states | states collected.");
      states_data <- exec$result
    }
  } else if (is_null(states_data)){
    gtx_debug("ideas_get_states | Querying states from: {db}.ideas_states")
    sql_statement <- 
      glue::glue_collapse(
        c(glue::glue("SELECT * FROM {db}.ideas_states WHERE"), 
          glue::glue_collapse(glue::glue("chrom = \"{unique(chrom)}\""), sep = " OR ")), 
        sep = " ")
    
    safely_get_query <- purrr::safely(implyr::dbGetQuery)
    exec <- safely_get_query(impala, sql_statement)
    if (!is_null(exec$error)){
      gtx_error("ideas_get_states | unable to query states (2), error:\n{exec$error}");
      gtx_error("Failed SQL query:\n {sql_statement}");
      stop();
    } else {
      gtx_debug("ideas_get_states | states collected.");
      states_data <- exec$result
    }
  }
  else {
    gtx_debug("ideas_get_states | Using states_data from 'states_data' parameter.");
  }
  
  rt = NULL
  rt$X = rt$P0 = NULL
  for (i in unique(chrom)) {
    gtx_debug("ideas_get_states | Harmonizing states across chr:{i}")
    t = which(chrom == i)
    if (length(t) == 0) {
      next;
    }
    
    vpst = as.integer((pos[t, 1] - 1)/200)
    vped = as.integer((pos[t, 2] - 1)/200)
    g <- states_data %>% dplyr::filter(chrom == i) %>% arrange(pos_start)
    gp <- as.integer(as.matrix(g$pos_start)/200)
    # Remove extra cols, leaving just cell type + data
    g <- g %>% dplyr::select(-id, -chrom, -pos_start, -pos_end, -posclass) %>% as.matrix() 
    celln = ncol(g)
    # Tabulate for each col
    p0 = apply(g, 2, function(x) {
      tabulate(x + 1, nbins = staten)
    })
    
    gtx_debug("ideas_get_states | Cross reference GWAS and states.")
    mst = match(vpst, gp)
    med = match(vped, gp)
    tt = which(is_null(mst) == T)
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
    tt = which(is_null(med) == T)
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
    tt = which(is_null(mst) == F & is_null(med) == F & mst <= med & mst >= 0 & med >= 0)
    A = array(0, dim = c(length(vpst), celln * staten))
    if (length(tt) > 0) {
      for (l in 0:max(med[tt] - mst[tt])) {
        ttt = tt[which(med[tt] >= mst[tt] + l)]
        apply(cbind(ttt, g[mst[ttt] + l, ]), 1, function(x) {
          a = x[-1] + 1 + (1:celln - 1) * staten
          # WTF is this <<-, seems dangerous ... @KBS
          A[x[1], a] <<- A[x[1], a] + 1
          return(NULL)
        })
      }
    }
    rt$X = rbind(rt$X, cbind(t, A))
    rt$P0 = rbind(rt$P0, c(p0))
  }
  gtx_debug("ideas_get_states | Loading & harmonizing all chromosome states complete.")
  return(rt)
}

#' ideas_sum_pp_cross_cell_types
#' 
#' Sum posterior probabilities across cell types
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @author Karl Guo \email{karl.x.guo@@gsk.com}
#' @param prior TBD
#' @param prob TBD
#' @param cellweight TBD
#' @param snpindex TBD
#' @param confidence [0.95]
#' @param cut [0]
#' @param layer [10]
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
#' @param chrom Chromosomes 
#' @param pos Positions
#' @param vindex TBD
#' @param vprior TBD
#' @param vlp TBD
#' @param states TBD
#' @param use_genome_background [TRUE]
#' @param inter [5]
#' @param alpha [0.99]
#' @param layer [10]
#' @param permute [FALSE]
#' @param impala [getOption("gtx.impala", NULL)] Implyr impala connection
#' @import gglasso
#' @import dplyr
#' @return TBD
ideas_lasso <- function(chrom, pos, vindex, vprior, vlp, states_data = NULL, 
                        use_genome_background = TRUE, inter = 5, alpha = 0.99, 
                        layer = 10, permute = FALSE, states = NULL, 
                        impala = getOption("gtx.impala", NULL), db = config_db()){
  
  gtx_debug("ideas_lasso | Starting ideas_lasso.");
  gtx_debug("ideas_lasso | Validate inputs - starting.");
  if(any(map_lgl(c(chrom, pos, vindex, vprior, vlp), missing))){
    gtx_error("ideas_lasso | Validate inputs - failure. Check: chrom, pos, vindex, vprior, vlp");
    stop("ideas_lasso failure.");
  }
  gtx_debug("ideas_lasso | Validate inputs - success.");  
  
  BB = inter
  
  impala <- validate_impala(impala = impala)
  para_tbl <- dplyr::tbl(impala, glue::glue("{db}.ideas_para"))
  staten <- para_tbl %>% dplyr::tally() %>% dplyr::pull() %>% as.numeric()
  gtx_debug("ideas_lasso | Number of states = {staten}.");
  
  gtx_debug("ideas_lasso | Querying ideas_metadata.");
  ideas_metadata <- implyr::dbGetQuery(impala, glue::glue("SELECT * FROM {db}.ideas_metadata"));
  # Hard code number of rows in metadata temporarily for debug
  if(nrow(ideas_metadata) == 127){
    gtx_debug("ideas_lasso | ideas_metadata collected.");  
  }
  
  L = nrow(pos)
  oindex = vindex
  chrom = as.matrix(chrom)
  if (ncol(pos) == 1) {
    pos = as.matrix(pos[, 1])
    pos = cbind(pos, pos)
  } else {
    pos = as.matrix(pos)
  }
  
  if (isTRUE(use_genome_background)) {
    gtx_debug("ideas_lasso | Create genome background.");
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
    gtx_debug("ideas_lasso | Create genome background - complete.");
  }
  
  gtx_info("ideas_predict | Preprocessing states.");
  gtx_debug("ideas_predict | ideas_lasso | Preprocessing with ideas_get_states()");
  trt = ideas_get_states(chrom  = chrom, 
                         pos    = pos, 
                         staten = staten, 
                         states = states, 
                         states_data = states_data)
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
  
  gtx_info("ideas_lasso | Predicting states")
  if (use_genome_background == FALSE && length(t) > 10) {
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
    rt = gglasso::cv.gglasso(cbind(xm[mm], score[mm, t]), 
                    y[mm], 
                    group = c(1, rep(1 + 1:celln, each = dim(score)[2]/celln)), 
                    pred.loss = "L2", 
                    loss = "ls", 
                    nfolds = 10)
    te = coef(rt, s = "lambda.min")[-1, 1]
    ee = rep(0, 1 + dim(score)[2])
    ee[1 + t] = te[-1]
    ee[1] = te[1]
    t = which(ee[-1] != 0)
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
  
  tt = ideas_sum_pp_cross_cell_types(rt$p0, 
                                     rt$p, 
                                     rt$cellweight, 
                                     oindex, 
                                     confidence = 0.95, 
                                     layer = layer)
  rt$n = tt$n
  rt$adjprob = tt$prob
  rt$adjprior = tt$prior
  rt$auc = tt$auc
  rt$ideas_metadata = ideas_metadata
  
  gtx_debug("ideas_lasso | complete.")
  
  return(rt)
}

#' ideas_preload_states
#' 
#' Preload all IDEA states across "all" chrom. This is ideal for when running ideas_predict across multiple GWAS;
#' however, it would be slower for single GWAS.
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param impala [getOption("gtx.impala", NULL)] Implyr impala connection
#' @param db [config_db()] Database to use for queries.
#' @param chrom [c(1:22, "X", "Y")] Specify which chromosomes to load. 
#' @return states_data to be used with \code{\link{ideas_predict}}
#' @export
#' @family ideas_predict
#' @import implyr
#' @import dplyr
#' @import purrr
ideas_preload_states <- function(impala = getOption("gtx.impala", NULL), db = config_db(), 
                                 chrom = c(1:22, "X", "Y")) {
  impala <- validate_impala(impala = impala)
  # Build SQL statement to load all chrom data
  gtx_debug("ideas_preload_states | Querying states from: {db}.ideas_states")
  sql_statement <- 
    glue::glue_collapse(
      c(glue::glue("SELECT * FROM {db}.ideas_states WHERE"), 
        glue::glue_collapse(glue::glue("chrom = \"{unique(chrom)}\""), sep = " OR ")), 
      sep = " ")
  # Safely query data
  safely_get_query <- purrr::safely(implyr::dbGetQuery)
  states_data <- safely_get_query(impala, sql_statement)
  # Check if we had an error
  if (!is.null(states_data$error)){
    gtx_error("ideas_preload_states | unable to query states (0), error:\n{states_data$error}")
    stop()
  } else {
    # Without an error, keep the results to return
    states_data <- states_data$result %>% arrange(chrom, pos_start)
  }
  gtx_debug("ideas_preload_states | Querying complete.")
  return(states_data)
} 

#' ideas_preload_states
#' 
#' Preload all IDEA states across "all" chrom. This is ideal for when running ideas_predict across multiple GWAS;
#' however, it would be slower for single GWAS.
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param analysis A string or vector of analysis ids.
#' @param ht_load [FALSE] TRUE = Load the IDEAS states for high throughput.
#' @param states_data Pass \code{\link{ideas_preload_states}} data instead of loading it. 
#' @param relaxed [3] Cell types count toward total if they are in the top # of 'relaxed' cell types. e.g., 3 = The top 3 cell types count.
#' @param group [FALSE] TRUE = group cell types together.
#' @param dbc [getOption("gtx.dbConnection", NULL)]
#' @param impala [getOption("gtx.impala", NULL)]
#' @param db [`gtx::config_db()`]
#' @export
#' @family ideas_predict
#' @import dplyr
#' @import purrr
ideas_wrapper <- function(analysis, ht_load = FALSE, states_data,
                          relaxed = 3, group = FALSE,
                          impala = getOption("gtx.impala", NULL),
                          dbc    = getOption("gtx.dbConnection", NULL), 
                          db     = config_db()){
  gtxdbcheck(dbc);
  # Confirm we have input
  if(missing(analysis)){
    gtx_error("ideas_wrapper | parameter 'analysis' is required but is missing.");
    stop();
  }
  
  # Check/do if we need to do high throughput loading states loading. 
  if(missing(states_data) & ht_load == TRUE){
    states_data <- ideas_preload_states(db = db);
  } else if(missing(states_data) & ht_load == FALSE){
    states_data = NULL;
  } 
  
  # Setup input tibble
  input <- dplyr::tibble("analysis" = analysis);
  
  # Runs ideas predict for each analysis id. 
  ret <- input %>%
    dplyr::group_by(analysis) %>%
    dplyr::mutate(ideas_make_dat    = purrr::map(analysis,          ideas_make)) %>%
    dplyr::mutate(ideas_predict_dat = purrr::map(ideas_make_dat,    ideas_predict, states_data = states_data, db = db)) %>%
    dplyr::mutate(ideas_gwas_dat    = purrr::map(ideas_predict_dat, ideas_gwas_summarize, group = group, relaxed = relaxed)) %>%
    dplyr::mutate(ideas_loci_dat    = purrr::map2(ideas_make_dat,   ideas_predict_dat, ideas_loci_summarize));
  
  return(ret);
}

#' ideas_gwas_summarize
#' 
#' Summarize most likely cell type per GWAS analysis id.
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param .data \code{\link{ideas_predict}} output
#' @param relaxed [3] Numeric value, 1:10, rank cutoff_gt for inclusion. e.g. 3 = Must be in top 3 cell types. 
#' @param group [FALSE] TRUE = Group similar cell types. 
#' @param n_cs_ge [10] min # (>=) of cred sets to GWAS summarize
#' @param n_snps_ge [10] min # (>=) of snps across all cred sets
#' @export
#' @family ideas_predict
#' @import dplyr
ideas_gwas_summarize <- function(.data, relaxed = 3, group = FALSE,
                                 n_cs_ge = 10, n_snps_ge = 10,
                                 impala = getOption("gtx.impala", NULL)){
  # Confirm data input exist
  if(missing(.data)){
    gtx_error("ideas_gwas_summarize | no input .data specified.");
    stop();
  } else if(is_null(.data)){
    gtx_warn("ideas_gwas_summarize | ideas_predict input .data is NULL.");
    return(NULL);
  }
  
  # Confirm we have enough data to summarize
  n_credsets <- .data %>% pluck("auc") %>% nrow()
  n_cs_snps  <- .data %>% pluck("sel") %>% length()
  if(n_credsets < n_cs_ge){
    gtx_warn("ideas_gwas_summarize | Skipping ... GWAS contains {n_credsets} cred sets but min for ideas_gwas_summarize() is {n_cs_ge}.");
    return(NULL);
  }
  if(n_cs_snps < n_snps_ge){
    gtx_warn("ideas_gwas_summarize | Skipping ... GWAS contains {n_cs_snps} snps but min for ideas_gwas_summarize() is {n_snps_ge}.");
    return(NULL);
  }
  # Confirm params look correct
  if(relaxed == 0){
    gtx_error("ideas_gwas_summarize | relaxed must be numeric, 1:10.")
  }
  if(!is_logical(group)){
    gtx_error("ideas_gwas_summarize | group must be TRUE or FALSE.");
    stop();
  }
  
  cell_info <- .data$ideas_metadata
  
  a = t(.data$n)
  manual_colnames <- c("first", "second", "third", "fourth", "fifth", 
                       "sixth", "seventh", "eigth", "ninth", "tenth")
  colnames(a) <- manual_colnames
  results <- cbind(cell_info, a)

  if (is.numeric(relaxed) == TRUE & relaxed != 1) {
    results$first = rowSums(results[, manual_colnames[1:relaxed]])
    results$first = results$first/sum(results$first)
  }
  results <- results[order(results$first), ]
  
  if (group == TRUE) {
    groups = unique(dplyr::select(results, group, color))
    first = c()
    for (i in 1:dim(groups)[1]) {
      sum = dplyr::filter(results, group == groups[i, 1])$first %>% sum()
      first = c(first, sum)
    }
    z = matrix(0, nrow = dim(groups)[1], ncol = 1)
    results = cbind.data.frame(z, z, z, z, groups, first)
    colnames(results) = c("a", "b", "c", "d", "group", "Coloc", "first")
    results = results[order(results$first), ]
  }
  
  if (group == FALSE) {
    labels = results$epigenome_mnemonic
  } else {
    labels = results$group
  }
  
  results_tbl = cbind.data.frame(labels, results$first)
  colnames(results_tbl) = c("cell", "proportion")
  results_tbl <- results_tbl %>% 
    dplyr::arrange(desc(proportion)) %>% 
    dplyr::mutate(cell = as.character(cell))
  
  # Attach cell_info aka ideas_metadata to ret obj
  attr(results_tbl, "ideas_metadata") <- cell_info
  attr(results_tbl, "group") <- group
  
  gtx_debug("ideas_gwas_summarize | complete.")
  
  return(results_tbl)
}

#' ideas_gwas_plot
#' 
#' Create a plot of the \code{\link{ideas_gwas_summarize}} data.
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @author Karl Guo \email{karl.x.guo@@gsk.com}
#' @param .data \code{\link{ideas_gwas_summarize}} output
#' @param path Path for output PDF. 
#' @param cutoff_gt [0.001] Plotting "Minimum proportion of genetic signal" predicted in cell types. Range = o:1.
#' @param legend [TRUE] TRUE = display cell type legend 
#' @export
#' @family ideas_predict
#' @import dplyr
ideas_gwas_plot <- function(.data, path, cutoff_gt = 0.001, legend = TRUE){
  if(missing(.data)){
    gtx_error("ideas_gwas_plot | no input ideas_gwas_dat '.data' specified.");
    stop();
  } else if(is_null(.data)){
    gtx_warn("ideas_gwas_plot | Skipping because input ideas_gwas_dat '.data' is NULL.");
    return();
  } else {
    results = .data %>% arrange(proportion);
  }
  if (!missing(path)) {
    if(is_character(path)){
      out_file <- paste(path, ".pdf", sep = "")
      gtx_debug("ideas_plot | output file:{out_file}")
      pdf(out_file, width = 16, height = 8) 
    } else {
      gtx_error("ideas_plot | supplied path is not valid characters:{path}");
      stop();
    }
  }
  # Load ideas_metadata from attr
  ideas_metadata <- attr(results, "ideas_metadata", exact = TRUE)
  if(length(ideas_metadata) == 0){
    gtx_error("ideas_gwas_plot | input `.data` missing attribute 'ideas_metadata'")
  }
  
  group <- attr(results, "group", exact = TRUE)
  if(length(group) == 0){
    gtx_error("ideas_gwas_plot | input `.data` missing attribute 'group'")
  }
  
  if (is.numeric(cutoff_gt) == TRUE) {
    results = dplyr::filter(results, proportion > cutoff_gt)
  } else {
    gtx_error("ideas_plot | cutoff_gt parameters is not numeric. must be 0:1.")
    stop();
  }
  
  if (group == FALSE) {
    results <- 
      inner_join(results, 
                 ideas_metadata, 
                 by = c("cell" = "epigenome_mnemonic"));
    labels2 = results %>% dplyr::distinct(group, color) %>% dplyr::pull(group)
    cc2     = results %>% dplyr::distinct(group, color) %>% dplyr::pull(color)
  }
  if (group == TRUE) {
    results <- 
      inner_join(results, 
                 dplyr::distinct(ideas_metadata, group, color), 
                 by = c("cell" = "group"));
  }   
  labels  = results$cell;
  
  cclr = col2rgb(as.matrix(results$color))
  hclr = apply(cclr, 2, function(x) {
    rgb2hsv(x[1], x[2], x[3])
  })
  cc = c(apply(hclr, 2, function(x) {
    hsv(x[1], x[2], x[3])
  }))
  par(mar = c(7, 4, 3, 3))
  x = barplot(results$proportion, col = cc, cex.names = 0.6, ylab = "Proportion of genetic signal")
  if (legend == TRUE & group == FALSE) {
    legend("topleft", inset = 0.02, cex = 0.7, fill = cc2, labels2, ncol = 2)
  }
  text(x, labels = labels, srt = 45, par("usr")[3], adj = c(1, 1.1), xpd = TRUE, cex = 0.7)
  
  if (!missing(path)){
    dev.off()
  }
}