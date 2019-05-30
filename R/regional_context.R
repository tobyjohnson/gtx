#' regional_context_query
#' 
#' Query pre-computed regional_context analyses
#' 
regional_context_query <- function(hgnc, hgncid){
  
}


#' regional_context_analysis
#' 
#' Calculate regional_context on-the-fly. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @param hgnc Will analyze all moderate/high impact VEP variants overlapping gene position (db.vep)
#' @param rs RSID(s) to analyze as a string or vector
#' @param chrom chromosome(s) to analyze as a string or vector
#' @param pos positions(s) to analyze as a single or vector
#' @param ref ref allele(s) to analyze as a string or vector. Optional w.r.t. chrom/pos.
#' @param alt alt allele(s) to analyze as a string or vector. Optional w.r.t. chrom/pos.
#' @param impala [getOption("gtx.impala", NULL)] Implyr impala connection
regional_context_analysis <- function(hgnc, hgncid, rs, rsid, chrom, pos, ref, alt, 
                                      dbc = getOption("gtx.dbConnection", NULL),
                                      impala = getOption("gtx.impala", NULL),
                                      ignore_ukb_neale = TRUE, ignore_ukb_cane = TRUE,
                                      ignore_qtls = TRUE, phewas_pval_le = 1e-7,
                                      cpu = 8){
  # --- Check inputs
  gtx_debug("regional_context_analysis | validating GTX odbc connection.")
  gtxdbcheck(dbc)
  
  gtx_debug("regional_context_analysis | validating impala connection.")
  impala_dbc <- validate_impala(impala = impala)
  
  if(missing(hgnc) & missing(hgncid) & missing(rs) & missing(rsid) & missing(chrom) & missing(pos)){
    gtx_fatal_stop("regional_context_analysis | missing input. Please check arguements and try again.")
  }
  
  # --- Harmonize inputs to a table reference in RDIP
  input_tbl <- int_input_tbl(hgnc   = hgnc,  hgncid = hgncid, 
                             rs     = rs,    rsid   = rsid,
                             chrom  = chrom, pos    = pos, 
                             ref    = ref,   alt    = alt, 
                             impala = impala_dbc)
  
  # --- Perform PheWAS on each input
  phewas_hits <- int_ht_phewas(input = input_tbl, impala = impala_dbc,
                               ignore_ukb_neale = ignore_ukb_neale, 
                               ignore_ukb_cane = ignore_ukb_cane, 
                               ignore_qtls = ignore_qtls, 
                               phewas_pval_le = phewas_pval_le)
  
  # --- Perform regional_context analysis for each sig PheWAS hit
  
  
}

#' int_input_tbl
#' 
#' Get an impala table reference for the various inputs harmonized to chrom/pos/ref/alt
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param hgnc Will analyze all moderate/high impact VEP variants overlapping gene position (db.vep)
#' @param rs RSID(s) to analyze as a string or vector
#' @param chrom chromosome(s) to analyze as a string or vector
#' @param pos positions(s) to analyze as a single or vector
#' @param ref ref allele(s) to analyze as a string or vector. Optional w.r.t. chrom/pos.
#' @param alt alt allele(s) to analyze as a string or vector. Optional w.r.t. chrom/pos.
#' @param impala [getOption("gtx.impala", NULL)] Implyr impala connection
#' @return implyr_tbl reference to a list of the input variants.
int_input_tbl <- function(hgnc, hgncid, rs, rsid, chrom, pos, ref, alt,
                          impala = getOption("gtx.impala", NULL)){
  gtx_debug("int_input_tbl | validating impala connection.")
  impala_dbc <- validate_impala(impala = impala)
  
  if(!missing(hgnc) | !missing(hgncid)){
    vep_tbl   <- tbl(impala_dbc, glue("{config_db()}.vep"))
    
    if(!missing(hgncid)){    input <- tibble(input = hgncid) }
    else if(!missing(hgnc)){ input <- tibble(input = hgnc)   }
    
    input <- input %>% mutate(where_clause = map_chr(input, ~gtxwhere(hgncid = .)))
    
    select_statement <-
      glue('SELECT * FROM {config_db()}.genes WHERE ',
           glue_collapse(input$where_clause, sep = " OR "))
    
    
    input_tbl <- 
      tbl(impala_dbc, sql(select_statement)) %>% 
      select(-genetype, -ensemblversion) %>% 
      inner_join(vep_tbl, by = c("chrom")) %>% 
      filter(pos >= pos_start & pos <= pos_end) %>% 
      select(input = hgncid, chrom, pos, ref, alt)
    
  } else if(!missing(rs) | !missing(rsid)){
    if(!missing(rsid))   { input <- tibble(input = rsid) }
    else if(!missing(rs)){ input <- tibble(input = rs)   }
    
    input <- input %>% mutate(where_clause = map_chr(input, ~gtxwhere(rs = .)))
    
    select_statement <-
      glue('SELECT * FROM {config_db()}.sites WHERE ',
           glue_collapse(input$where_clause, sep = " OR "))
    
    input_tbl <- 
      tbl(impala_dbc, sql(select_statement)) %>% 
      mutate(input = as.character(rsid)) %>% 
      mutate(input = concat("rs", input)) %>% 
      select(input, chrom, pos, ref, alt)
    
  } else if(!missing(ref)){
    input <- tibble(chrom = chrom, pos = pos, ref = ref, alt = alt)
    
    input <- 
      input %>% 
      mutate(where_clause = pmap_chr(list(chrom, pos, ref, alt),
                                     ~gtxwhere(chrom = ..1, pos = ..2, 
                                               ref   = ..3, alt = ..4)))
    
    select_statement <-
      glue('SELECT * FROM {config_db()}.sites WHERE ',
           glue_collapse(input$where_clause, sep = " OR "))
    
    input_tbl <- 
      tbl(impala_dbc, sql(select_statement)) %>% 
      mutate(input = as.character(pos)) %>% 
      mutate(input = concat("chr", chrom, "_", input, "_", ref, "_", alt)) %>% 
      select(input, chrom, pos, ref, alt)
    
  } else if(!missing(chrom)) { 
    input <- tibble(chrom = chrom, pos = pos)
    
    input <- 
      input %>% 
      mutate(where_clause = map2_chr(chrom, pos, ~gtxwhere(chrom = .x, pos = .y)))
    
    select_statement <-
      glue('SELECT * FROM {config_db()}.sites WHERE ',
           glue_collapse(input$where_clause, sep = " OR "))
    
    input_tbl <- 
      tbl(impala_dbc, sql(select_statement)) %>% 
      mutate(input = as.character(pos)) %>% 
      mutate(input = concat("chr", chrom, "_", input)) %>% 
      select(input, chrom, pos, ref, alt)
  } else {
    gtx_fatal_stop("int_input_tbl | unable to determine which input(s) to use.")
  }
  
  attr(input_tbl, "int_input_tbl") <- TRUE
  return(input_tbl)
}

#' int_ht_phewas
#' 
#' High-throughput PheWAS of an int_input_tbl.
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param input int_input_tbl object for PheWAS
#' @param ignore_ukb_neale [Default = TRUE] TRUE = ignore Neale UKB GWAS in PheWAS
#' @param ignore_ukb_cane [Default = TRUE] TRUE = ignore Canelaxandri UKB GWAS in PheWAS
#' @param ignore_qtls [Default = TRUE] TRUE = ignore QTLs in PheWAS
#' @param impala [getOption("gtx.impala", NULL)] Implyr impala connection
#' @return tibble with all PheWAS hits
int_ht_phewas <- function(input, ignore_ukb_neale = TRUE, ignore_ukb_cane = TRUE,
                          ignore_qtls = TRUE, phewas_pval_le = 1e-7, 
                          impala = getOption("gtx.impala", NULL)){
  # --- Verify inputs
  gtx_debug("int_ht_phewas | validating input as int_input_tbl.")
  if(attr(input, "int_input_tbl", exact = TRUE) != TRUE){
    gtx_fatal_stop("int_ht_phewas | input is not an int_input_tbl.")
  }
  
  gtx_debug("int_ht_phewas | validating impala connection.")
  impala_dbc <- validate_impala(impala = impala)
  
  # --- Make references to both GWAS tables for marginal and CLEO results
  marg_tbl  <- 
    tbl(impala_dbc, glue("{config_db()}.gwas_results")) %>% 
    filter(pval <= phewas_pval_le) 
  
  cleo_tbl <- 
    tbl(impala_dbc, glue("{config_db()}.gwas_results_cond")) %>% 
    filter(pval_cond <= phewas_pval_le) 
  
  # --- Remove any GWAS we don't want
  if(isTRUE(ignore_qtls)){
    marg_tbl <- marg_tbl %>% filter(is.na(entity))
    # cleo_tbl <- cleo_tbl %>% filter(is.na(entity)) # Will need to undo this in the future?
  }
  
  if(isTRUE(ignore_ukb_neale)){
    marg_tbl <- marg_tbl %>% filter(!regexp_like(analysis, "neale"))
    cleo_tbl <- cleo_tbl %>% filter(!regexp_like(analysis, "neale"))
  }
  
  if(isTRUE(ignore_ukb_cane)){
    marg_tbl <- marg_tbl %>% filter(!regexp_like(analysis, "canelaxandri17"))
    cleo_tbl <- cleo_tbl %>% filter(!regexp_like(analysis, "canelaxandri17"))
  }
  
  # --- PheWAS on marginal GWAS results
  gtx_debug("int_ht_phewas | starting PheWAS on marginal GWAS results . . .")
  tictoc::tic.clearlog();
  tictoc::tic();
  phewas_marg <- 
    inner_join(marg_tbl, 
               input, 
               by = c("chrom", "pos", "ref", "alt")) %>% 
    collect() %>% 
    mutate(signal = NA_integer_) %>% 
    select(input, chrom, pos, ref, alt, signal, analysis, pval)
  tictoc::toc(log = TRUE, quiet = TRUE)
  gtx_debug("int_ht_phewas | PheWAS on marginal GWAS results: {tictoc::tic.log(format = TRUE)}.")
  
  # --- PheWAS on conditional GWAS results
  gtx_debug("int_ht_phewas | starting PheWAS on CLEO GWAS results . . .")
  tictoc::tic.clearlog();
  tictoc::tic();
  phewas_cleo <- 
    inner_join(cleo_tbl, 
               input, 
               by = c("chrom", "pos", "ref", "alt")) %>% 
    collect() %>% 
    select(input, chrom, pos, ref, alt, signal, analysis, pval = pval_cond)
  tictoc::toc(log = TRUE, quiet = TRUE)
  gtx_debug("int_ht_phewas | PheWAS on CLEO GWAS results: {tictoc::tic.log(format = TRUE)}.")
  
  # Union the two data streams and return
  return(union(phewas_marg, phewas_cleo))
}

#' int_ht_regional_context
#' 
#' INTernal High-Throughput regional_context
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param input int_ht_phewas object 
#' @param chrom TODO
#' @param pos TODO
#' @param ref TODO
#' @param alt TODO
#' @param analysis TODO
#' @return tibble 
int_ht_regional_context <- function(input, chrom, pos, ref, alt, analysis, cpu = 8){
  # --- Verify inputs
  if(missing(input) & (missing(chrom) & missing(pos) & missing(ref) & missing(alt) & missing(analysis))){
    gtx_fatal_stop("int_ht_regional_context_analysis | missing input arguement(s).")
  }
  
  if(!missing(input)){
    required_cols <- c("chrom", "pos", "ref", "alt", "analysis")
    if(!all(required_cols %in% (names(input)))){
      gtx_fatal_stop('int_ht_regional_context_analysis | input is missing required cols. ',
                     'Required cols include: {paste(required_cols, collapse = ", ")}')
    }
  } else if(!(missing(analysis) & !missing(chrom) & !missing(pos))){
    input = tibble(analysis = analysis, chrom = chrom, pos = pos, ref = ref, alt = alt)
  }
  
  if(nrow(input) == 0){
    gtx_fatal_stop("int_ht_regional_context_analysis | input has no data.")
  }
  
  # --- Split marginal and conditional GWAS results for analysis
  marg_hits <- input %>% filter(is.na(signal))
  cleo_hits <- input %>% filter(!is.na(signal))
  
  marg_rc <- int_ht_regional_context_analysis(input = marg_hits, style = 'signal',  cpu = cpu)
  cleo_rc <- int_ht_regional_context_analysis(input = cleo_hits, style = 'signals', cpu = cpu)
  
  ret <- union(marg_rc, cleo_rc)
  
  return(ret)
}

#' int_ht_regional_context_analysis
#' 
#' High-throughput regional_context_analysis
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param input int_ht_phewas object 
#' @param style TODO
#' @param cpu [Default=8]
#' @return tibble 
int_ht_regional_context_analysis <- function(input, style, cpu = 8, ...){
  if(missing(input) | nrow(input) == 0){
    gtx_fatal_stop("int_ht_regional_context_analysis | missing input arguement(s) or data.")
  }
  
  if(missing(style) | (style != 'signal' & style != 'signals')){
    gtx_fatal_stop("int_ht_regional_context_analysis | Must use style = 'signal' OR 'signals'.")
  }
  
  # --- Remove database connections if making multi-threaded
  if(nrow(input) > 1 & cpu > 1){
    options(gtx.dbConnection = NULL);  
    gtxcache(disconnect = TRUE);
    plan(multiprocess, workers = as.integer(cpu))
  } else {
    plan(sequential)
  }
  futile.logger::flog.threshold(ERROR)
  
  if(style == 'signal'){
    ret <- 
      input %>% 
      mutate(cs = future_pmap(list(analysis, chrom, pos), 
                              ~int_ht_cred_set_wrapper(analysis = ..1, 
                                                       chrom    = ..2, 
                                                       pos      = ..3,
                                                       style    = 'signal')))
    
  } else if(style == 'signals'){
    ret <- 
      input %>% 
      mutate(cs = future_pmap(list(analysis, chrom, pos, signal), 
                              ~int_ht_cred_set_wrapper(analysis = ..1, 
                                                       chrom    = ..2, 
                                                       pos      = ..3,
                                                       style    = 'signals',
                                                       signal   = ..4))) 
  } else {
    gtx_fatal_stop("int_ht_regional_context_analysis | Unable to determine style type.")
  }
  
  ret <- 
    ret %>% 
    mutate(var_in_cs  = pmap_lgl(list(cs, chrom, pos, ref, alt), 
                                 ~filter(..1, chrom == ..2 & pos == ..3 & 
                                           ref == ..4 & alt == ..5 & cs_signal == TRUE) %>% 
                                   nrow >= 1)) %>% 
    mutate(var_pp     = pmap_dbl(list(cs, chrom, pos, ref, alt), 
                                 ~filter(..1, chrom == ..2 & pos == ..3 &
                                           ref == ..4 & alt == ..5 ) %>% 
                                   pull(pp_signal))) %>% 
    mutate(n_in_cs     = map_dbl(cs, nrow)) %>% 
    mutate(cs_th_pp   = map_dbl(cs, ~max(.$pp_signal))) %>% 
    mutate(cs_th_cord = map(cs, ~top_n(., 1, pp_signal) %>% select(pos, ref, alt))) %>% 
    unnest(cs_th_cord) %>% 
    rename(th_pos = pos1, th_ref = ref1, th_alt = alt1) %>% 
    select(-cs)
  
  futile.logger::flog.threshold(INFO)
  return(ret)
}


#' int_ht_cred_set_wrapper
#' 
#' int_ht_cred_set_wrapper
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param chrom TODO
#' @param pos TODO
#' @param analysis TODO
#' @param ... TODO
#' @return tibble 
int_ht_cred_set_wrapper <- function(analysis, chrom, pos, style, ...){
  if(missing(chrom) | missing(pos) | missing(analysis) | missing(style)){
    gtx_fatal_stop("int_ht_cred_set_wrapper | missing input arguement(s).")
  }
  gtxconnect(use_database = config_db(), cache = FALSE)
  
  if(style == 'signal'){
    ret <- 
      regionplot.data(analysis = analysis, chrom = chrom, pos = pos, style = style, ...)
  } else if(style == 'signals'){
    ret <- 
      regionplot.data(analysis = analysis, chrom = chrom, pos = pos, style = style, ...) %>% 
      select(-matches("signal")) %>% 
      rename(pp_signal = pp_cleo, cs_signal = cs_cleo)
  }
  
  ret <- ret %>% filter(cs_signal == TRUE | (chrom == !! chrom & pos == !! pos))
  
  return(ret)
}
