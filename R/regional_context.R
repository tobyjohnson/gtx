#' regional_context_query
#' 
#' Query pre-computed regional_context analyses
#' 
regional_context_query <- function(hgnc, hgncid){
  
}


#' regional_context_analysis
#' 
#' Calculate regional_context on-the-fly. Run PheWAS for all variants, identifying GWAS at
#' pval <= phewas_pval_le (1e-7 by default). For each variant+GWAS pair, generate credible sets
#' for each GWAS and determine if the variant is in the credible set and gather others stats on each 
#' credible set. Perform this analysis on both marginal GWAS and CLEO GWAS results.
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @param hgnc Analyze all moderate/high impact VEP variants overlapping gene position.
#' @param rs RSID(s) to analyze as a string or vector
#' @param chrom chromosome(s) to analyze as a string or vector
#' @param pos positions(s) to analyze as a single or vector
#' @param ref ref allele(s) to analyze as a string or vector. Optional w.r.t. chrom/pos.
#' @param alt alt allele(s) to analyze as a string or vector. Optional w.r.t. chrom/pos.
#' @param cpu [Default = 8] Number of cpus to use for regional context analysis (i.e. regionplot.data's)
#' @param ignore_ukb_neale [Default = TRUE] TRUE = ignore Neale UKB GWAS in PheWAS
#' @param ignore_ukb_cane [Default = TRUE] TRUE = ignore Canelaxandri UKB GWAS in PheWAS
#' @param ignore_qtls [Default = TRUE] TRUE = ignore QTLs in PheWAS (gwas_results with entity)
#' @param dbc [Default = getOption("gtx.dbConnection", NULL)] gtxconnection. 
#' @return tibble (data.frame)
regional_context_analysis <- function(hgnc, hgncid, rs, rsid, chrom, pos, ref, alt, 
                                      dbc = getOption("gtx.dbConnection", NULL),
                                      ignore_ukb_neale = TRUE, 
                                      ignore_ukb_cane  = TRUE,
                                      phewas_pval_le = 1e-7,
                                      cpu = 8, 
                                      drop_cs = TRUE){
  # --- Check gtx db connection
  gtx_debug("regional_context_analysis | validating gtx DB connection.")
  gtxdbcheck(dbc)
  
  # --- Check inputs
  if(missing(hgnc) & missing(hgncid) & missing(rs) & missing(rsid) & missing(chrom) & missing(pos)){
    gtx_fatal_stop("regional_context_analysis | missing input. Please check arguements and try again.")
  }
  
  # --- Perform PheWAS on each input
  input_tbl_nvars <- tally(input_tbl) %>% collect() %>% pull(n);
  gtx_info("regional_context_analysis | Performing PheWAS on {input_tbl_nvars} variants~GWAS pairs.")
  phewas_hits <- int_ht_phewas(hgnc   = hgnc,  hgncid = hgncid, 
                               rs     = rs,    rsid   = rsid,
                               chrom  = chrom, pos    = pos, 
                               ref    = ref,   alt    = alt,
                               ignore_ukb_neale = ignore_ukb_neale, 
                               ignore_ukb_cane  = ignore_ukb_cane, 
                               phewas_pval_le   = phewas_pval_le)
  
  # --- Perform regional_context analysis for each PheWAS hit
  gtx_info("regional_context_analysis | Performing regional context analysis on {nrow(phewas_hits)} variants.")
  region_context_data <- int_ht_regional_context(input = phewas_hits, cpu = cpu, drop_cs = drop_cs)
  
  return(region_context_data)
}

#' int_ht_phewas
#' 
#' High-throughput PheWAS
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param hgnc Will analyze all moderate/high impact VEP variants overlapping gene position (db.vep)
#' @param rs RSID(s) to analyze as a string or vector
#' @param chrom chromosome(s) to analyze as a string or vector
#' @param pos positions(s) to analyze as a single or vector
#' @param ref ref allele(s) to analyze as a string or vector. Optional w.r.t. chrom/pos.
#' @param alt alt allele(s) to analyze as a string or vector. Optional w.r.t. chrom/pos.
#' @param ignore_ukb_neale [Default = TRUE] TRUE = ignore Neale UKB GWAS in PheWAS
#' @param ignore_ukb_cane [Default = TRUE] TRUE = ignore Canelaxandri UKB GWAS in PheWAS
#' @param ignore_qtls [Default = TRUE] TRUE = ignore QTLs in PheWAS (gwas_results with entity)
#' @param dbc [Default = getOption("gtx.dbConnection", NULL)] gtxconnection. 
#' @return tibble with all PheWAS hits
int_ht_phewas <- function(hgnc, hgncid, rs, rsid, chrom, pos, ref, alt,
                          ignore_ukb_neale = TRUE, 
                          ignore_ukb_cane = TRUE,
                          ignore_qtls = TRUE, 
                          phewas_pval_le = 1e-7, 
                          dbc = getOption("gtx.dbConnection", NULL)){
  gtxdbcheck(dbc)
  # --- Verify inputs
  gtx_debug("int_ht_phewas | Checking for input params.")
  if(missing(hgnc) & missing(hgncid) & missing(rs) & missing(rsid) & missing(chrom)){
    gtx_fatal_stop("int_ht_phewas | Missing input param(s).")
  }
  
  # --- Start building SQL statements for high-throughput PheWAS
  # --- First, get SELECT sql for variants to PheWAS
  vars_sql <- int_input_2_select_snps_sql(hgnc  = hgnc,  hgncid = hgncid, 
                                          rs    = rs,    rsid   = rsid, 
                                          chrom = chrom, pos    = pos, 
                                          ref   = ref,   alt    = alt)
  
  # --- Second, create select sql for GWAS
  marg_sql <- 
    glue::glue('SELECT * FROM gwas_results
            WHERE pval <= {phewas_pval_le} AND pval IS NOT NULL')
  
  cleo_sql <- 
    glue::glue('SELECT * FROM gwas_results_cond
            WHERE pval_cond <= {phewas_pval_le} AND pval_cond IS NOT NULL')
  
  # --- Remove any GWAS we don't want
  if(isTRUE(ignore_qtls)){
    gtx_debug('int_ht_phewas | Ignoring eQTLs in PheWAS')
    marg_sql <- glue::glue('{marg_sql} AND entity IS NULL')
    # cleo_sql <- glue::glue(cleo_sql, " AND entity IS NULL") # Will need to undo this in the future?
  }
  
  if(isTRUE(ignore_ukb_neale)){
    gtx_debug('int_ht_phewas | Ignoring Neale UKB GWAS in PheWAS')
    marg_sql <- glue::glue('{marg_sql}\n AND !regexp_like(analysis, "neale")')
    cleo_sql <- glue::glue('{cleo_sql}\n AND !regexp_like(analysis, "neale")')
  }
  
  if(isTRUE(ignore_ukb_cane)){
    gtx_debug("int_ht_phewas | Ignoring canelaxandri17 UKB GWAS in PheWAS")
    marg_sql <- glue::glue('{marg_sql}\n AND !regexp_like(analysis, "canelaxandri17")')
    cleo_sql <- glue::glue('{cleo_sql}\n AND !regexp_like(analysis, "canelaxandri17")')
  }
  
  # --- PheWAS on marginal GWAS results
  gtx_debug("int_ht_phewas | starting PheWAS on marginal GWAS results . . .")
  tictoc::tic.clearlog();
  tictoc::tic();
  
  phewas_marg_sql <- glue::glue('
    WITH 
    -- Define vars_q
      vars_q AS ({vars_sql}),
    -- Define gwas_q
      gwas_q AS ({marg_sql})
    -- PheWAS by inner joining gwas_results and variants
    SELECT vars_q.input, gwas_q.chrom, gwas_q.pos, gwas_q.ref,
           gwas_q.alt, cast(NULL as integer) as signal, gwas_q.analysis, gwas_q.pval
      FROM gwas_q
      INNER JOIN /* +BROADCAST */ vars_q
      USING (chrom, pos, ref, alt)')
  
  phewas_marg <- sqlWrapper(dbc, phewas_marg_sql, uniq = FALSE, zrok = FALSE)
  
  tictoc::toc(log = TRUE, quiet = TRUE)
  gtx_debug("int_ht_phewas | PheWAS on marginal GWAS results: {tictoc::tic.log(format = TRUE)}.")
  
  # --- PheWAS on conditional GWAS results
  gtx_debug("int_ht_phewas | starting PheWAS on CLEO GWAS results . . .")
  tictoc::tic.clearlog();
  tictoc::tic();
  
  phewas_cleo_sql <- glue::glue('
    WITH 
    -- Define vars_q
      vars_q AS ({vars_sql}),
    -- Define gwas_q
      gwas_q AS ({cleo_sql})
    -- PheWAS by inner joining gwas_results_cond and variants
    SELECT vars_q.input, gwas_q.chrom, gwas_q.pos, gwas_q.ref, 
           gwas_q.alt, gwas_q.signal, gwas_q.analysis, gwas_q.pval_cond as pval
      FROM gwas_q
      INNER JOIN /* +BROADCAST */ vars_q
      USING (chrom, pos, ref, alt)')
  
  phewas_cleo <- sqlWrapper(dbc, phewas_cleo_sql, uniq = FALSE, zrok = FALSE)
  
  tictoc::toc(log = TRUE, quiet = TRUE)
  gtx_debug("int_ht_phewas | PheWAS on CLEO GWAS results: {tictoc::tic.log(format = TRUE)}.")
  
  # Union the two data streams and return
  return(as_tibble(union(phewas_marg, phewas_cleo)))
}

#' int_ht_regional_context
#' 
#' High-Throughput regional_context.
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param input data.frame with chrom, pos, ref, alt, analysis, signal; likely from phewas. 
#' @param chrom chromosome
#' @param pos position
#' @param ref reference allele
#' @param alt alt allele
#' @param analysis analysis ID
#' @param signal CLEO signal
#' @return tibble 
int_ht_regional_context <- function(input, chrom, pos, ref, alt, analysis, signal, cpu = 8, drop_cs = TRUE){
  # --- Verify inputs
  if(missing(input) & (missing(chrom) & missing(pos) & missing(ref) & missing(alt) & missing(analysis))){
    gtx_fatal_stop("int_ht_regional_context_analysis | missing input arguement(s).")
  }
  
  if(!missing(input)){
    required_cols <- c("chrom", "pos", "analysis")
    if(!all(required_cols %in% (names(input)))){
      gtx_fatal_stop('int_ht_regional_context_analysis | input is missing required cols. ',
                     'Required cols include: {paste(required_cols, collapse = ", ")}')
    }
  } else if(!(missing(analysis) & !missing(chrom) & !missing(pos) & missing(signal))){
    input = dplyr::tibble(analysis = analysis, chrom = chrom, pos = pos, ref = ref, alt = alt, signal = NA)
  } else if(!(missing(analysis) & !missing(chrom) & !missing(pos) & !missing(signal))){
    input = dplyr::tibble(analysis = analysis, chrom = chrom, pos = pos, ref = ref, alt = alt, signal = signal)
  }
  
  if(nrow(input) == 0){
    gtx_fatal_stop("int_ht_regional_context_analysis | input has no data.")
  }
  
  # --- Split marginal and conditional GWAS results for analysis
  marg_hits <- input %>% dplyr::filter(is.na(signal))
  cleo_hits <- input %>% dplyr::filter(!is.na(signal))
  
  if(nrow(marg_hits) > 0){
    marg_rc <- int_ht_regional_context_analysis(input = marg_hits, style = 'signal',  cpu = cpu, drop_cs = drop_cs)  
  } else {
    marg_rc = NULL
  }
  
  if(nrow(cleo_hits) > 0){
    cleo_rc <- int_ht_regional_context_analysis(input = cleo_hits, style = 'signals', cpu = cpu, drop_cs = drop_cs)
  } else {
    cleo_rc = NULL
  }
  
  if(!is.null(marg_rc) & !is.null(cleo_rc)){
    ret <- dplyr::union(marg_rc, cleo_rc)
  } else if(!is.null(marg_rc) & is.null(cleo_rc)){
    ret <- marg_rc
  } else if(is.null(marg_rc) & !is.null(cleo_rc)){
    ret <- cleo_rc
  } else {
    gtx_fatal_stop("int_ht_regional_context_analysis | Unable to determine how to properly merge results.")
  }
  
  # Re-establish DB connection
  gtxconnect(use_database = config_db())
  
  return(as_tibble(ret))
}

#' int_ht_regional_context_analysis
#' 
#' High-throughput regional_context_analysis
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param input int_ht_phewas object 
#' @param style CLEO style {'signal'|'signals'}
#' @param cpu [Default=8]
#' @return tibble 
int_ht_regional_context_analysis <- function(input, style, cpu = 8, drop_cs = TRUE, ...){
  if(missing(input) | nrow(input) == 0){
    gtx_fatal_stop("int_ht_regional_context_analysis | missing input arguement(s) or data.")
  }
  
  if(missing(style) | (style != 'signal' & style != 'signals')){
    gtx_fatal_stop("int_ht_regional_context_analysis | Must use style = 'signal' OR 'signals'.")
  }
  
  # --- Remove database connections
  futile.logger::flog.threshold(ERROR);
  options("gtx.dbConnection" = NULL);
  gtxcache(disconnect = TRUE);
  futile.logger::flog.threshold(INFO);
  if(nrow(input) > 1 & cpu > 1){
    future::plan(multisession, workers = as.integer(cpu))
  } else {
    future::plan(sequential)
  }
  
  if(style == 'signal'){
    ret <- 
      input %>% 
      dplyr::mutate(cs = furrr::future_pmap(list(analysis, chrom, pos), 
                                     ~int_ht_cred_set_wrapper(analysis = ..1, 
                                                              chrom    = ..2, 
                                                              pos      = ..3,
                                                              style    = 'signal')))
    
  } else if(style == 'signals'){
    ret <- 
      input %>% 
      dplyr::mutate(cs = furrr::future_pmap(list(analysis, chrom, pos, signal), 
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
    dplyr::mutate(var_in_cs  = purrr::pmap_lgl(list(cs, chrom, pos, ref, alt), 
                                 ~dplyr::filter(..1, chrom == ..2 & pos == ..3 & 
                                           ref == ..4 & alt == ..5 & cs_signal == TRUE) %>% 
                                   nrow >= 1)) %>% 
    dplyr::mutate(var_pp     = purrr::pmap_dbl(list(cs, chrom, pos, ref, alt), 
                                 ~dplyr::filter(..1, chrom == ..2 & pos == ..3 &
                                           ref == ..4 & alt == ..5 ) %>% 
                                   dplyr::pull(pp_signal))) %>% 
    dplyr::mutate(n_in_cs    = purrr::map_dbl(cs, nrow)) %>% 
    dplyr::mutate(cs_th_pp   = purrr::map_dbl(cs, ~max(.$pp_signal))) %>% 
    dplyr::mutate(cs_th_pval = purrr::map_dbl(cs, ~min(.$pval, na.rm = TRUE))) %>% 
    dplyr::mutate(cs_th_cord = purrr::map(cs, ~top_n(., 1, pp_signal) %>% dplyr::select(pos, ref, alt))) %>% 
    tidyr::unnest(cs_th_cord) %>% 
    dplyr::rename(th_pos = pos1, th_ref = ref1, th_alt = alt1) 
  
  # Drop the cred sets. TRUE by default
  if(isTRUE(drop_cs)){
    ret = ret %>% dplyr::select(-cs)
  }

  return(ret)
}


#' int_ht_cred_set_wrapper
#' 
#' int_ht_cred_set_wrapper
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param chrom chromosome
#' @param pos position
#' @param analysis analysis ID
#' @param ... Additional arguements for regionplot.data
#' @return tibble 
int_ht_cred_set_wrapper <- function(analysis, chrom, pos, style, ...){
  if(missing(chrom) | missing(pos) | missing(analysis) | missing(style)){
    gtx_fatal_stop("int_ht_cred_set_wrapper | missing input arguement(s).")
  }
  
  futile.logger::flog.threshold(ERROR) # Turn off gtxconnect & regionplot INFO msgs
  
  if(is.null(getOption("gtx.dbConnection", NULL))){
    gtx_debug("int_ht_cred_set_wrapper | Establishing gtxconnection.")
    gtxconnect(use_database = config_db(), cache = FALSE)
  } else {
    gtx_debug("int_ht_cred_set_wrapper | Using pre-established gtxconnection.")
  }
  gtxdbcheck(getOption("gtx.dbConnection", NULL))
  
  if(style == 'signal'){
    ret <- 
      regionplot.data(analysis = analysis, chrom = chrom, pos = pos, style = style, ...)
  } else if(style == 'signals'){
    ret <- 
      regionplot.data(analysis = analysis, chrom = chrom, pos = pos, style = style, ...) %>% # Pass signal in ...
      dplyr::select(-dplyr::matches("signal")) %>% 
      dplyr::rename(pp_signal = pp_cleo, cs_signal = cs_cleo)
  }
  
  ret <- ret %>% dplyr::filter(cs_signal == TRUE | (chrom == !! chrom & pos == !! pos))
  
  return(ret)
}


#' int_input_2_select_snps_sql
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
#' @return string to select variants
int_input_2_select_snps_sql <- function(hgnc, hgncid, rs, rsid, chrom, pos, ref, alt){
  if(!missing(hgnc) | !missing(hgncid)){
    gtx_debug("int_input_tbl | hgnc symbol input - finding all VEP variants for gene(s).")
    
    if(!missing(hgncid)){    input <- dplyr::tibble(input = hgncid) }
    else if(!missing(hgnc)){ input <- dplyr::tibble(input = hgnc)   }
    
    input <- input %>% dplyr::mutate(where_clause = purrr::map_chr(input, ~gtxwhere(hgncid = .)))
    
    select_statement <- glue::glue('
      WITH 
        t1 as (SELECT * FROM genes 
               WHERE {glue::glue_collapse(input$where_clause, sep = " OR ")})
      SELECT hgncid AS input, vep.chrom, vep.pos, vep.ref, vep.alt 
        FROM vep
        INNER JOIN /* +BROADCAST */ t1
        USING (chrom)
        WHERE pos >= pos_start AND pos <= pos_end')
    
  } else if(!missing(rs) | !missing(rsid)){
    gtx_debug("int_input_tbl | rs ID input - finding all genomic coordinates.")
    if(!missing(rsid))   { input <- dplyr::tibble(input = rsid) }
    else if(!missing(rs)){ input <- dplyr::tibble(input = rs)   }
    
    input <- input %>% dplyr::mutate(where_clause = purrr::map_chr(input, ~gtxwhere(rs = .)))
    
    select_statement <- glue::glue('
      SELECT CONCAT("rs", cast(rsid as string)) AS input, chrom, pos, ref, alt
        FROM sites
        WHERE {glue::glue_collapse(input$where_clause, sep = " OR ")}')
    
  } else if(!missing(ref)){
    gtx_debug("int_input_tbl | chrom,pos,ref,alt input - confirming genomic coordinates.")
    input <- 
      dplyr::tibble(chrom = chrom, pos = pos, ref = ref, alt = alt) %>% 
      dplyr::mutate(where_clause = purrr::pmap_chr(list(chrom, pos, ref, alt),
                                     ~gtxwhere(chrom = ..1, pos = ..2, 
                                               ref   = ..3, alt = ..4)))
    
    select_statement <- glue::glue('
      SELECT CONCAT("chr", chrom, "_", cast(pos as string), "_", ref, "_", alt) AS input, chrom, pos, ref, alt
        FROM sites
        WHERE {glue::glue_collapse(input$where_clause, sep = " OR ")}')
    
  } else if(!missing(chrom)) { 
    gtx_debug("int_input_tbl | chrom,pos,ref - ID all genomic coordinates.")
    input <- 
      dplyr::tibble(chrom = chrom, pos = pos) %>% 
      dplyr::mutate(where_clause = purrr::map2_chr(chrom, pos, ~gtxwhere(chrom = .x, pos = .y)))
    
    select_statement <- glue::glue('
      SELECT CONCAT("chr", chrom, "_", cast(pos as string)) AS input, chrom, pos, ref, alt
        FROM sites
        WHERE {glue::glue_collapse(input$where_clause, sep = " OR ")}')
    
  } else {
    gtx_fatal_stop("int_input_tbl | unable to determine which input(s) to use.")
  }
  
  return(select_statement)
}