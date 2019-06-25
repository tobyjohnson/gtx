#' regional_context_query
#' 
#' Query pre-computed regional_context analyses
#' 
regional_context_query <- function(hgnc, hgncid){
  
}


#' regional_context_analysis
#' 
#' Calculate regional_context on-the-fly. 
#' *Critical - DO NOT run this function in a loop, use vector inputs*.
#' Run PheWAS for all variants, identifying GWAS at the significance threshold
#' (1e-7 by default). If a gene is used as input, the function will find all moderate/high impact
#' variants from the VEP table and use all the VEP variants for PheWAS. 
#' In the Phewas, the default behavior is to also ignore Ben Neale v1 & Canelaxandri17 UKB,
#' and all xQTL data. For each variant+GWAS pair, generate credible sets and determine if the variant is in the credible set and gather 
#' others stats on each credible set. Perform this analysis on both marginal GWAS and CLEO GWAS results. Due to the 
#' high-throughput nature of this function, it is suggested to reduce logging info when running this function using: 
#' futile.logger::flog.threshold(ERROR). This function will also use multiple cores if available.
#' 
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
#' @param ... regionplot.data params
#' @return tibble (data.frame)
regional_context_analysis <- function(hgnc, hgncid, rs, rsid, chrom, pos, ref, alt, 
                                      ignore_ukb_neale = TRUE, 
                                      ignore_ukb_cane  = TRUE,
                                      ignore_qtls      = TRUE, 
                                      phewas_pval_le = 1e-7,
                                      cpu = 8, drop_cs = TRUE, 
                                      dbc = getOption("gtx.dbConnection", NULL), ...){
  # --- Check gtx db connection
  gtx_debug("regional_context_analysis | validating gtx DB connection.")
  gtxdbcheck(dbc)
  
  # --- Check inputs
  if(missing(hgnc) & missing(hgncid) & missing(rs) & missing(rsid) & missing(chrom) & missing(pos)){
    gtx_fatal_stop("regional_context_analysis | missing input. Please check arguements and try again.")
  }
  
  # --- Perform PheWAS on each input
  phewas_hits <- 
    int_ht_phewas_wrapper(hgnc   = hgnc,  hgncid = hgncid, 
                          rs     = rs,    rsid   = rsid,
                          chrom  = chrom, pos    = pos, 
                          ref    = ref,   alt    = alt,
                          ignore_ukb_neale = ignore_ukb_neale, 
                          ignore_ukb_cane  = ignore_ukb_cane,
                          ignore_qtls      = ignore_qtls, 
                          phewas_pval_le   = phewas_pval_le, 
                          dbc = dbc)
  
  # --- Perform regional_context analysis for each PheWAS hit
  gtx_info("Performing regional context analysis on {nrow(phewas_hits)} variants+GWAS pairs.")
  region_context_data <- int_ht_regional_context(input = phewas_hits, cpu = cpu, dbc = dbc, drop_cs = drop_cs, ...)
  
  return(region_context_data)
}

#' int_ht_phewas_wrapper
#' 
#' High-throughput PheWAS - high level wrapper pointing to int_ht_phewas. However, due to 
#' limitations of impala, need to chunk the large phewas joins into smaller chunks. 
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
int_ht_phewas_wrapper <- function(hgnc, hgncid, rs, rsid, chrom, pos, ref, alt,
                                  ignore_ukb_neale = TRUE, 
                                  ignore_ukb_cane = TRUE,
                                  ignore_qtls = TRUE, 
                                  phewas_pval_le = 1e-7, 
                                  chunk_size = 500,
                                  dbc = getOption("gtx.dbConnection", NULL)){
  # --- Verify inputs
  gtx_debug("int_ht_phewas | Checking for input params.")
  if(missing(hgnc) & missing(hgncid) & missing(rs) & missing(rsid) & missing(chrom)){
    gtx_fatal_stop("int_ht_phewas | Missing input param(s).")
  }
  
  if(!missing(hgnc) | !missing(hgncid)){
    if(!missing(hgncid)){    input <- dplyr::tibble(hgncid = hgncid) }
    else if(!missing(hgnc)){ input <- dplyr::tibble(hgncid = hgnc)   }
    
    if(nrow(input) > chunk_size){
      gtx_debug("Spliting input into chunks for PheWAS.")
      ret <- 
        input %>% 
        mutate(tile = ntile(row_number(), ceiling(nrow(.)/chunk_size))) %>% 
        group_by(tile) %>% 
        nest(.key = "data") %>% 
        mutate(phewas_data = map(data, ~int_ht_phewas(hgncid = .$hgncid, 
                                                      ignore_ukb_neale = ignore_ukb_neale, 
                                                      ignore_ukb_cane  = ignore_ukb_cane,
                                                      ignore_qtls      = ignore_qtls, 
                                                      phewas_pval_le   = phewas_pval_le))) %>% 
        select(-data, -tile) %>% 
        unnest(phewas_data)
      
    } else {
      ret <- int_ht_phewas(hgncid = input$hgncid,
                           ignore_ukb_neale = ignore_ukb_neale, 
                           ignore_ukb_cane  = ignore_ukb_cane,
                           ignore_qtls      = ignore_qtls, 
                           phewas_pval_le   = phewas_pval_le)
    }
  } else if(!missing(rs) | !missing(rsid)){
    if(!missing(rsid))   { input <- dplyr::tibble(rsid = rsid) }
    else if(!missing(rs)){ input <- dplyr::tibble(rsid = rs)   }
    
    if(nrow(input) > chunk_size){
      gtx_debug("Spliting input into chunks for PheWAS.")
      ret <- 
        input %>% 
        mutate(tile = ntile(row_number(), ceiling(nrow(.)/chunk_size))) %>% 
        group_by(tile) %>% 
        nest(.key = "data") %>% 
        mutate(phewas_data = map(data, ~int_ht_phewas(rsid = .$rsid, 
                                                      ignore_ukb_neale = ignore_ukb_neale, 
                                                      ignore_ukb_cane  = ignore_ukb_cane,
                                                      ignore_qtls      = ignore_qtls, 
                                                      phewas_pval_le   = phewas_pval_le))) %>% 
        select(-data, -tile) %>% 
        unnest(phewas_data)
      
    } else {
      ret <- int_ht_phewas(rsid = input$rsid, 
                           ignore_ukb_neale = ignore_ukb_neale, 
                           ignore_ukb_cane  = ignore_ukb_cane,
                           ignore_qtls      = ignore_qtls, 
                           phewas_pval_le   = phewas_pval_le)
    } 
  } else if(!missing(ref)){
    input <- dplyr::tibble(chrom = chrom, pos = pos, ref = ref, alt = alt) 
    if(nrow(input) > chunk_size){
      gtx_debug("Spliting input into chunks for PheWAS.")
      ret <- 
        input %>% 
        mutate(tile = ntile(row_number(), ceiling(nrow(.)/chunk_size))) %>% 
        group_by(tile) %>% 
        nest(.key = "data") %>% 
        mutate(phewas_data = map(data, ~int_ht_phewas(chrom = .$chrom, pos = .$pos, 
                                                      ref = .$ref, alt = .$alt,
                                                      ignore_ukb_neale = ignore_ukb_neale, 
                                                      ignore_ukb_cane  = ignore_ukb_cane,
                                                      ignore_qtls      = ignore_qtls, 
                                                      phewas_pval_le   = phewas_pval_le))) %>% 
        select(-data, -tile) %>% 
        unnest(phewas_data)
      
    } else {
      ret <- int_ht_phewas(chrom = input$chrom, pos = input$pos, 
                           ref = input$ref, alt = input$alt, 
                           ignore_ukb_neale = ignore_ukb_neale, 
                           ignore_ukb_cane  = ignore_ukb_cane,
                           ignore_qtls      = ignore_qtls, 
                           phewas_pval_le   = phewas_pval_le)
    }
  } else if(!missing(chrom)) { 
    gtx_debug("int_input_tbl | chrom,pos,ref - ID all genomic coordinates.")
    input <- dplyr::tibble(chrom = chrom, pos = pos) 
    if(nrow(input) > chunk_size){
      gtx_debug("Spliting input into chunks for PheWAS.")
      ret <- 
        input %>% 
        mutate(tile = ntile(row_number(), ceiling(nrow(.)/chunk_size))) %>% 
        group_by(tile) %>% 
        nest(.key = "data") %>% 
        mutate(phewas_data = map(data, ~int_ht_phewas(chrom = .$chrom, pos = .$pos,
                                                      ignore_ukb_neale = ignore_ukb_neale, 
                                                      ignore_ukb_cane  = ignore_ukb_cane,
                                                      ignore_qtls      = ignore_qtls, 
                                                      phewas_pval_le   = phewas_pval_le))) %>% 
        select(-data, -tile) %>% 
        unnest(phewas_data)
      
    } else {
      ret <- int_ht_phewas(chrom = input$chrom, 
                           pos   = input$pos, 
                           ignore_ukb_neale = ignore_ukb_neale, 
                           ignore_ukb_cane  = ignore_ukb_cane,
                           ignore_qtls      = ignore_qtls, 
                           phewas_pval_le   = phewas_pval_le)
    }
  }
  return(ret)
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
#' @return tibble with all PheWAS hits!
int_ht_phewas <- function(hgnc, hgncid, rs, rsid, chrom, pos, ref, alt,
                          ignore_ukb_neale = TRUE, 
                          ignore_ukb_cane = TRUE,
                          ignore_qtls = TRUE, 
                          phewas_pval_le = 1e-7, 
                          dbc = getOption("gtx.dbConnection", NULL)){
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
  
  phewas_vars <- sqlWrapper(dbc, vars_sql, uniq = FALSE, zrok = FALSE)
  gtx_info("{nrow(phewas_vars)} variant(s) identified for PheWAS.")
  # Impala isn't smart enough to limit the partitions searched for only the actual chroms used in
  # the "join on (chrom)". Manually build WHERE clause to explicitly limit partitions searched.
  phewas_chroms <- phewas_vars %>% distinct(chrom)
  phewas_chroms_where_clause <- 
    glue::glue('(',
               glue::glue_collapse(glue::glue("(chrom=\"{phewas_chroms$chrom}\")"), sep = " OR "),
               ')')
  
  # --- Second, create select sql for GWAS
  marg_sql <- glue::glue('
    SELECT * FROM gwas_results
      LEFT SEMI JOIN analyses
      USING (analysis)
      WHERE pval <= {phewas_pval_le} 
      AND pval IS NOT NULL
      AND {phewas_chroms_where_clause}')
  
  cleo_sql <- glue::glue('
    SELECT * FROM gwas_results_cond
      LEFT SEMI JOIN analyses
      USING (analysis)
      WHERE pval_cond <= {phewas_pval_le} 
      AND pval_cond IS NOT NULL
      AND {phewas_chroms_where_clause}')
  
  # --- Remove any GWAS we don't want
  if(isTRUE(ignore_qtls)){
    gtx_debug('int_ht_phewas | Ignoring eQTLs in PheWAS')
    marg_sql <- glue::glue('{marg_sql} AND entity IS NULL')
    # cleo_sql <- glue::glue(cleo_sql, " AND entity IS NULL") # Will need to undo this in the future?
  }
  
  if(isTRUE(ignore_ukb_neale)){
    gtx_debug('int_ht_phewas | Ignoring Neale UKB GWAS in PheWAS')
    marg_sql <- glue::glue('{marg_sql} AND !regexp_like(analysis, "neale")')
    cleo_sql <- glue::glue('{cleo_sql} AND !regexp_like(analysis, "neale")')
  }
  
  if(isTRUE(ignore_ukb_cane)){
    gtx_debug("int_ht_phewas | Ignoring canelaxandri17 UKB GWAS in PheWAS")
    marg_sql <- glue::glue('{marg_sql} AND !regexp_like(analysis, "canelaxandri17")')
    cleo_sql <- glue::glue('{cleo_sql} AND !regexp_like(analysis, "canelaxandri17")')
  }
  
  # --- PheWAS on marginal GWAS results
  gtx_info("Starting PheWAS of marginal GWAS results.")
  tictoc::tic.clearlog();
  tictoc::tic();
  
  phewas_marg_sql <- glue::glue('
      WITH 
      -- Define vars_q
      vars_q AS ({vars_sql}),
      -- Define gwas_q
      gwas_q AS ({marg_sql}),
      -- PheWAS by inner joining gwas_results and variants subqueries
      phewas_q AS 
        (SELECT vars_q.input, gwas_q.chrom, gwas_q.pos, gwas_q.ref, gwas_q.alt,
                cast(NULL as integer) as signal, gwas_q.analysis, gwas_q.pval, gwas_q.beta
         FROM gwas_q
         INNER JOIN vars_q
         USING (chrom, pos, ref, alt))
      -- Append GWAS top hits.
      SELECT phewas_q.input, phewas_q.chrom, phewas_q.pos, phewas_q.ref, phewas_q.alt,
             phewas_q.signal, phewas_q.analysis, phewas_q.pval, phewas_q.beta, 
             gwas_results_top_hits.pos_index AS th_pos, 
             gwas_results_top_hits.pos_start, 
             gwas_results_top_hits.pos_end
        FROM phewas_q 
        INNER JOIN gwas_results_top_hits
        USING (chrom, analysis)
        WHERE pos >= pos_start AND pos <= pos_end')
  
  phewas_marg <- 
    sqlWrapper(dbc, phewas_marg_sql, uniq = FALSE, zrok = TRUE) %>% 
    dplyr::as_tibble() %>% 
    dplyr::select(-pos_start, -pos_end)
  
  tictoc::toc(log = TRUE, quiet = TRUE)
  gtx_debug("int_ht_phewas | PheWAS on marginal GWAS results: {tictoc::tic.log(format = TRUE)}.")
  
  # --- PheWAS on conditional GWAS results
  gtx_info("Starting PheWAS of CLEO GWAS results.")
  tictoc::tic.clearlog();
  tictoc::tic();
  
  phewas_cleo_sql <- glue::glue('
    WITH 
    -- Define vars_q
    vars_q AS ({vars_sql}),
    -- Define gwas_q
    gwas_q AS ({cleo_sql})
    -- PheWAS by inner joining gwas_results_cond and variants
    SELECT vars_q.input, gwas_q.chrom, gwas_q.pos, gwas_q.ref, gwas_q.alt, 
           gwas_q.signal, gwas_q.analysis, gwas_q.pval_cond as pval, gwas_q.beta_cond as beta
      FROM gwas_q
      INNER JOIN vars_q
      USING (chrom, pos, ref, alt)')
  
  phewas_cleo <- 
    sqlWrapper(dbc, phewas_cleo_sql, uniq = FALSE, zrok = TRUE) %>% 
    dplyr::as_tibble() %>% 
    dplyr::mutate(th_pos = pos)
  
  tictoc::toc(log = TRUE, quiet = TRUE)
  gtx_debug("int_ht_phewas | PheWAS on CLEO GWAS results: {tictoc::tic.log(format = TRUE)}.")
  
  # Union the two data streams and return
  return(dplyr::union(phewas_marg, phewas_cleo))
}

#' int_ht_regional_context
#' 
#' High-Throughput regional_context.
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param input data.frame with chrom, pos, analysis, signal; likely & suggested input from int_ht_phewas. Other cols will be kept. 
#' @param chrom chromosome (Mandatory)
#' @param pos position (Mandatory)
#' @param analysis analysis ID (Mandatory)
#' @param signal CLEO signal OR NA. (Mandatory)
#' @param ... Params to pass to regionplot.data (e.g. surround = 1e6)
#' @param dbc [Default = getOption("gtx.dbConnection", NULL)] gtxconnection. 
#' @return tibble 
int_ht_regional_context <- function(input, chrom, pos, analysis, signal, 
                                    cpu = 8, drop_cs = TRUE, 
                                    dbc = getOption("gtx.dbConnection", NULL), ...){
  # --- Verify inputs
  # confirm input params are declared
  if(missing(input) & ((missing(chrom) | missing(pos) | missing(analysis) | missing(signal)))){
    gtx_fatal_stop("int_ht_regional_context_analysis | missing input arguement(s).")
  }
  # if input param used, verify it has the proper cols
  if(!missing(input)){
    required_cols <- c("chrom", "pos", "analysis", "signal")
    if(!all(required_cols %in% (names(input)))){
      gtx_fatal_stop('int_ht_regional_context_analysis | input is missing required cols. ',
                     'Required cols include: {paste(required_cols, collapse = ", ")}')
    }
  }
  # else create input table
  else if(!(missing(analysis) & !missing(chrom) & !missing(pos) & missing(signal))){
    input = dplyr::tibble(analysis = analysis, chrom = chrom, pos = pos, ref = ref, alt = alt, signal = NA)
  } else {
    gtx_fatal_stop("int_ht_regional_context | unable to process input properly.") 
  }
  if(nrow(input) == 0){
    gtx_warn("int_ht_regional_context_analysis | input has no data.")
    return(NA)
  }
  
  # --- Setup multi-threading
  current_database <- sqlWrapper(dbc, 'SELECT current_database() AS current_table;') %>% pull(current_table)
  if(nrow(input) > 3 & cpu > 1){
    # Remove global variables for dbc 
    options("gtx.dbConnection" = NULL);
    gtxcache(disconnect = TRUE)
    DBI::dbDisconnect(dbc)
    future::plan(future::multisession, workers = as.integer(cpu))
  } else {
    future::plan(future::sequential)
  }
  # Split input into uniq cred sets that we need to gather
  cs2get <- input %>% dplyr::distinct(analysis, chrom, th_pos, signal)
  
  # Get cred sets
  gtx_info("Starting regional context analysis for all GWAS results.")
  cs2get_ret <- 
    cs2get %>% 
    dplyr::mutate(cs = 
                    furrr::future_pmap(list(analysis, chrom, th_pos, signal), 
                                       ~int_ht_rpd_wrapper(analysis = ..1, chrom = ..2, pos = ..3, 
                                                           signal = ..4, use_database = current_database)))
    
  cs2get_ret <- 
    dplyr::inner_join(input, cs2get_ret, by = c("analysis", "chrom", "th_pos", "signal")) %>%
    dplyr::select(-th_pos)
  
  # Annotate critical data (e.g. variant in the cred set)
  ret <- 
    cs2get_ret %>% 
    dplyr::mutate(cs_annots = furrr::future_pmap(list(cs, chrom, pos, ref, alt), 
                  ~int_ht_regional_annot(cs = ..1, chrom = ..2, pos = ..3, 
                                         ref = ..4, alt = ..5)))
  
  if(isTRUE(drop_cs)){
    ret <- ret %>% dplyr::select(-cs)
  } else {
    ret <- 
      ret %>% 
      dplyr::mutate(cs = purrr::pmap(list(cs, chrom, pos),
                                     ~dplyr::filter(..1, cs_signal == TRUE | (chrom == ..2 & pos == ..3))))
  }
  
  # reconnect to DB
  if(nrow(input) > 3 & cpu > 1){
    options("gtx.dbConnection" = NULL);
    gtxcache(disconnect = TRUE)
    DBI::dbDisconnect(dbc)
    Sys.sleep(1)
    gtxconnect(use_database = current_database, cache = TRUE)
    Sys.sleep(1)
  }
  
  return(ret)
}

#' int_ht_rpd_wrapper
#' 
#' int_ht_rpd_wrapper - internal high-throughput RegionPlot.Data wrapper.
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param chrom chromosome
#' @param pos position
#' @param analysis analysis ID
#' @param signal CLEO signal or NA
#' @param use_database Database for gtxconnect
#' @param ... Additional arguements for regionplot.data
#' @examples 
#' ret <- tibble(analysis = "GSK500KV3lin_HES_Asthma", chrom = "9", pos = 6255967, signal = 4) %>%  
#'        mutate(cs = future_pmap(list(analysis, chrom, pos, signal), 
#'                                ~int_ht_rpd_wrapper(analysis = ..1, chrom = ..2, pos = ..3, 
#'                                                         signal = ..4, use_database = config_db())))
#' @return tibble with the credible set AND the variant queried
int_ht_rpd_wrapper <- function(analysis, chrom, pos, signal, use_database, ...){
  if(missing(chrom) | missing(pos) | missing(analysis) | missing(signal) | missing(use_database)){
    gtx_fatal_stop("int_ht_rpd_wrapper | missing input arguement(s).")
  }
  
  if(is.null(getOption("gtx.dbConnection", NULL))){
    gtx_debug("int_ht_rpd_wrapper | Establishing gtxconnection.")
    gtxconnect(use_database = use_database, cache = TRUE)
    gtxdbcheck(getOption("gtx.dbConnection"))
  } else {
    gtx_debug("int_ht_rpd_wrapper | Using pre-established gtxconnection.")
  }
  
  if(is.na(signal)){
    ret <- 
      regionplot.data(analysis = analysis, chrom = chrom, pos = pos, style = 'signal', ...) %>% 
      dplyr::as_tibble() %>% 
      dplyr::select(-freq, -rsq)
  } else if(!is.na(signal)){
    ret <- 
      # TODO look into using fm_cleo.data instead of regionplot.data here. 
      regionplot.data(analysis = analysis, chrom = chrom, pos = pos, style = 'signals', signal = signal, ...) %>% 
      dplyr::as_tibble() %>% 
      dplyr::filter(signal == !! signal) %>%       
      dplyr::select(-pp_signal, -cs_signal) %>% 
      dplyr::rename(pp_signal = pp_cleo, cs_signal = cs_cleo)
  } else {
    gtx_fatal_stop("int_ht_rpd_wrapper | Unable to determine type of signal to process.")
  }

  return(ret)
}

#' int_ht_regional_annot
#' 
#' High-Throughput regional_contex annot of cs.
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param input data.frame with chrom, pos, analysis, signal; likely & suggested input from int_ht_phewas. Other cols will be kept. 
#' @param chrom chromosome (Mandatory)
#' @param pos position (Mandatory)
#' @param analysis analysis ID (Mandatory)
#' @param signal CLEO signal OR NA. (Mandatory)
#' @param ... Params to pass to regionplot.data (e.g. surround = 1e6)
#' @param dbc [Default = getOption("gtx.dbConnection", NULL)] gtxconnection. 
#' @return tibble 
int_ht_regional_annot <- function(cs, chrom, pos, ref, alt){
  if(missing(cs) | missing(chrom) | missing(pos) | missing(ref) | missing(alt)){
    gtx_fatal_stop("int_ht_regional_annot missing input param.")
  }
  
  # queried variant in the cred set
  var_in_cs <-
    case_when(dplyr::filter(cs, chrom == !! chrom & pos == !! pos & ref == !! ref & alt == !! alt &
                              cs_signal == TRUE) %>% nrow(.) >= 1 ~ TRUE,
              TRUE ~ FALSE)

  # queried variant posterior probability
  var_pp <-
    case_when(dplyr::filter(cs, chrom == !! chrom & pos == !! pos & ref == !! ref & alt == !! alt) %>% nrow(.) >= 1 ~
                dplyr::filter(cs, chrom == !! chrom & pos == !! pos & ref == !! ref & alt == !! alt) %>% dplyr::pull(pp_signal),
              TRUE ~ as.double(NA))

  n_in_cs = dplyr::filter(cs, cs_signal == TRUE) %>% nrow()

  # cred set top hit posterior probability
  cs_th_pp <- max(cs$pp_signal)
  cs_th_data <-
    dplyr::top_n(cs, 1, pp_signal) %>%
    dplyr::select(pos, ref, alt, pval) %>%
    dplyr::rename(th_pos = pos, th_ref = ref, th_alt = alt, th_pval = pval)

  return(dplyr::tibble(var_in_cs, var_pp, n_in_cs, 
                       th_pos = cs_th_data$th_pos,
                       th_ref = cs_th_data$th_ref,
                       th_alt = cs_th_data$th_alt,
                       th_pval = cs_th_data$th_pval))
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
    
    input <- 
      input %>% 
      dplyr::mutate(where_clause = purrr::map_chr(input, ~gtxwhere(hgncid = .))) %>% 
      dplyr::mutate(where_clause = glue::glue("({where_clause})"))
    
    select_statement <- glue::glue('
      WITH 
        t1 as (SELECT * FROM genes 
               WHERE {glue::glue_collapse(input$where_clause, sep = " OR ")})
      SELECT hgncid AS input, vep.chrom, vep.pos, vep.ref, vep.alt 
        FROM vep
        INNER JOIN t1
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

# __END__
# # queried variant in the cred set
# var_in_cs <- 
#   case_when(dplyr::filter(cs, chrom == chrom & pos == pos &
#                             ref == ref & alt == alt & 
#                             cs_signal == TRUE) %>% nrow(.) >= 1 ~ TRUE, 
#             TRUE ~ FALSE)
# 
# # queried variant posterior probability
# var_pp <- 
#   case_when(dplyr::filter(cs, chrom == chrom & pos == pos & ref == ref & alt == alt) %>% nrow(.) >= 1 ~ 
#               dplyr::filter(cs, chrom == chrom & pos == pos & ref == ref & alt == alt) %>% dplyr::pull(pp_signal),
#             TRUE ~ NA_complex_)
# 
# n_in_cs = dplyr::filter(cs, cs_signal == TRUE) %>% nrow()
# 
# # cred set top hit posterior probability
# cs_th_pp <- max(cs$pp_signal)
# cs_th_data <- 
#   dplyr::top_n(cs, 1, pp_signal) %>% 
#   dplyr::select(pos, ref, alt, pval) %>% 
#   dplyr::rename(th_pos = pos, th_ref = ref, th_alt = alt, th_pval = pval) 
# 
# return(dplyr::tibble(var_in_cs, var_pp, n_in_cs, cs_th_pp, cs_th_data))


# ret <- 
#   input %>% 
#   dplyr::mutate(var_in_cs  = purrr::pmap_lgl(list(cs, chrom, pos, ref, alt), 
#                                              ~dplyr::filter(..1, chrom == ..2 & pos == ..3 & ref == ..4 & 
#                                                               alt == ..5 & cs_signal == TRUE) %>% 
#                                                nrow(.) >= 1)) %>% 
#   # queried variant posterior probability
#   dplyr::mutate(var_pp     = purrr::pmap(list(cs, chrom, pos, ref, alt), 
#                                          ~dplyr::filter(..1, chrom == ..2 & pos == ..3 & 
#                                                           ref == ..4 & alt == ..5 ) %>% 
#                                            dplyr::pull(pp_signal))) %>% 
#   # tidyr::unnest(var_pp, .drop = FALSE) %>% 
#   # number of snps in the cred set
#   dplyr::mutate(n_in_cs    = purrr::map_dbl(cs, ~dplyr::filter(., cs_signal == TRUE) %>% nrow)) %>% 
#   # cred set top hit posterior probability
#   dplyr::mutate(cs_th_pp   = purrr::map_dbl(cs, ~max(.$pp_signal))) %>% 
#   dplyr::mutate(cs_th_data = purrr::map(cs, ~dplyr::top_n(., 1, pp_signal) %>% 
#                                           dplyr::select(pos, ref, alt, pval))) %>% 
#   dplyr::select(-cs, -chrom, -pos, -ref, -alt) %>%
#   tidyr::unnest(cs_th_data, .drop = FALSE) %>% 
#   # dplyr::rename(th_pos = pos, th_ref = ref, th_alt = alt, th_pval = pval) 