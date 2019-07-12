#' nearest_gene
#' 
#' Calculate nearest_gene V2G to identify targets~indications.
#' 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @param analysis GWAS analysis ID
#' @param ensemblid Ensembl gene ID.
#' @param hgnc Analyze all moderate/high impact VEP variants overlapping gene position.
#' @param rs RSID(s) to analyze as a string or vector
#' @param chrom chromosome(s) to analyze as a string or vector
#' @param pos positions(s) to analyze as a single or vector
#' @param nearest_method  [Default = 'nearest_tss'] String(s) specifying the method for calculating nearest gene. 
#' {'nearest_tss', 'nearest_dist_pruned', 'nearest_localized'}
#' @param localized_dist_ge [Default = 1e5] For nearest_method 'nearest_localized', this is the distance that the second nearest gene must be from the TH. 
#' @param include_negative_results [Default = FALSE] Return non-nearest gene within surround distance.
#' @param surround [Default = 1e6] Max distance (bp) to non-nearest gene. 
#' @param ignore_ukb_neale [Default = TRUE] TRUE = ignore Neale UKB GWAS in PheWAS
#' @param ignore_ukb_cane [Default = TRUE] TRUE = ignore Canelaxandri UKB GWAS in PheWAS
#' @param ignore_qtls [Default = TRUE] TRUE = ignore QTLs in PheWAS (gwas_results with entity)
#' @param dbc [Default = getOption("gtx.dbConnection", NULL)] gtxconnection. 
#' @return tibble (data.frame)
nearest_gene <- function(analysis, hgnc, hgncid, ensemblid, 
                         rs, rsid, chrom, pos,
                         nearest_method = c('nearest_tss', 'nearest_dist_pruned', 'nearest_localized'),
                         localized_dist_ge = 1e5,
                         # max_dist2tss = NULL, #TODO
                         protein_coding_only = TRUE,
                         include_negative_results = FALSE,
                         surround = 1e6,
                         gwas_pval_le = 5e-8,
                         calc_th = FALSE, gwas_args,
                         ignore_ukb_neale = TRUE, 
                         ignore_ukb_cane  = TRUE,
                         dbc = getOption("gtx.dbConnection", NULL)){
  # --- Check gtx db connection
  gtx_debug("regional_context_analysis | validating gtx DB connection.")
  gtxdbcheck(dbc)
  
  th_gwas_where_clause <- 
    dplyr::tribble(~limits,            ~value,           ~sql_statement,
                   "ignore_NULL",      TRUE,             'analysis IS NOT NULL',
                   "ignore_ukb_neale", ignore_ukb_neale, '!regexp_like(analysis, "neale")',
                   "ignore_ukb_cane",  ignore_ukb_cane,  '!regexp_like(analysis, "canelaxandri17")') %>% 
           #"ignore_qtls",      ignore_qtls,      'entity IS NULL')
    dplyr::filter(value == TRUE) %>% 
    dplyr::pull(sql_statement) %>% 
    glue::glue_collapse(sep = " AND ")
  
  # do stuff
  if(!missing(hgncid) | !missing(ensemblid)){
    nearest_genes <- 
      int_nearest_gene_instr_on_gene(hgncid    = hgncid, ensemblid = ensemblid, 
                                      nearest_method = nearest_method, localized_dist_ge = localized_dist_ge,
                                      include_negative_results = include_negative_results,
                                      surround = surround, gwas_pval_le = gwas_pval_le,
                                      protein_coding_only  = protein_coding_only,
                                      th_gwas_where_clause = th_gwas_where_clause, 
                                      dbc = dbc)
  }
  
  
  # Order results
  nearest_genes <- 
    nearest_genes %>% 
    arrange(chrom, pos_index, signal, analysis) %>% 
    select(-start_gene, -end_gene, -strand_gene) %>% 
    select(-n_gene_overlap_th, -matches("nearest"), 
           everything(), 
           n_gene_overlap_th, matches("nearest"))
  
  # return
  return(nearest_genes)
}

int_nearest_gene_instr_on_gene <- function(hgncid, ensemblid,
                                           nearest_method,
                                           localized_dist_ge,
                                           include_negative_results = FALSE,
                                           surround = 1e6, gwas_pval_le = 5e-8,
                                           th_gwas_where_clause, protein_coding_only = FALSE, 
                                           dbc = getOption("gtx.dbConnection", NULL)){
  # Verify inputs
   if(missing(hgncid) & missing(ensemblid)){
     gtx_fata_stop("Missing input param(s). Use hgncid, ensemblid, chrom, pos. ")
   }
  
  if(!missing(hgncid)){
    input <- 
      dplyr::tibble(input = hgncid) %>% 
      dplyr::mutate(where_clause = purrr::map_chr(input, ~gtxwhere(hgncid = .)))
  } else if (!missing(ensemblid)) {
    input <- 
      dplyr::tibble(input = ensemblid) %>% 
      dplyr::mutate(where_clause = purrr::map_chr(input, ~gtxwhere(ensemblid = .)))
  } 
  # Build SQL to get genes of interest
  genes_sql <- glue::glue('
    SELECT * 
      FROM t_genes 
      WHERE {glue::glue_collapse(input$where_clause, sep = " OR ")}')
  
  genes_res <- sqlWrapper(dbc, genes_sql, uniq = FALSE, zrok = FALSE)
  
  # We will need to limit our searches across gwas_results_top_hits by chrom based on the input
  genes_chroms <- 
    genes_res %>% 
    dplyr::group_by(chrom) %>% 
    dplyr::summarize(chrom_min_pos_start = max(pos_start), chrom_max_pos_end = max(pos_end)) %>% 
    dplyr::mutate(chrom_min_pos_start = chrom_min_pos_start - surround) %>% 
    dplyr::mutate(chrom_max_pos_end = chrom_max_pos_end + surround) %>% 
    dplyr::mutate(th_where_clause = glue::glue('chrom = "{chrom}" AND pos_index >= {chrom_min_pos_start} \\
                                                AND pos_index <= {chrom_max_pos_end}'))
    
  th_cp_where_clause <- 
    glue::glue('(',
         glue::glue_collapse(glue::glue("({genes_chroms$th_where_clause})"), sep = " OR "),
         ')')
  
  genes_tbl_c_where_clause <- 
    glue::glue('(',
               glue::glue_collapse(glue::glue('(chrom="{distinct(genes_chroms, chrom) %>% pull()}")'), sep = " OR "),
               ')')
  
  # Build SQL to query and join ALL genes within ${surround} of all marginal GWAS top hits based on chrom
  gtx_debug("int_nearest_gene_instr_on_gene | JOIN gwas_results_top_hits with genes.")
  marg_nearest_sql <- glue::glue('
    WITH 
      gwas_q AS 
        (SELECT analysis, chrom, cast(NULL as integer) as signal, pos_index, ref_index,
           alt_index, pval_index, beta_index, pos_start AS start_index, pos_end AS end_index
         FROM gwas_results_top_hits 
         WHERE {th_gwas_where_clause} 
           AND {th_cp_where_clause} 
           AND pval_index <= {gwas_pval_le}),
      genes_q AS 
        (SELECT ensemblid, hgncid, pos_start AS start_gene, pos_end AS end_gene, 
                strand AS strand_gene, tss AS tss_gene, genetype, chrom 
         FROM t_genes 
         WHERE {genes_tbl_c_where_clause})
    SELECT analysis, gwas_q.chrom, signal, pos_index, ref_index, alt_index, 
           pval_index, beta_index, start_index, end_index, 
           ensemblid, hgncid, genetype, start_gene, end_gene, strand_gene, 
           tss_gene, ABS(genes_q.tss_gene - gwas_q.pos_index) AS dist2tss
    FROM gwas_q
    INNER JOIN genes_q
    ON (gwas_q.chrom = genes_q.chrom)
    WHERE ABS(genes_q.tss_gene - gwas_q.pos_index) <= {surround}')
  
  marg_near <- sqlWrapper(dbc, marg_nearest_sql, uniq = FALSE, zrok = FALSE) %>% as_tibble()
  
  # Build SQL to query and join ALL genes within ${surround} of all CLEO GWAS JOINT hits based on chrom
  gtx_debug("int_nearest_gene_instr_on_gene | JOIN gwas_results_joint with genes.")
  th_cp_where_clause <- stringr::str_replace_all(th_cp_where_clause, "pos_index", "pos")
  cleo_nearest_sql <- glue::glue('
    WITH 
      gwas_q AS 
        (SELECT analysis, chrom, signal, pos AS pos_index, ref AS ref_index, 
                alt AS alt_index, pval_joint AS pval_index, beta_joint AS beta_index, 
                cast(NULL as integer) as start_index, cast(NULL as integer) as end_index
         FROM gwas_results_joint 
         WHERE {th_gwas_where_clause} 
          AND {th_cp_where_clause} 
          AND pval_joint <= {gwas_pval_le}),
      genes_q AS 
        (SELECT ensemblid, hgncid, pos_start AS start_gene, pos_end AS end_gene, 
                strand AS strand_gene, tss AS tss_gene, genetype, chrom 
         FROM t_genes 
         WHERE {genes_tbl_c_where_clause})
    SELECT analysis, gwas_q.chrom, signal, pos_index, ref_index, alt_index, 
           pval_index, beta_index, start_index, end_index, 
           ensemblid, hgncid, genetype, start_gene, end_gene, strand_gene, 
           tss_gene, ABS(genes_q.tss_gene - gwas_q.pos_index) AS dist2tss
    FROM gwas_q
    INNER JOIN genes_q
    ON (gwas_q.chrom = genes_q.chrom)
    WHERE ABS(genes_q.tss_gene - gwas_q.pos_index) <= {surround}')
  
  cleo_near = sqlWrapper(dbc, cleo_nearest_sql, uniq = FALSE, zrok = FALSE) %>% as_tibble()
  
  # Merge multiple data streams
  if(nrow(marg_near) >= 1 & nrow(cleo_near) >= 1){
    gtx_debug("Merging marginal and CLEO nearest genes.")
    all_near <- dplyr::union(marg_near, cleo_near)
  } else if(nrow(marg_near) >= 1 & nrow(cleo_near) == 0){
    gtx_debug("Only using marginal nearest genes.")
    all_near <- marg_near
  } else if(nrow(marg_near) == 0 & nrow(cleo_near) >= 1){
    gtx_debug("Only using CLEO nearest genes.")
    all_near <- cleo_near
  } else {
    gtx_fatal_stop("Did not find ANY genes near GWAS top hits. Something has likely gone terribly wrong.")
  }
  
  # Remove non-protein_coding transcripts if param = TRUE
  if(protein_coding_only == TRUE){
    gtx_warn("Filtering for protein_coding genes only.")
    all_near <- all_near %>% dplyr::filter(genetype == "protein_coding")
  }
  
  # From all the near genes, ID the nearest genes 
  all_nearest <- int_calc_nearest_from_near(input = all_near, 
                                            nearest_method = nearest_method, 
                                            localized_dist_ge = localized_dist_ge,
                                            include_negative_results = include_negative_results)
  
  # Mark loci with input genes and keep entire loci containing input hits.
  if(!missing(hgncid)){
    ret <- 
      inner_join(input %>% 
                 select(-where_clause) %>% 
                 rename(hgncid = input),
               all_nearest %>% 
                 dplyr::filter_at(.vars = dplyr::vars(dplyr::matches("nearest")), 
                                  .vars_predicate = dplyr::any_vars(. == TRUE)), 
               by = "hgncid") %>% 
      distinct(analysis, chrom, signal, pos_index, ref_index, alt_index) %>% 
      mutate(input_locus = TRUE) %>% 
      left_join(all_nearest, 
                by = c("analysis", "chrom", "signal", "pos_index", "ref_index", "alt_index")) %>% 
      filter(input_locus == TRUE) %>% 
      select(-input_locus)
  } else if(!missing(ensemblid)){
    ret <- 
      inner_join(input %>% 
                   select(-where_clause) %>% 
                   rename(hgncid = input),
                 all_nearest %>% filter(nearest_gene == TRUE), 
                 by = "hgncid") %>% 
      distinct(analysis, chrom, signal, pos_index, ref_index, alt_index) %>% 
      mutate(input_locus = TRUE) %>% 
      left_join(all_nearest, 
                by = c("analysis", "chrom", "signal", "pos_index", "ref_index", "alt_index")) %>% 
      filter(input_locus == TRUE) %>% 
      select(-input_locus)
  } else {
    gtx_fatal_stop("Uncertain how to filter final results.")
  }
  
  return(ret)
}

int_nearest_gene_instr_on_chrom_pos <- function(chrom, pos,
                                       include_negative_results = FALSE,
                                       surround = 1e5, gwas_pval_le = 5e-8,
                                       th_gwas_where_clause, protein_coding_only = FALSE, 
                                       dbc = getOption("gtx.dbConnection", NULL)){
  # Verify inputs
  if(missing(chrom) & missing(pos)){
    gtx_fata_stop("Missing input param(s). Use chrom & pos.")
  }
  
  input <- 
    dplyr::tibble(chrom = chrom, pos = pos) %>% 
    dplyr::mutate(input = glue::glue("chr{chrom}_{pos}")) %>% 
    dplyr::mutate(where_clause = purrr::map2_chr(chrom, ~gtxwhere(chrom = .)))
  
  # Build SQL to get genes of interest
  genes_sql <- glue::glue('
                          SELECT * 
                          FROM t_genes 
                          WHERE {glue::glue_collapse(input$where_clause, sep = " OR ")}')
  
  if(is.null(all_near)){
    gtx_fatal_stop("Did not find ANY genes near GWAS top hits. Something has likely gone terribly wrong.")
  } else {
    all_near <- 
      all_near %>% 
      as_tibble() %>% 
      mutate(dist2tss = as.integer(dist2tss))
  }
  
  if(protein_coding_only == TRUE){
    gtx_warn("Filtering for protein_coding genes only.")
    all_near <- all_near %>% filter(genetype == "protein_coding")
  }
  
  all_nearest <- int_calc_nearest_from_near(input = all_near, include_negative_results = include_negative_results)
  
  if(!missing(hgncid)){
    # Mark loci with input genes and keep entire loci containing input hits.
    ret <- 
      inner_join(input %>% 
                   select(-where_clause) %>% 
                   rename(hgncid = input),
                 all_nearest %>% filter(nearest_gene == TRUE), 
                 by = "hgncid") %>% 
      distinct(analysis, chrom, signal, pos, ref, alt) %>% 
      mutate(input_locus = TRUE) %>% 
      left_join(all_nearest, 
                by = c("analysis", "chrom", "signal", "pos", "ref", "alt")) %>% 
      filter(input_locus == TRUE) %>% 
      select(-input_locus)
    
    return(ret)
  } else if(!missing(ensemblid)){
    ret <- 
      inner_join(input %>% 
                   select(-where_clause) %>% 
                   rename(hgncid = input),
                 all_nearest %>% filter(nearest_gene == TRUE), 
                 by = "hgncid") %>% 
      distinct(analysis, chrom, signal, pos, ref, alt) %>% 
      mutate(input_locus = TRUE) %>% 
      left_join(all_nearest, 
                by = c("analysis", "chrom", "signal", "pos", "ref", "alt")) %>% 
      filter(input_locus == TRUE) %>% 
      select(-input_locus)
    
    return(ret)
  } else {
    gtx_fatal_stop("Uncertain how to filter final results.")
  }
}

int_nearest_gene_instr_on_gwas <- function(analysis, aba_gwas_results = FALSE, calc_th = FALSE, gwas_args,
                                           include_negative_results = FALSE,
                                           surround = 1e6){
  # Verify inputs
  if(missing(analysis)){
    gtx_fata_stop("Missing input param: analysis.")
  }
  
  if(!missing(analysis) & aba_gwas_results == FALSE){
    input <- 
      dplyr::tibble(input = analysis) %>% 
      dplyr::mutate(where_clause = glue::glue('analysis = "{analysis}"'))
    
    th_analysis_where_clause <-
      glue::glue('(',
                 glue::glue_collapse(glue::glue("({input$where_clause})"), sep = " OR "),
                 ')')
    
    th_chroms_sql <- 
      glue::glue('SELECT distinct(chrom) 
                  FROM gwas_results_top_hits
                  WHERE {th_analysis_where_clause}
                  AND pval_index IS NOT NULL')
    
    th_chroms_res <- sqlWrapper(dbc, th_chroms_sql, uniq = FALSE, zrok = FALSE)
    
  } else if(missing(analysis) & aba_gwas_results == TRUE){
    
  } else {
    gtx_fatal_stop("int_nearest_gene_instr_on_gwas | Unable to determine which GWAS to instrument on.")
  }
  
  
 
  # Build SQL to get genes of interest
  genes_sql <- glue::glue('
                          SELECT * 
                          FROM t_genes 
                          WHERE {glue::glue_collapse(input$where_clause, sep = " OR ")}')
  
  genes_res <- sqlWrapper(dbc, genes_sql, uniq = FALSE, zrok = FALSE)
  
  # We will need to limit our searches across gwas_results_top_hits by chrom based on the input
  genes_chroms <- 
    genes_res %>% 
    dplyr::group_by(chrom) %>% 
    dplyr::summarize(chrom_min_pos_start = max(pos_start), chrom_max_pos_end = max(pos_end)) %>% 
    dplyr::mutate(chrom_min_pos_start = chrom_min_pos_start - surround) %>% 
    dplyr::mutate(chrom_max_pos_end = chrom_max_pos_end + surround) %>% 
    dplyr::mutate(th_where_clause = glue::glue('chrom = "{chrom}" AND pos_index >= {chrom_min_pos_start} \\
                                               AND pos_index <= {chrom_max_pos_end}'))
  
  th_cp_where_clause <- 
    glue::glue('(',
               glue::glue_collapse(glue::glue("({genes_chroms$th_where_clause})"), sep = " OR "),
               ')')
  
  genes_tbl_c_where_clause <- 
    glue::glue('(',
               glue::glue_collapse(glue::glue('(chrom="{distinct(genes_chroms, chrom) %>% pull()}")'), sep = " OR "),
               ')')
  
  # Build SQL to query and join ALL genes within ${surround} of all marginal GWAS top hits based on chrom
  gtx_debug("int_nearest_gene_instr_on_gene | JOIN gwas_results_top_hits with genes.")
  marg_nearest_sql <- glue::glue('
                                 WITH 
                                 gwas_q AS 
                                 (SELECT analysis, chrom, cast(NULL as integer) as signal, pos_index, ref_index,
                                 alt_index, pval_index, beta_index, pos_start AS start_index, pos_end AS end_index
                                 FROM gwas_results_top_hits 
                                 WHERE {th_gwas_where_clause} 
                                 AND {th_cp_where_clause} 
                                 AND pval_index <= {gwas_pval_le}),
                                 genes_q AS 
                                 (SELECT ensemblid, hgncid, pos_start AS start_gene, pos_end AS end_gene, 
                                 strand AS strand_gene, tss AS tss_gene, genetype, chrom 
                                 FROM t_genes 
                                 WHERE {genes_tbl_c_where_clause})
                                 SELECT analysis, gwas_q.chrom, signal, pos_index, ref_index, alt_index, 
                                 pval_index, beta_index, start_index, end_index, 
                                 ensemblid, hgncid, genetype, start_gene, end_gene, strand_gene, 
                                 tss_gene, ABS(genes_q.tss_gene - gwas_q.pos_index) AS dist2tss
                                 FROM gwas_q
                                 INNER JOIN genes_q
                                 ON (gwas_q.chrom = genes_q.chrom)
                                 WHERE ABS(genes_q.tss_gene - gwas_q.pos_index) <= {surround}')
  
  marg_near <- sqlWrapper(dbc, marg_nearest_sql, uniq = FALSE, zrok = FALSE) %>% as_tibble()
  
  # Build SQL to query and join ALL genes within ${surround} of all CLEO GWAS JOINT hits based on chrom
  gtx_debug("int_nearest_gene_instr_on_gene | JOIN gwas_results_joint with genes.")
  th_cp_where_clause <- stringr::str_replace_all(th_cp_where_clause, "pos_index", "pos")
  cleo_nearest_sql <- glue::glue('
                                 WITH 
                                 gwas_q AS 
                                 (SELECT analysis, chrom, signal, pos AS pos_index, ref AS ref_index, 
                                 alt AS alt_index, pval_joint AS pval_index, beta_joint AS beta_index, 
                                 cast(NULL as integer) as start_index, cast(NULL as integer) as end_index
                                 FROM gwas_results_joint 
                                 WHERE {th_gwas_where_clause} 
                                 AND {th_cp_where_clause} 
                                 AND pval_joint <= {gwas_pval_le}),
                                 genes_q AS 
                                 (SELECT ensemblid, hgncid, pos_start AS start_gene, pos_end AS end_gene, 
                                 strand AS strand_gene, tss AS tss_gene, genetype, chrom 
                                 FROM t_genes 
                                 WHERE {genes_tbl_c_where_clause})
                                 SELECT analysis, gwas_q.chrom, signal, pos_index, ref_index, alt_index, 
                                 pval_index, beta_index, start_index, end_index, 
                                 ensemblid, hgncid, genetype, start_gene, end_gene, strand_gene, 
                                 tss_gene, ABS(genes_q.tss_gene - gwas_q.pos_index) AS dist2tss
                                 FROM gwas_q
                                 INNER JOIN genes_q
                                 ON (gwas_q.chrom = genes_q.chrom)
                                 WHERE ABS(genes_q.tss_gene - gwas_q.pos_index) <= {surround}')
  
  cleo_near = sqlWrapper(dbc, cleo_nearest_sql, uniq = FALSE, zrok = FALSE) %>% as_tibble()
  
  # Merge multiple data streams
  if(nrow(marg_near) >= 1 & nrow(cleo_near) >= 1){
    gtx_debug("Merging marginal and CLEO nearest genes.")
    all_near <- dplyr::union(marg_near, cleo_near)
  } else if(nrow(marg_near) >= 1 & nrow(cleo_near) == 0){
    gtx_debug("Only using marginal nearest genes.")
    all_near <- marg_near
  } else if(nrow(marg_near) == 0 & nrow(cleo_near) >= 1){
    gtx_debug("Only using CLEO nearest genes.")
    all_near <- cleo_near
  } else {
    gtx_fatal_stop("Did not find ANY genes near GWAS top hits. Something has likely gone terribly wrong.")
  }
  
  # Remove non-protein_coding transcripts if param = TRUE
  if(protein_coding_only == TRUE){
    gtx_warn("Filtering for protein_coding genes only.")
    all_near <- all_near %>% dplyr::filter(genetype == "protein_coding")
  }
  
  # From all the near genes, ID the nearest genes 
  all_nearest <- int_calc_nearest_from_near(input = all_near, 
                                            nearest_method = nearest_method, 
                                            localized_dist_ge = localized_dist_ge,
                                            include_negative_results = include_negative_results)
  
  # Mark loci with input genes and keep entire loci containing input hits.
  if(!missing(hgncid)){
    ret <- 
      inner_join(input %>% 
                   select(-where_clause) %>% 
                   rename(hgncid = input),
                 all_nearest %>% 
                   dplyr::filter_at(.vars = dplyr::vars(dplyr::matches("nearest")), 
                                    .vars_predicate = dplyr::any_vars(. == TRUE)), 
                 by = "hgncid") %>% 
      distinct(analysis, chrom, signal, pos_index, ref_index, alt_index) %>% 
      mutate(input_locus = TRUE) %>% 
      left_join(all_nearest, 
                by = c("analysis", "chrom", "signal", "pos_index", "ref_index", "alt_index")) %>% 
      filter(input_locus == TRUE) %>% 
      select(-input_locus)
  } else if(!missing(ensemblid)){
    ret <- 
      inner_join(input %>% 
                   select(-where_clause) %>% 
                   rename(hgncid = input),
                 all_nearest %>% filter(nearest_gene == TRUE), 
                 by = "hgncid") %>% 
      distinct(analysis, chrom, signal, pos_index, ref_index, alt_index) %>% 
      mutate(input_locus = TRUE) %>% 
      left_join(all_nearest, 
                by = c("analysis", "chrom", "signal", "pos_index", "ref_index", "alt_index")) %>% 
      filter(input_locus == TRUE) %>% 
      select(-input_locus)
  } else {
    gtx_fatal_stop("Uncertain how to filter final results.")
  }
  
  return(ret)
}


int_calc_nearest_from_near <- function(input, nearest_method = c('nearest_tss', 'nearest_dist_pruned', 'nearest_localized'), 
                                       include_negative_results = FALSE, localized_dist_ge = 1e5){
  if(missing(input)){
    gtx_fatal_stop("Missing input into int_calc_nearest_from_near.")
  }
  
  if(!missing(input)){
    # Make sure the input has all the required columns
    required_cols <- c("analysis", "signal", "chrom", "pos_index", "ref_index", "alt_index", 
                       "dist2tss", "start_index", "end_index", "start_gene", "end_gene")
    if(!all(required_cols %in% (names(input)))){
      gtx_fatal_stop('int_ht_regional_context_analysis | input is missing required cols. 
                      Required cols include: {paste(required_cols, collapse = ", ")}')
    }
  }
  if(!is.character(nearest_method)){
    gtx_fatal_stop("nearest_method is not a type.") #TODO Fix documentation here. 
  }
  
  input <- 
    input %>% 
    dplyr::group_by(analysis, chrom, signal, pos_index, ref_index, alt_index) %>% 
    dplyr::mutate(dist2tss = as.integer(dist2tss)) %>% 
    dplyr::arrange(dist2tss)
  
  if(any(nearest_method == "nearest_tss")){
    input <- 
      input %>% 
      dplyr::mutate(nearest_gene = (dist2tss == min(dist2tss, na.rm = TRUE)))
  } 
  
  if(any(nearest_method == "nearest_dist_pruned")){
    input <- 
      input %>% 
      dplyr::mutate(next_dist2tss = dplyr::lead(dist2tss)) %>% 
      dplyr::mutate(nearest_dist_pruned = (dist2tss == min(dist2tss, na.rm = TRUE) & next_dist2tss >= localized_dist_ge))
  } 
  
  if(any(nearest_method == "nearest_localized")){
    input <- 
      input %>% 
      dplyr::mutate(next_dist2tss = dplyr::lead(dist2tss)) %>% 
      tidyr::nest() %>% 
      dplyr::mutate(n_gene_overlap_th = purrr::map_dbl(data, 
        ~dplyr::filter(., (start_gene >= start_index & start_gene <= end_index) |
                             (end_gene >= start_index & end_gene <= end_index) |
                             (start_index < start_gene & end_index > end_gene)) %>% 
              dplyr::n_distinct(hgncid, na.rm = TRUE))) %>% 
      tidyr::unnest() %>% 
      dplyr::group_by(analysis, chrom, signal, pos_index, ref_index, alt_index) %>% 
      dplyr::mutate(nearest_localized = (dist2tss == min(dist2tss, na.rm = TRUE) & 
                                         next_dist2tss >= localized_dist_ge & 
                                         ((start_index >= start_gene & start_index <= end_gene) | 
                                          (end_index   >= start_gene & end_index   <= end_gene) | 
                                          (start_index >= start_gene & end_index   <= end_gene) |
                                          (start_index <= start_gene & end_index   >= end_gene)) & 
                                           n_gene_overlap_th == 1)) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(nearest_localized = replace(nearest_localized, !is.na(signal), NA)) %>% 
      dplyr::mutate(n_gene_overlap_th = replace(n_gene_overlap_th, !is.na(signal), NA))
    
  } 
  
  ret <- input %>% dplyr::ungroup()
  
  if(include_negative_results == FALSE) {
    # Filter for only positive nearest gene results
    ret <- ret %>% dplyr::filter_at(.vars = dplyr::vars(dplyr::matches("nearest")), 
                                    .vars_predicate = dplyr::any_vars(. == TRUE))
  }
  
  return(ret)
}

