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
#' @param include_negative_results [Default = FALSE] Return non-nearest gene within surround distance.
#' @param surround [Default = 1e5] Max distance (bp) to non-nearest gene. 
#' @param ignore_ukb_neale [Default = TRUE] TRUE = ignore Neale UKB GWAS in PheWAS
#' @param ignore_ukb_cane [Default = TRUE] TRUE = ignore Canelaxandri UKB GWAS in PheWAS
#' @param ignore_qtls [Default = TRUE] TRUE = ignore QTLs in PheWAS (gwas_results with entity)
#' @param dbc [Default = getOption("gtx.dbConnection", NULL)] gtxconnection. 
#' @return tibble (data.frame)
nearest_gene <- function(analysis, hgnc, hgncid, ensemblid, 
                         rs, rsid, chrom, pos,
                         protein_coding_only = TRUE,
                         calc_th = FALSE, gwas_args,
                         include_negative_results = FALSE,
                         surround = 1e5,
                         gwas_pval_le = 5e-8,
                         ignore_ukb_neale = TRUE, 
                         ignore_ukb_cane  = TRUE,
                         ignore_qtls      = TRUE, 
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
      int_nearest_gene_gene(hgncid    = hgncid, 
                            ensemblid = ensemblid, 
                            include_negative_results = include_negative_results,
                            surround     = surround, 
                            gwas_pval_le = gwas_pval_le,
                            protein_coding_only = protein_coding_only,
                            th_gwas_where_clause = th_gwas_where_clause, 
                            dbc = dbc)
  }
  # Order results
  nearest_genes <- nearest_genes %>% arrange(chrom, pos, signal, analysis)
  
  # return
  return(nearest_genes)
}

int_nearest_gene_gwas <- function(analysis, calc_th = FALSE, gwas_args,
                                  include_negative_results = FALSE,
                                  surround = 1e5,
                                  th_gwas_where_clause){
  
}

int_nearest_gene_gene <- function(hgncid, ensemblid,
                                  include_negative_results = FALSE,
                                  surround = 1e5, gwas_pval_le = 5e-8,
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
  
  genes_res <- sqlWrapper(dbc, genes_sql)
  # We will need to limit our searches across gwas_results_top_hits by chrom based on the input
  genes_chroms <- genes_res %>% dplyr::distinct(chrom)
  genes_chrom_where_clause <- 
    glue::glue('(',
         glue::glue_collapse(glue::glue("(chrom=\"{genes_chroms$chrom}\")"), sep = " OR "),
         ')')
  
  # Build SQL to query and join ALL genes within ${surround} of all marginal GWAS top hits based on chrom
  gtx_debug("int_nearest_gene_gene | JOIN gwas_results_top_hits with genes.")
  marg_nearest_sql <- glue('
    WITH 
      gwas_q AS 
        (SELECT analysis, chrom, cast(NULL as integer) as signal, pos_index AS pos, ref_index AS ref,
           alt_index AS alt, pval_index AS pval, beta_index AS beta
         FROM gwas_results_top_hits 
         WHERE {th_gwas_where_clause} AND
           {genes_chrom_where_clause} AND 
           pval_index <= {gwas_pval_le}),
      genes_q AS 
        (SELECT ensemblid, hgncid, pos_start, pos_end, strand, tss, genetype, chrom 
         FROM t_genes 
         WHERE {genes_chrom_where_clause})
    SELECT analysis, gwas_q.chrom, signal, pos, ref, alt, pval, beta, 
           ensemblid, hgncid, genetype, pos_start, pos_end, strand, tss, ABS(genes_q.tss - gwas_q.pos) AS dist2tss
      FROM gwas_q
      INNER JOIN genes_q
      ON (gwas_q.chrom = genes_q.chrom)
      WHERE ABS(genes_q.tss - gwas_q.pos) <= {surround}')
  
  marg_near = sqlWrapper(dbc, marg_nearest_sql, uniq = FALSE, zrok = FALSE)
  
  # Build SQL to query and join ALL genes within ${surround} of all CLEO GWAS JOINT hits based on chrom
  gtx_debug("int_nearest_gene_gene | JOIN gwas_results_joint with genes.")
  cleo_nearest_sql <- glue('
    WITH 
      gwas_q AS 
        (SELECT analysis, signal, chrom, pos, ref, alt, pval_joint AS pval, beta_joint AS beta
         FROM gwas_results_joint 
         WHERE {th_gwas_where_clause} AND
           {genes_chrom_where_clause} AND
           pval_joint <= {gwas_pval_le}),
      genes_q AS 
        (SELECT ensemblid, hgncid, pos_start, pos_end, strand, tss, genetype, chrom 
         FROM t_genes 
         WHERE {genes_chrom_where_clause})
    SELECT analysis, gwas_q.chrom, signal, pos, ref, alt, pval, beta, 
           ensemblid, hgncid, genetype, pos_start, pos_end, strand, tss, ABS(genes_q.tss - gwas_q.pos) AS dist2tss
      FROM gwas_q
      INNER JOIN genes_q
      ON (gwas_q.chrom = genes_q.chrom)
      WHERE ABS(genes_q.tss - gwas_q.pos) <= {surround}')
  
  cleo_near = sqlWrapper(dbc, cleo_nearest_sql, uniq = FALSE, zrok = FALSE)
  
  # Merge multiple data streams
  if(nrow(marg_near) >= 1 & nrow(cleo_near) >= 1){
    gtx_debug("Merging marginal and CLEO near genes.")
    all_near <- union(marg_near, cleo_near)
  } else if(nrow(marg_near) >= 1 & nrow(cleo_near) == 0){
    gtx_debug("Only using marginal near genes.")
    all_near <- marg_near
  } else if(nrow(marg_near) == 0 & nrow(cleo_near) >= 1){
    gtx_debug("Only using CLEO near genes.")
    all_near <- cleo_near
  } else {
    gtx_fatal_stop("Did not find ANY genes near GWAS top hits. Something has likely gone terribly wrong.")
  }
  
  all_near <- 
    all_near %>% 
    as_tibble() %>% 
    mutate(dist2tss = as.integer(dist2tss))
  
  # Remove non-protein_coding transcripts if param = TRUE
  if(protein_coding_only == TRUE){
    gtx_warn("Filtering for protein_coding genes only.")
    all_near <- all_near %>% filter(genetype == "protein_coding")
  }
  
  # From all the near genes, ID the nearest genes 
  all_nearest <- int_calc_nearest_from_near(input = all_near, include_negative_results = include_negative_results)
  
  # Mark loci with input genes and keep entire loci containing input hits.
  if(!missing(hgncid)){
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

int_nearest_gene_chrom_pos <- function(chrom, pos,
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

int_calc_nearest_from_near <- function(input, include_negative_results = FALSE){
  if(missing(input)){
    gtx_fatal_stop("Missing input into int_calc_nearest_from_near.")
  }
  
  if(!missing(input)){
    required_cols <- c("analysis", "signal", "chrom", "pos", "ensemblid")
    if(!all(required_cols %in% (names(input)))){
      gtx_fatal_stop('int_ht_regional_context_analysis | input is missing required cols. ',
                     'Required cols include: {paste(required_cols, collapse = ", ")}')
    }
  }
  
  ret <- 
    input %>% 
    group_by(analysis, chrom, signal, pos, ref, alt) %>% 
    mutate(min_dist2tss = min(dist2tss, na.rm = TRUE)) %>% 
    mutate(nearest_gene = min_dist2tss == dist2tss) %>% 
    select(-min_dist2tss) %>% 
    ungroup()
  
  if(include_negative_results == TRUE){
    # This includes the "negative" non-nearest genes calculted within {surround}
    return(ret) 
  } else if(include_negative_results == FALSE) {
    # Filter for only positive nearest gene results
    return(ret %>% filter(nearest_gene == TRUE))
  } else {
    gtx_fatal_stop("include_negative_results is not TRUE/FALSE.")
  }
}
