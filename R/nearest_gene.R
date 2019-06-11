#' nearest_gene
#' 
#' Calculate nearest_gene V2G to identify targets~indications.
#' 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @param analysis GWAS analysis ID
#' @param hgnc Analyze all moderate/high impact VEP variants overlapping gene position.
#' @param rs RSID(s) to analyze as a string or vector
#' @param chrom chromosome(s) to analyze as a string or vector
#' @param pos positions(s) to analyze as a single or vector
#' @param include_negative_results [Default = FALSE] Return non-nearest gene within surround distance.
#' @param surround [Default = 1e6] Max distance to non-nearest gene
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
                         surround = 1e6,
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
}

int_nearest_gene_gwas <- function(analysis, calc_th = FALSE, gwas_args,
                                  include_negative_results = FALSE,
                                  surround = 1e6,
                                  th_gwas_where_clause){
  
}

# Combine this with chrom pos? #TODO
int_nearest_gene_gene <- function(hgncid, ensemblid, chrom, pos,
                                  include_negative_results = FALSE,
                                  surround = 1e6, th_gwas_where_clause){
   if(missing(hgncid) & missing(ensemblid) & (missing(chrom) | missing(pos))){
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
  
  genes_sql <- glue::glue('
    SELECT * 
      FROM t_genes 
      WHERE {glue::glue_collapse(input$where_clause, sep = " OR ")}')
  
  genes_res <- sqlWrapper(dbc, genes_sql)
  genes_chroms <- genes_res %>% dplyr::distinct(chrom)
  genes_chrom_where_clause <- 
    glue::glue('(',
         glue::glue_collapse(glue::glue("(chrom=\"{genes_chroms$chrom}\")"), sep = " OR "),
         ')')
  
  marg_nearest_sql <- glue('
    WITH 
      gwas_q AS 
        (SELECT analysis, cast(NULL as integer) as signal, chrom, pos_index AS pos, ref_index AS ref,
           alt_index AS alt, pval_index AS pval, beta_index AS beta
         FROM gwas_results_top_hits 
         WHERE {th_gwas_where_clause} AND
           {genes_chrom_where_clause}),
      genes_q AS 
        (SELECT ensemblid, hgncid, pos_start, pos_end, strand, tss, genetype, chrom 
         FROM t_genes 
         WHERE {genes_chrom_where_clause})
    SELECT *, ABS(genes_q.tss - gwas_q.pos) AS `dist2tss`
      FROM gwas_q
      INNER JOIN genes_q
      USING (chrom)
      WHERE ABS(genes_q.tss - gwas_q.pos) <= {surround}')
  
}
