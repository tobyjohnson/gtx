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
#' @param ignore_ukb_neale [Default = TRUE] TRUE = ignore Neale UKB GWAS in PheWAS
#' @param ignore_ukb_cane [Default = TRUE] TRUE = ignore Canelaxandri UKB GWAS in PheWAS
#' @param ignore_qtls [Default = TRUE] TRUE = ignore QTLs in PheWAS (gwas_results with entity)
#' @param dbc [Default = getOption("gtx.dbConnection", NULL)] gtxconnection. 
#' @return tibble (data.frame)
nearest_gene <- function(analysis, hgnc, hgncid, ensemblid, 
                         rs, rsid, chrom, pos, ref, alt, 
                         calc_th = FALSE, gwas_args,
                         include_negative_results = FALSE,
                         ignore_ukb_neale = TRUE, 
                         ignore_ukb_cane  = TRUE,
                         ignore_qtls      = TRUE, 
                         dbc = getOption("gtx.dbConnection", NULL)){
  # --- Check gtx db connection
  gtx_debug("regional_context_analysis | validating gtx DB connection.")
  gtxdbcheck(dbc)
  
  
  
}

int_nearest_gene_gwas <- function(analysis, calc_th = FALSE, gwas_args){
  
}

int_nearest_gene_gene <- function(){
   
}

int_nearest_gene_pos <- function(){
  
}