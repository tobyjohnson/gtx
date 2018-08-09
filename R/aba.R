#' Coloc all-by-all functions
#' 
#' \strong{aba.query_locus() - Query a locus for coloc all-by-all results}
#' Query all colocalization data in a region. Typical use 
#' is to query a gene ID (\code{ensemblid} or \code{hgncid}) 
#' and to pull all coloc data in the surrounding the queried 
#' gene (+/- \code{surround} bp). Alternatively, can specifiy 
#' the region coordinates using: \code{chrom}, \code{pos_start}, 
#' & \code{pos_end}. The rate limiting step is the \code{\link{aba.query_locus}}, 
#' after which \code{\link{aba.filter}} & \code{\link{aba.plot}} are very fast. 
#' This means \code{\link{aba.query_locus}} results should be saved to an object,
#' but the following steps can be piped together to quickly iterate 
#' filtering and ploting. \code{\link{aba.query_locus}} takes ~ 7 min.
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param hgncid HGNC gene ID
#' @param ensemblid Ensembl gene ID 
#' @param surround [Default = 1e6] Define region as \code{surround} +/- gene TSS. 
#' @param rsid SNP rsid
#' @param chrom chromosome - Used to define a specific region
#' @param pos_start start position - Used to define a specific region, overrides surround
#' @param pos_end end position - Used to define a specific region, overrides surround
#' @param sc [Default =  getOption("gtx.sc")]. Spark connection. 
#' @param flog_level [Default = getOption("gtx.flog_level", "WARN")] \code{\link{futile.logger}} \code{\link{flog.threshold}} [INFO, WARN, ERROR]
#' @return data.frame with all coloc results in the region
#' @examples 
#' \strong{Query aba colocs:}
#' \code{colocs <- aba.query_locus(hgncid = "HMGCR")}
#' \code{colocs <- aba.query_locus(chrom = 1, pos_start = 2e3, pos_end = 4e5, sc = sc)}
#' 
#' \strong{To establish a spark connection:}
#' \code{sc <- spark_connect(master     = "yarn-client",
#'                           spark_home = "/opt/cloudera/parcels/SPARK2/lib/spark2",
#'                           version    = "2.1")}
#' \strong{Set global gtx.sc option:}               
#' \code{sc <- spark_connect(master = "yarn-client", spark_home = "/opt/cloudera/parcels/SPARK2/lib/spark2", version = "2.1")
#'      options(gtx.sc = sc)
#'      colocs <- aba.query_locus(hgncid = "HMGCR")}
#' @export
#' @import dplyr
#' @importFrom tidyr crossing
#' @import stringr
#' @import sparklyr
#' @import purrr
#' @import futile.logger
#' @import glue 
aba.query_locus <- function(hgncid, ensemblid, surround = 1e6, 
                            chrom, pos_start, pos_end, rsid,
                            sc         = getOption("gtx.sc", NULL),
                            flog_level = getOption("gtx.flog_level", "WARN")){
  if(flog_level != "INFO" & flog_level != "WARN" & flog_level != "ERROR"){
    flog.error(glue("aba.query_locus - flog_level defined as: {flog_level}. Needs to be [INFO, WARN, or ERROR]."))
  }
  flog.threshold(flog_level)
  flog.info("aba.query_locus - validating input.")
  # Validate an input gene or region was defined
  if(missing(hgncid) & missing(ensemblid) & missing(rsid) & (missing(chrom) & missing(pos_start) & missing(pos_end))){
    flog.error("aba.query_locus - must specify either: hgncid, ensemblid, or [chrom, pos_start, & pos_end]. Skipping.")
    next;
  }
  # Check we have a spark connection
  flog.info("aba.query_locus - validating Spark connection.")
  abaQ_sc_created = FALSE
  if(is.null(sc) | !any(str_detect(class(sc), "spark_connection"))){ 
    flog.info("aba.query_locus - no spark connection was passed to aba.query_locus. Will try to establish a connection.")
    # Try to establish a spark connection 
    safe_spark_connect <- safely(spark_connect)
    safe_sc <- safe_spark_connect(master     = "yarn-client",
                                  spark_home = "/opt/cloudera/parcels/SPARK2/lib/spark2",
                                  version    = "2.1")
    # If there was no error, use the connection
    if(is.null(safe_sc$error)){ 
      sc <- safe_sc$result 
      # Record that we made a connection in aba.query_locus
      abaQ_sc_created = TRUE
      flog.info("aba.query_locus - Spark connection established.")
    }
    # Otherwise advise user to manually create and pass the spark connection
    else if (!is.null(safe_sc$error)){
      flog.fatal(glue('aba.query_locus - unable to make the spark connection. Try initiating the spark connection manually and passing the connection to through option sc. To build a spark connection run:\n
\tsc <- spark_connect(master = {double_quote("yarn-client")}, \n\t\tspark_home = {double_quote("/opt/cloudera/parcels/SPARK2/lib/spark2")}, \n\t\tversion = {double_quote("2.1")})\n
After establishing the spark connection (above), you can also set the gtx options to use this connection by default running:\n
\t\toptions(gtx.sc = sc)'))
      next;
    }
  }
  
  # quote surround
  surround <- enquo(surround)
  
  # Define tables we need to use
  flog.info("aba.query_locus - establishing connection to database tables.")
  colocs_tbl   <- tbl(sc, "ukbiobank.coloc_results")
  gwas_th_tbl  <- tbl(sc, "ukbiobank.gwas_results_top_hits")
  genes_tbl    <- tbl(sc, "ukbiobank.genes")
  analyses_tbl <- tbl(sc, "ukbiobank.analyses")
  sites_tbl    <- tbl(sc, "ukbiobank.sites_ukb_500kv3")
  
  # Collect gene information to define region
  flog.info("aba.query_locus - collecting information to define region to look for GWAS top hits.")
  if(!missing(hgncid)){
    hgncid <- enquo(hgncid)
    # Verify the hgncid is reasonable
    if(!str_detect(quo_name(hgncid), "\\w+")){
      flog.error(glue("aba.query_locus - this hgncid ({tmp}) does not have characters ~ wrong.", 
                      tmp = quo_name(hgncid)))
      next;
    }
    # Query gene table with hgncid
    gene_info <- 
      genes_tbl %>% 
      filter(hgncid == !! hgncid) %>% 
      collect() 
    
    # Verify output
    if(nrow(gene_info) == 0){ 
      flog.warn(glue("aba.query_locus - could not find this hgncid ({tmp}) in the gene table.",  
                     tmp = quo_name(hgncid)))
      next;
    }  
    if(nrow(gene_info) >  1){ 
      flog.warn(glue("aba.query_locus - found more than one gene for this hgncid ({tmp}) in gene table, unsure which to use.",  
                     tmp = quo_name(hgncid)))
      next;
    }  
  }
  else if (!missing(ensemblid)){
    # Verify the ensemblid is reasonable
    ensemblid <- enquo(ensemblid)
    if(!str_detect(quo_name(ensemblid), "ENSG\\d+")){
      flog.error(glue("aba.query_locus - this ensemblid ({tmp}) doesn't look right.", 
                      tmp = quo_name(ensemblid)))
      next;
    }
    # Query gene table with ensemblid
    gene_info <- 
      genes_tbl %>% 
      filter(ensemblid == !! ensemblid) %>% 
      collect()
    
    # Verify output
    if(nrow(gene_info) == 0){ 
      flog.warn(glue("aba.query_locus - could not find this ensemblid ({tmp}) in the gene table.",
                     tmp = quo_name(ensemblid)))
      next;
    }  
    if(nrow(gene_info) >  1){ 
      flog.warn(glue("aba.query_locus - found more than one ensemblid ({tmp}) in gene table, unsure which to use.",
                     tmp = quo_name(ensemblid)))
      next;
    } 
  }
  else if(!missing(rsid)){
    # Verify the rsid is reasonable
    rsid <- enquo(rsid)
    if(!str_detect(quo_name(rsid), "[rs]?\\d+")){
      flog.error(glue("aba.query_locus - this rsid ({tmp}) doesn't look right.",
                      tmp = quo_name(rsid)))
      next;
    }
    # Remove "rs" from start of rsid
    rsid_clean <- str_replace(quo_name(rsid), "^rs", "")
    
    # Query gene table with hgncid
    gene_info <- 
      sites_tbl %>% 
      filter(rs == rsid_clean) %>% 
      select(chrom, pos_start = "pos") %>% 
      collect()
    
    # Verify output
    if(nrow(gene_info) == 0){ 
      flog.warn(glue("aba.query_locus - could not find this rsid ({tmp}) in the sties table.",
                     tmp = quo_name(rsid)))
      next;
    }  
    if(nrow(gene_info) >  1){ 
      flog.warn(glue("aba.query_locus - found more than one rsid ({tmp}) in sites table, unsure which to use.",
                     tmp = quo_name(rsid)))
      next;
    }  
  }
  else if(!missing(chrom) & !missing(pos_start) & !missing(pos_end)){
    # Verify inputs, note: this may change in the future (e.g. M, alt chrom, etc)
    if(!str_detect(chrom, "^[\\dXY]+$")){
      flog.error(glue("aba.query_locus - this chrom ({chrom}) doesn't match 1-22, X, or Y."))
      next;
    }

    # Set gene info
    gene_info <- tibble(
      chrom     = as.integer(chrom),
      pos_start = as.integer(pos_start),
      pos_end   = as.integer(pos_end))
  }
  else {
    flog.error("aba.query_locus - unable to properly handle input.")
    next;
  }
  
  # Determine filter region depending on input
  if(!missing(hgncid) | !missing(ensemblid) | !missing(rsid)){
    gene_info <- 
      gene_info %>% 
      mutate(pos_start_filter = pos_start - !! surround) %>% 
      mutate(pos_end_filter   = pos_start + !! surround)
  }
  else {
    gene_info <- 
      gene_info %>% 
      mutate(pos_start_filter = pos_start) %>% 
      mutate(pos_end_filter   = pos_end)
  }

  # Using calculated regions to find all GWAS top hits in the region
  flog.info(glue("aba.query_locus - pulling all GWAS top hits in the region: {gene_info$chrom}:{gene_info$pos_start_filter}-{gene_info$pos_end_filter}."))
  gwas_th <- 
    filter(gwas_th_tbl,
           chrom     == gene_info$chrom & 
           pos_start >= gene_info$pos_start_filter & 
           pos_end   <= gene_info$pos_end_filter &
           !is.na(min_pval) &
           !is.na(pval_index)) %>% 
    collect()

  # Verify GWAS top hits
  if(nrow(gwas_th) == 0){
    flog.warn(glue("aba.query_locus - no GWAS top hits found in this region: {gene_info$chrom}:{gene_info$pos_start_filter}-{gene_info$pos_end_filter}."))
    next;
  }
  
  # Summarize min and max GWAS top hits
  gwas_th_sum <- 
    gwas_th %>% 
    distinct(chrom, pos_index) %>% 
    summarize(min_pos = min(pos_index), 
              max_pos = max(pos_index))
  
  # Pull in all genes that could overlap the GWAS top hits
  flog.info(glue("aba.query_locus - expanding region based on genes with cis-windows that can overlap GWAS top hits."))
  locus_genes <- 
    genes_tbl %>% 
    filter(chrom     == gene_info$chrom & 
           pos_start >= gwas_th_sum$min_pos - !! surround &
           pos_end   <= gwas_th_sum$max_pos + !! surround) %>% 
    collect()
  
  # Verify locus genes
  if(nrow(locus_genes) == 0){
    flog.warn("aba.query_locus - no genes found to overlap the GWAS top hits in this region.\n")
    next;
  }
  
  # Pull coloc results for all locus genes + gwas top hits traits
  flog.info("aba.query_locus - pulling all colocs (GWAS + gene pairs) in the region.")
  dat2pull_tbl <-
    sdf_copy_to(sc,
                crossing(entity    = locus_genes$ensemblid,
                         analysis2 = gwas_th$analysis))
  
  coloc_region <- 
    left_join(dat2pull_tbl,
              colocs_tbl,
              by   = c("analysis2", "entity")) %>%
    left_join(., 
              genes_tbl %>% select(ensemblid, gene_start = pos_start, genetype, hgncid), 
              by = c("entity" = "ensemblid")) %>% 
    left_join(.,
              analyses_tbl %>% select(analysis, description, ncase, ncontrol, ncohort),
              by = c("analysis2" = "analysis")) %>% 
    left_join(.,
              gwas_th %>% select(analysis, gwas_chrom = chrom, gwas_th_pos = pos_index, gwas_th_pval = pval_index),
              by   = c("analysis2" = "analysis"),
              copy = TRUE) %>% 
    left_join(.,
              sites_tbl %>% select(gwas_th_pos = pos, gwas_chrom = chrom, rs),
              by = c("gwas_chrom", "gwas_th_pos")) %>%
    collect() %>% 
    mutate(log10_min_pval_colocGWAS_over_TopHitGWAS = log10(gwas_th_pval / minpval2)) %>% 
    mutate(rsid = paste0("rs", rs)) %>% 
    mutate(minpval1_sci = formatC(minpval1, format = "e", digits = 2)) %>% 
    mutate(minpval2_sci = formatC(minpval2, format = "e", digits = 2))
  
  # Validate we got data back
  if(nrow(coloc_region) == 0){
    flog.warn("aba.query_locus - no colocs found in this region.")
    next;
  }
  
  # If we created a spark connection inside aba.query_locus, disconnect it
  if( isTRUE(abaQ_sc_created) ){ spark_disconnect(sc) }
  
  # Return colocs found in the region
  flog.info("aba.query_locus - COMPLETED")
  return(coloc_region) 
}

#' \strong{aba.filter() - Recommended filtering for aba results}
#' 
#' The aba results and region queries will return colocs
#' that are not robust and therefore need further filtering.
#' This function provdies recommended minimal filtering to 
#' improve the robustness of the aba results. In addition, 
#' the results will fill in the matrix of positive colocs
#' based on the filtering to enable complete plotting. The
#' ability to fill in the matrix can be toggled.
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param .data \code{\link{aba.query_locus}} results object to filter
#' @param fill_matrix [Default = TRUE] Fill in the matrix based filtering. \emph{required} for \code{\link{aba.plot}}
#' @param p12_ge [Default >= 0.80] This is the "H4" posterior probability cutoff 
#' @param minpval1_lt [Default <= 1e-5] Min pval seen in the eQTL                
#' @param minpval2_lt [Default <= 5e-8] Min pval seen in the GWAS data           
#' @param log10_min_pval_colocGWAS_over_TopHitGWAS_lt [Default < 0.5] Log10 ratio of coloc-gwas/top-hit-gwas pval. 
#' This ensures that the coloc gwas and gwas top hits are 1/2 order of magnitude apart and therefore shared.
#' @param gwas2gene_dist_lt [Default = 1e6] Max distance between coloc gene and gwas top hit. 
#' @param ncase_lt [Default >= 200] Minimum ncases for traits                   
#' @param ncohort_lt [Default >= 200] Minimum ncohort for traits
#' @param protein_coding [Default = TRUE] Filter only for protein coding transcripts
#' @param neale_only [Default = FALSE] Filter onyl for Neale traits, reduces redundancy b/w GSK & Neale data.
#' @param gsk_only [Default = FALSE] Filter only for GSK traits, reduces redundancy b/w GSK & Neale data.
#' @param flog_level [Default = getOption("gtx.flog_level", "WARN")] \code{\link{futile.logger}} \code{\link{flog.threshold}} [INFO, WARN, ERROR]
#' @return data.frame with input \code{\link{aba.query_locus}} reasonably filtered
#' @examples 
#' \strong{Basic use:}
#' \code{colocs_filtered <- aba.filter(colocs)}
#' \code{colocs_filtered <- aba.filter(colocs, p12 = 0.9, minpval1 = 5e-8)}
#' 
#' \strong{String together}
#' \code{colocs_filtered <- aba.query_locus(hgncid = "HMGCR") %>% aba.filter()}
#' @export
#' @import dplyr
#' @importFrom tidyr crossing
#' @import futile.logger
#' @import glue
aba.filter <- function(.data, p12_ge = 0.80, 
                       minpval1_le   = 1e-5, minpval2_le = 5e-8,
                       log10_min_pval_colocGWAS_over_TopHitGWAS_lt = 0.5, 
                       gwas2gene_dist_lt = 1e6,
                       ncase_ge = 200, ncohort_ge = 200, 
                       protein_coding_only = TRUE, 
                       neale_only  = FALSE, 
                       gsk_only    = FALSE, 
                       fill_matrix = TRUE, 
                       flog_level  = getOption("gtx.flog_level", "WARN")){
  if(flog_level != "INFO" & flog_level != "WARN" & flog_level != "ERROR"){
    flog.error(glue("aba.query_locus - flog_level defined as: {flog_level}. Needs to be [INFO, WARN, or ERROR]."))
  }
  flog.threshold(flog_level)
  # Verify input
  flog.info("aba.filter - validating input")
  input <- .data
  required_cols <- c("p12", "minpval1", "minpval2", "log10_min_pval_colocGWAS_over_TopHitGWAS", 
                     "gwas_th_pos", "gene_start", "ncase", "ncohort", "genetype", "analysis2")
  if(!all(required_cols %in% (names(input)))){
    flog.error(paste0("aba.filter - input is missing required cols. Required cols include:", paste(required_cols, collapse = ", ")))
    next;
  }
  if(!is.numeric(p12_ge) | !is.numeric(gwas2gene_dist_lt) | 
     !is.numeric(minpval1_le) | !is.numeric(minpval2_le) |
     !is.numeric(log10_min_pval_colocGWAS_over_TopHitGWAS_lt) | 
     !is.numeric(ncase_ge) | !is.numeric(ncohort_ge)){
    flog.error("aba.filter - input parameter is not numeric.")
    next;
  }
  
  # Perform basic uniform input
  flog.info("aba.filter - filtering data")
  input_filtered <-
    input %>% 
    filter(p12      >= !! p12_ge & # filter H4 >= 0.80
           minpval1 <= !! minpval1_le  & # min eQTL pval
           minpval2 <= !! minpval2_le  & # min gwas pval
           abs(log10_min_pval_colocGWAS_over_TopHitGWAS) < !! log10_min_pval_colocGWAS_over_TopHitGWAS_lt & # coloc min gwas pval must be 1/2 order of magnitude from gwas top hits
           abs(gwas_th_pos - gene_start)  < !! gwas2gene_dist_lt & # & # make sure the gene and GWAS hit are < 1 MB apart
           (ncase   >= !! ncase_ge | ncohort >= !! ncohort_ge))
  
  # Filter for protein coding transcripts
  if(isTRUE(protein_coding_only)){
    input_filtered <- 
      input_filtered %>% 
      filter(genetype == "protein_coding")
  }

  # Or filter for GSK data only
  if(isTRUE(gsk_only)){
    input_filtered <- 
      input_filtered %>% 
      filter(str_detect(analysis2, "GSK"))
  }
  # Filter for Neale data only
  else if(isTRUE(neale_only)){
    input_filtered <- 
      input_filtered %>% 
      filter(str_detect(analysis2, "neale"))
  }
  
  if(nrow(input_filtered) == 0){
    flog.warn("After filtering, no results are remaining.")
    next;
  }
  
  # Fill in the matrix based on the filtered results
  flog.info("aba.filter - filling in matrix of filtered results")
  if(isTRUE(fill_matrix)){
    input_filtered <- 
      crossing(entity    = input_filtered$entity,
               analysis2 = input_filtered$analysis2,
               analysis1 = input_filtered$analysis1) %>% 
      inner_join(., 
                 input,
                 by = c("entity", "analysis2", "analysis1")) %>% 
      collect()
  }
  
  flog.info("aba.filter - COMPLETED")
  return(input_filtered)
}

#' \strong{aba.plot() - Plot aba results}
#' 
#' This function tries to create a bubble plot to 
#' illustrate the coloc region. This function returns
#' a ggplot object which can be viewed or saved. It 
#' is \emph{highly} recommended to save the plot to pdf 
#' for best visualization of the plot. 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param p12_ge [Default = 0.80] Values below this will be colored as non-significant, while above have color for direction. 
#' @param max_dot_size [Default = 5] Adjust max size of dots.
#' @param title Set the plot title. NULL will remove title. 
#' @param flog_level [Default = getOption("gtx.flog_level", "WARN")] \code{\link{futile.logger}} \code{\link{flog.threshold}} [INFO, WARN, ERROR]
#' @return ggplot2 object for viz and export.
#' @examples 
#' \strong{Basic use:}
#' coloc_plot <- \code{\link{aba.plot}}(colocs_filtered)
#' coloc_plot <- \code{\link{aba.plot}}(colocs_filtered, max_dot_size = 3, title = "Cool plot")
#' 
#' \strong{View the returned object:}
#' \code{coloc_plot}
#' 
#' \strong{Save the object to PDF:}
#' \code{ggsave(filename = "hgmcr_region_coloc.pdf", 
#'        plot     = coloc_plot, 
#'        width    = 11, 
#'        height   = 8.5, 
#'        dpi      = 300)}
#'        
#' \strong{Advanced use: query data, remove death related traits, filter, and then plot}
#' \code{coloc_plot <- aba.query_locus(hgncid = "HMGCR", sc = sc) %>% 
#'                 filter(!str_detect(description, "death")) %>% 
#'                 aba.filter() %>% 
#'                 aba.plot(title = "Best plot ever")}
#' @export
#' @import dplyr
#' @import stringr
#' @import ggplot2
#' @import futile.logger
#' @import glue
aba.plot <- function(.data, p12_ge = 0.80, max_dot_size = 5, title = NULL,
                     flog_level = getOption("gtx.flog_level", "WARN")){
  if(flog_level != "INFO" & flog_level != "WARN" & flog_level != "ERROR"){
    flog.error(glue("aba.query_locus - flog_level defined as: {flog_level}. Needs to be [INFO, WARN, or ERROR]."))
  }
  flog.threshold(flog_level)
  # Verify input
  flog.info("aba.plot - validating input")
  input <- .data
  required_cols <- c("analysis1", "analysis2", "description", "hgncid", "p12", "alpha21")
  if(!all(required_cols %in% (names(input)))){
    flog.error(paste0("aba.plot - input is missing required cols. Required cols include:", paste(required_cols, collapse = ", ")))
    next;
  }
  if(!is.numeric(max_dot_size)){
    flog.warn("aba.plot - max_dot_size parameter is not numeric.")
  }
  if(!is.null(title) & !is.character(title)){
    flog.warn("aba.plot - title is neither NULL or character string.")
  }
  
  # Clean data (e.g. tissue & description), add simple direction of effect
  flog.info("aba.plot - cleaing data for plotting")
  fig_dat <- 
    input %>% 
    # clean tissue names
    mutate(analysis1 = str_replace_all(analysis1, "_", " ")) %>%
    mutate(analysis1 = str_replace_all(analysis1, "gtex7", "")) %>%
    mutate(analysis1 = str_to_title(analysis1)) %>% 
    mutate(analysis1 = str_trim(analysis1)) %>% 
    # clean GWAS description names
    mutate(description = str_replace_all(description, "_", " ")) %>% 
    # Fill in missing hgncdids with ensemblids
    mutate(hgncid = case_when(!str_detect(hgncid, "\\w+") ~ entity,
                              is.na(hgncid)               ~ entity,
                              TRUE                        ~ hgncid)) %>% 
    # Remove Broad in description b/c only using broad data here
    mutate(description = str_replace(description, "\\(UKB Broad\\)", "")) %>% 
    mutate(description_wrap = str_wrap(description,  width = 60, exdent = 8)) %>%
    # Adjust alpha21 so that non-significant (H4 < 0.80) = NA 
    mutate(direction = case_when(p12 >= !! p12_ge & alpha21 >  0        ~ "Increased", 
                                 p12 >= !! p12_ge & alpha21 <  0        ~ "Decreased", 
                                 p12 >= !! p12_ge & alpha21 == 0        ~ "None",
                                 p12 <  !! p12_ge & is.numeric(alpha21) ~ "Not Sig.")) %>% 
    mutate(rs_pos = paste0(rsid, "\n", "chr", gwas_chrom, ":", gwas_th_pos ))
  
  flog.info("aba.plot - ordering results by chromosome - position coordinates")
  # Order genes (cols) and rows (rs ids) by position
  fig_dat <- within(fig_dat, {
    hgncid <- reorder(hgncid, gene_start)
    rs_pos <- reorder(rs_pos, gwas_th_pos) 
  })
  
  flog.info("aba.plot - plotting data")
  fig <-
    fig_dat %>% 
    ggplot(aes(x = analysis1, y = description_wrap)) + 
    geom_point(aes(size = p12, fill = direction), shape = 21, color = "black", stroke = 0.2) +
    facet_grid(rs_pos ~ hgncid, scales = "free", space = "free", drop = TRUE, margins = FALSE) +
    scale_fill_manual(values = c("Decreased" = "blue", 
                                 "None"      = "white",
                                 "Increased" = "red",
                                 "Not Sig."  = "grey")) +
    scale_radius(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.00),
                 range  = c(1, max_dot_size)) + 
    guides(fill  = guide_legend(title = "Change in gene expression\nwith increased GWAS trait")) + 
    theme_bw() + 
    theme(legend.position  = "bottom",
          legend.text      = element_text(size = 8),
          legend.title     = element_text(face = "bold"),
          axis.text.x      = element_text(size = 6, angle = 35, hjust = 1),
          axis.text.y      = element_text(size = 6, angle = 0 ),
          axis.line        = element_line(color = "black"),
          axis.title.x     = element_blank(),
          axis.title.y     = element_blank(),
          panel.grid.major = element_line(size = 0.25),
          panel.grid.minor = element_line(size = 0.25),
          strip.text.x     = element_text(face = "bold"),
          strip.text.y     = element_text(size = 5, angle = 0),
          plot.title       = element_text(hjust = 0.5)) +
    guides(size  = guide_legend(title = "Posterior probability\nof colocalization"))
  
  if(!is.null(title)){
    fig <- fig + ggtitle(title)
  }
  
  flog.info("aba.plot - COMPLETED")
  return(fig)
}

#' \strong{aba.query_locus_wrapper() - single function to wrap around: [\code{\link{aba.query_locus}}, \code{\link{aba.filter}}, \code{\link{aba.plot}}]}
#' 
#' This function is a simple wrapper around \code{\link{aba.query_locus}}, 
#' \code{\link{aba.filter}}, & \code{\link{aba.plot}}. This function 
#' is ideal for high throughput generation of many plots. This 
#' fucntion is \emph{not} ideal for deep dives into a locus.
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param hgncid HGNC gene ID
#' @param ensemblid Ensembl gene ID 
#' @param rsid SNP rsid
#' @param surround [Default = 1e6] Define region as \code{surround} +/- input gene TSS. 
#' @param chrom chromosome - Used to define a specific region
#' @param pos_start start position - Used to define a specific region, overrides surround
#' @param pos_end end position - Used to define a specific region, overrides surround
#' @param sc Spark connection
#' @param ... Pass filtering options to \code{\link{aba.filter}}
#' @param flog_level [Default = getOption("gtx.flog_level", "WARN")] \code{\link{futile.logger}} \code{\link{flog.threshold}} [INFO, WARN, ERROR]
#' @return ggplot2 object for viz and export.
#' @examples  
#' \strong{Basic use:}
#' \code{hmgcr_fig <- aba.wrapper(hgncid = "HMGCR")}
#' \strong{Advanced use:}
#' map( wrapper)
#' map( export)
#' @export
#' @import dplyr
#' @import stringr
#' @import sparklyr
#' @import ggplot2
#' @import glue
#' @import purrr
#' @import futile.logger
aba.query_locus_wrapper <- function(hgncid, ensemblid, surround = 1e6, 
                        chrom, pos_start, pos_end, rsid,
                        sc = getOption("gtx.sc", NULL), 
                        flog_level = getOption("gtx.flog_level", "WARN"), 
                        ...){
  if(flog_level != "INFO" & flog_level != "WARN" & flog_level != "ERROR"){
    flog.error(glue("aba.query_locus - flog_level defined as: {flog_level}. Needs to be [INFO, WARN, or ERROR]."))
  }
  flog.threshold(flog_level)
  flog.info("aba.wrapper - validating Spark connection")
  # Check we have a spark connection
  abaW_sc_created = FALSE
  if(is.null(sc) | !any(str_detect(class(sc), "spark_connection"))){ 
    flog.info("aba.wrapper - did not find Spark connection, trying to establish a connection now.")
    # Try to establish a spark connection 
    safe_spark_connect <- safely(spark_connect)
    safe_sc <- safe_spark_connect(master     = "yarn-client",
                                  spark_home = "/opt/cloudera/parcels/SPARK2/lib/spark2",
                                  version    = "2.1")
    # If there was no error, use the connection
    if(is.null(safe_sc$error)){ 
      sc <- safe_sc$result
      # Record that we made a connection in aba.query_locus
      abaW_sc_created = TRUE
      flog.info("aba.wrapper - Spark connection established.")
    }
    # Otherwise advise user to manually create and pass the spark connection
    else if (!is.null(safe_sc$error)){
      stop(glue('aba.query_locus - unable to make the spark connection. Try initiating the spark connection manually and passing the connection to through option sc. To build a spark connection run:\n
\tsc <- spark_connect(master = {double_quote("yarn-client")}, \n\t\tspark_home = {double_quote("/opt/cloudera/parcels/SPARK2/lib/spark2")}, \n\t\tversion = {double_quote("2.1")})\n
                After establishing the spark connection (above), you can also set the gtx options to use this connection by default running:\n
                \t\toptions(gtx.sc = sc)'),
           call. = FALSE)
    }
  }
  
  flog.info("aba.wrapper - validating input & making query of colocs in region")
  surround <- enquo(surround)
  # Determine which input to use and make query
  if(!missing(hgncid)){
    hgncid <- enquo(hgncid)
    # Verify the hgncid is reasonable
    if(!str_detect(quo_name(hgncid), "\\w+")){
      flog.warn(glue("aba.wrapper - this hgncid ({tmp}) doesn't look right.",  
                     tmp = quo_name(hgncid)))
    }
    # Query region using hgncid
    colocs <- aba.query_locus(hgncid = !! hgncid, surround = !! surround, sc = sc)
    wrap_title <- glue("Colocalizations around {tmp} (+/- {tmp2})", 
                       tmp  = quo_name(hgncid),
                       tmp2 = quo_name(surround))
  }
  else if(!missing(ensemblid)){
    ensemblid <- enquo(ensemblid)
    # Verify the ensemblid is reasonable
    if(!str_detect(quo_name(ensemblid), "ENSG\\d+")){
      flog.warn("aba.wrapper - this ensemblid doesn't look right.")
      next;
    }
    # Query region using ensemblid
    colocs <- aba.query_locus(ensemblid = !! ensemblid, surround = !! surround, sc = sc)
    wrap_title <- glue("Colocalizations around {tmp} (+/- {tmp2})", 
                       tmp  = quo_name(ensemblid), 
                       tmp2 = quo_name(surround))
  }
  else if(!missing(rsid)){
    rsid <- enquo(rsid)
    # Verify the rsid is reasonable
    if(!str_detect(quo_name(rsid), "[rs]?\\d+")){
      flog.warn(glue("aba.wrapper - this rsid ({tmp}) doesn't look right.", 
                     tmp = quo_name(rsid)))
    }
    # Query region using rsid
    colocs <- aba.query_locus(rsid = !! rsid, surround = !! surround, sc = sc)
    wrap_title <- glue("Colocalizations around {tmp} (+/- {tmp2})", 
                       tmp  = quo_name(rsid), 
                       tmp2 = quo_name(surround))
  }
  else if(!missing(chrom) & !missing(pos_start) & !missing(pos_end)){
    # Verify inputs, note: this may change in the future (e.g. M, alt chrom, etc)
    if(!str_detect(chrom, "^[\\dXY]+$")){
      flog.warn("aba.wrapper - this chrom doesn't match 1-22, X, or Y.")
      next;
    }
    # Query region using hgncid
    colocs    <- aba.query_locus(chrom     =  chrom, 
                                 pos_start =  pos_start, pos_end = pos_end, 
                                 surround  =  surround,  sc = sc)
    wrap_title <- glue("Colocalizations between {chrom}:{pos_start}-{pos_end}")
  }
  # Run all the aba functions
  fig <- 
    colocs %>% 
    aba.filter(...) %>% 
    aba.plot(title = wrap_title)

  # If we created a spark connection inside aba.wrapper, disconnect it
  if( isTRUE(abaW_sc_created) ){ spark_disconnect(sc) }
  
  flog.info("aba.wrapper - COMPLETED")
  return(fig)
}

#' \strong{aba.query_traits() - Query gwas traits for coloc all-by-all results}
#' Query a list of GWAS traits for eQTL coloc results.
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param analysis_ids ukbiobank analysis id(s), single string or vector. 
#' @param flatten_data [Default = TRUE] Condense coloc data into single gene per row. 
#' @param surround [Default = 1e6] Distance from GWAS top hit to TSS. 
#' @param protein_coding_only [Default = TRUE] Filtering to only return protein coding transcripts. 
#' @param progressive_validation [Default = FALSE] TRUE = iteratively check data is obtained from the multiple queries. This dramatically slows down analysis, but is useful for trouble shooting.
#' @param aba_filter [Default = TRUE] Apply \code{\link{aba.filter}} to help validate colocs.
#' @param ... Options to pass to \code{\link{aba.filter}}
#' @param sc [Default =  getOption("gtx.sc")]. Spark connection. 
#' @param flog_level [Default = getOption("gtx.flog_level", "WARN")] \code{\link{futile.logger}} \code{\link{flog.threshold}} [INFO, WARN, ERROR]
#' @return a list with 2 data.frames. "flat" = flatten_data compact, 1 gene / row. "full" = full coloc results (after filtering)
#' @examples 
#' traits <- c("aid1", "aid2", "aid3")
#' traits_colocs <- aba.query_traits(analysis_ids = traits)
#' traits_colocs$flat %>% filter(hgncid == "Gene1")
#' traits_colocs$full %>% filter(hgncid == "Gene1")
#' @export
#' @import dplyr
#' @importFrom tidyr crossing
#' @import stringr
#' @import sparklyr
#' @import purrr
#' @import futile.logger
#' @import glue 
aba.query_traits <- function(analysis_ids,
                             aba_filter = TRUE,
                             flatten_data = TRUE,
                             surround = 1e6,
                             protein_coding_only = TRUE,
                             progressive_validation = FALSE,
                             sc         = getOption("gtx.sc", NULL),
                             flog_level = getOption("gtx.flog_level", "WARN"),
                             ...){
  if(flog_level != "INFO" & flog_level != "WARN" & flog_level != "ERROR"){
    flog.error(glue("aba.query_traits - flog_level defined as: {flog_level}. Needs to be [INFO, WARN, or ERROR]."))
  }
  flog.threshold(flog_level)
  flog.info("aba.query_traits - validating input.")
  # Validate an input gene or region was defined
  if(missing(analysis_ids)){
    flog.error("aba.query_traits - must specify either: analysis_ids")
    next;
  }
  # Check we have a spark connection
  flog.info("aba.query_traits - validating Spark connection.")
  abaQ_sc_created = FALSE
  if(is.null(sc) | !any(str_detect(class(sc), "spark_connection"))){ 
    flog.info("aba.query_traits - no spark connection was passed to aba.query_locus. Will try to establish a connection.")
    # Try to establish a spark connection 
    safe_spark_connect <- safely(spark_connect)
    safe_sc <- safe_spark_connect(master     = "yarn-client",
                                  spark_home = "/opt/cloudera/parcels/SPARK2/lib/spark2",
                                  version    = "2.1")
    # If there was no error, use the connection
    if(is.null(safe_sc$error)){ 
      sc <- safe_sc$result 
      # Record that we made a connection in aba.query_locus
      abaQ_sc_created = TRUE
      flog.info("aba.query_traits - Spark connection established.")
    }
    # Otherwise advise user to manually create and pass the spark connection
    else if (!is.null(safe_sc$error)){
      flog.fatal(glue('aba.query_traits - unable to make the spark connection. Try initiating the spark connection manually and passing the connection to through option sc. To build a spark connection run:\n
                      \tsc <- spark_connect(master = {double_quote("yarn-client")}, \n\t\tspark_home = {double_quote("/opt/cloudera/parcels/SPARK2/lib/spark2")}, \n\t\tversion = {double_quote("2.1")})\n
                      After establishing the spark connection (above), you can also set the gtx options to use this connection by default running:\n
                      \t\toptions(gtx.sc = sc)'))
      next;
    }
  }
  flog.info("aba.query_traits - establishing connection to database tables.")
  analyses_tbl <- tbl(sc, "ukbiobank.analyses")
  gwas_th_tbl  <- tbl(sc, "ukbiobank.gwas_results_top_hits")
  genes_tbl    <- tbl(sc, "ukbiobank.genes")
  colocs_tbl   <- tbl(sc, "ukbiobank.coloc_results")
  sites_tbl    <- tbl(sc, "ukbiobank.sites_ukb_500kv3")
  
  # Make a tmp table with our input traits
  traits     <- crossing(analysis_ids)
  traits_tbl <- sdf_copy_to(sc, traits, overwrite = TRUE)
  
  # Make sure we can find these analysis ids
  flog.info("aba.query_traits - pulling GWAS info about analysis ids.")
  traits_tbl <- 
    inner_join(
      traits_tbl,
      analyses_tbl %>% select(analysis, description, ncase, ncohort),
      by = c("analysis_ids" = "analysis"))
  
  # Validate we got GWAS top hits back
  if(isTRUE(progressive_validation)){
    if(traits_tbl %>% collect() %>% nrow() == 0){
      flog.warn("aba.query_traits - no GWAS analysis ids found to match input.")
      next;
    }
  }
  
  # Identify GWAS top hits in the input traits
  flog.info("aba.query_traits - pulling GWAS top hits.")
  traits_th_tbl <- 
    inner_join(
      traits_tbl,
      gwas_th_tbl,
      by = c("analysis_ids" = "analysis")) %>% 
    select(analysis2 = analysis_ids, description, ncase, ncohort, gwas_th_chrom = chrom, gwas_num_var = num_variants,
           gwas_th_pos = pos_index, gwas_th_pval = pval_index) %>% 
    left_join(.,
              sites_tbl %>% select(pos, chrom, rs),
              by = c("gwas_th_pos" = "pos", "gwas_th_chrom" = "chrom"))

  # Validate we got GWAS top hits back
  if(isTRUE(progressive_validation)){
    if(traits_th_tbl %>% collect() %>% nrow() == 0){
      flog.warn("aba.query_traits - no GWAS top hits found.")
      next;
    }
  }
  
  # Find all the genes that overlap each GWAS top hit
  flog.info("aba.query_traits - pulling genes that overlap GWAS top hits")
  traits_th_genes_tbl <- 
    left_join(
      traits_th_tbl,
      genes_tbl,
      by = c("gwas_th_chrom" = "chrom")) %>% 
    # Make sure the GWAS top hit is +/- {surround} from gene TSS (cis-window for GTEx)
    filter(gwas_th_pos >= pos_start - surround & 
           gwas_th_pos <= pos_start + surround)
  
  if(isTRUE(protein_coding_only)){
    flog.info("aba.query_traits - selecting only protein coding transcripts")
    traits_th_genes_tbl <- 
      traits_th_genes_tbl %>% 
      filter(genetype == "protein_coding")
  }
  
  # Validate we got genes overlapping GWAS top hits
  if(isTRUE(progressive_validation)){
    if(traits_th_genes_tbl %>% collect() %>% nrow() == 0){
      flog.error("aba.query_traits - did not find any genes that overlap the GWAS top hits.")
      next;
    }
  }
  
  # Pull coloc data 
  flog.info("aba.query_traits - pulling colocs for input GWAS & overlapping genes.")
  traits_colocs <-
    left_join(traits_th_genes_tbl,
              colocs_tbl,
              by = c("analysis2", "ensemblid" = "entity")) %>% 
    collect() %>% 
    mutate(log10_min_pval_colocGWAS_over_TopHitGWAS = log10(gwas_th_pval / minpval2)) %>% 
    mutate(rsid = paste0("rs", rs)) %>% 
    rename(gene_start = "pos_start", gene_end = "pos_end")
  
  # Validate we got genes overlapping GWAS top hits
  if(traits_colocs %>% nrow() == 0){
    flog.warn("aba.query_traits - did not find any colocs for these traits.")
    next;
  }
  
  # Filter for significant results
  if(isTRUE(aba_filter)){
    flog.info("aba.query_traits - performing aba.filter() on the coloc data.")
    traits_colocs <-
      traits_colocs %>% 
      aba.filter(fill_matrix = FALSE, ...)
  }

  # Validate we got results back after filtering
  if(traits_colocs %>% nrow() == 0){
    flog.warn("aba.query_traits - after aba.filter no colocs remain.")
    next;
  }
  
  # Filter for significant results
  if(isTRUE(flatten_data)){
    flog.info("aba.query_traits - flattening data")
    flat_dat <-
      traits_colocs %>% 
      group_by(rsid) %>% 
      mutate(locus_n_genes   = n_distinct(ensemblid)) %>% 
      mutate(locus_n_tissue  = n_distinct(analysis1)) %>% 
      group_by(ensemblid) %>%
      mutate(gene_n_top_hits = n_distinct(rsid)) %>%
      mutate(gene_n_traits   = n_distinct(analysis2)) %>% 
      mutate(gene_n_tissue   = n_distinct(analysis1)) %>% 
      mutate(loci_max_n_genes  = case_when(gene_n_top_hits >  1 ~ max(locus_n_genes) %>% as.integer(), 
                                           gene_n_top_hits <= 1 ~ NA_integer_)) %>% 
      mutate(loci_min_n_genes  = case_when(gene_n_top_hits >  1 ~ min(locus_n_genes %>% as.integer()), 
                                           gene_n_top_hits <= 1 ~ NA_integer_)) %>% 
      mutate(loci_max_n_tissue = case_when(gene_n_top_hits >  1 ~ max(locus_n_tissue %>% as.integer()), 
                                           gene_n_top_hits <= 1 ~ NA_integer_)) %>% 
      mutate(loci_min_n_tissue = case_when(gene_n_top_hits >  1 ~ min(locus_n_tissue %>% as.integer()), 
                                           gene_n_top_hits <= 1 ~ NA_integer_)) %>% 
      ungroup() %>% 
      mutate(locus_n_genes     = case_when(gene_n_top_hits >  1 ~ NA_integer_,
                                           gene_n_top_hits <= 1 ~ locus_n_genes %>% as.integer())) %>% 
      mutate(locus_n_tissue    = case_when(gene_n_top_hits >  1 ~ NA_integer_,
                                           gene_n_top_hits <= 1 ~ locus_n_tissue %>% as.integer())) %>% 
      select(ensemblid, hgncid, 
             gene_n_traits, gene_n_tissue, gene_n_top_hits,
             locus_n_genes,     locus_n_tissue, 
             loci_max_n_genes,  loci_min_n_genes,
             loci_max_n_tissue, loci_min_n_tissue) %>% 
      distinct() %>% 
      arrange(gene_n_top_hits, gene_n_tissue, locus_n_genes) 
  }
  
  # Validate we got results back after filtering
  if(traits_colocs %>% nrow() == 0){
    flog.warn("aba.query_traits - after flattening, no colocs remain. Likely error/bug.")
    next;
  }
  
  flog.info("aba.query_traits - COMPLETED")
  if(isTRUE(flatten_data)){
    return(list(flat = flat_dat, full = traits_colocs))
  }
  else {
    return(list(full = traits_colocs))
  }
} 

#' \strong{aba.info_table() - description of the cols from \code{\link{aba.query_locus}}}
#' 
#' Running this function will return a data frame with 
#' descriptions for each column returned by \code{\link{aba.query_locus}}}. 
#' These column names are consistent with other RDIP tables for joins. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @return data.frame
#' @examples
#' \code{info <- aba.info_table()}
#' @export
#' @import dplyr
aba.info_table <- function(){
  ret <- tibble(
    col_name = c("analysis1", 
                 "entity", 
                 "p0", 
                 "p1", 
                 "p2", 
                 "p1c2", 
                 "p12",
                 "alpha21", 
                 "minpval1", 
                 "minpval2", 
                 "analysis2", 
                 "gene_start",                  
                 "genetype", 
                 "hgncid", 
                 "description", 
                 "ncase", 
                 "ncontrol",
                 "ncohort", 
                 "gwas_chrom", 
                 "gwas_pos_index", 
                 "gwas_th_pval", 
                 "rs", 
                 "log10_min_pval_colocGWAS_over_TopHitGWAS", 
                 "minpval1_sci", 
                 "minpval2_sci"),
     description = c("eQTL tissue name", 
                     "transcript ensembld ID", 
                     "coloc H0 - No association", 
                     "coloc H1 - One variant associated with phenotype 1 only", 
                     "coloc H2 - One variant associated with phenotype 2 only", 
                     "coloc H3 - Two variants separately associated with phenotypes 1 and 2", 
                     "coloc H4 - One variant associated with phenotypes 1 and 2",
                     "direction of effect - average betas weighted posterior probability of the GWAS over the eQTL", 
                     "min eQTL pval seen in coloc", 
                     "min GWAS pval seen in coloc", 
                     "GWAS analysis ID", 
                     "entity start position", 
                     "type of entity e.g., protein coding, miRNA, etc", 
                     "Gene symbol", 
                     "description of GWAS trait", 
                     "number of cases in GWAS trait", 
                     "number of controls in GWAS trait",
                     "number of cohort in GWAS trait", 
                     "GWAS top hit chromosome", 
                     "GWAS top hit index position", 
                     "GWAS top hit index pval",
                     "GWAS top hit RSID",
                     "log10 of the minpval2 / gwas_th_pval - Used to ensure the coloc GWAS signal overlaps primary top hit signal",
                     "minpval1 formated as scientific string", 
                     "minpval2 formated as scientific string"))
  
  return(ret)
}


#' \strong{aba.info_opt() - info about initiating and setting global opts for aba}
#' 
#' Running this function will return a data frame with 
#' descriptions for each column returned by \code{\link{aba.query_locus}}}. 
#' These column names are consistent with other RDIP tables for joins. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @return data.frame
#' @examples
#' \code{info <- aba.info_opt()}
#' @export
#' @import dplyr
aba.info_opt <- function(){
  ret <- tibble(
         opt = c("sc",
                 "flog_level"),
         info = c("Spark connection", 
                  "futile logger threshold"),
         init = c(glue('sc <- spark_connect(master = {double_quote("yarn-client")}, spark_home = {double_quote("/opt/cloudera/parcels/SPARK2/lib/spark2")}, version = {double_quote("2.1")})'),
                  glue('options(gtx.flog_level = {double_quote("WARN")})')),
         global_opt = c(glue('options(gtx.sc = sc)'),
                        glue('options(gtx.flog_level = {double_quote("WARN")})')))
  
  return(ret)
}