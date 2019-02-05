#' aba.query() - Query aba results for all positive results
#' 
#' Query all colocalization data in a region or set of traits. Typical use 
#' is to query a gene ID (\code{ensemblid} or \code{hgncid}) 
#' and to pull all coloc data in the surrounding queried 
#' region (+/- \code{surround} bp). A region can also be specified
#' using: \code{chrom}, \code{pos_start}, & \code{pos_end}. Alternatively,
#' a set of traits can be queried instead of the entire aba to help
#' generate hypothesis before looking at specific loci. The 
#' \code{\link{aba.query()}} will return only positive results based on
#' filtering (using params below). 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param analysis_ids ukbiobank analysis id(s), single string or vector.
#' @param hgncid HGNC symbol. Single string or vector. 
#' @param ensemblid Ensembl gene ID. Single string or vector.
#' @param rsid SNP rsid. Single string or vector.
#' @param surround [Default = 1e6] Distance to pull in TH + respective genes. 
#' @param chrom chromosome - Used to define a specific region
#' @param pos position that will have +/- \code{surround} for bounds to pull in colocs.
#' @param pos_start start position - Used to define a specific region, overrides surround
#' @param pos_end end position - Used to define a specific region, overrides surround
#' @param p12_ge [Default >= 0.80] This is the "H4" posterior probability cutoff 
#' @param minpval1_le [Default <= 1e-4] Min pval seen in the eQTL                
#' @param minpval2_le [Default <= 5e-8] Min pval seen in the GWAS data           
#' @param ncase_ge [Default >= 200] Minimum ncases for traits.                   
#' @param ncohort_ge [Default >= 200] Minimum ncohort for traits.
#' @param protein_coding_only [Default = FALSE] Filter only for protein coding transcripts
#' @param ttam_only [Default = FALSE] Filter only for 23andMe (TTAM) traits.
#' @param gsk_only [Default = FALSE] Filter only for GSK traits, reduces redundancy b/w GSK & Neale data.
#' @param neale_only [Default = FALSE] Filter only for Neale traits, reduces redundancy b/w GSK & Neale data.
#' @param db [Defalt = gtx::config_db()] Database to pull data from. 
#' @param impala [getOption("gtx.impala", NULL)] Implyr impala connection
#' @return data.frame with all coloc results in the region
#' @examples 
#' Query aba colocs:
#' colocs <- aba.query(hgncid = "HMGCR")
#' 
#' colocs <- aba.query(analysis_ids = "ukb_cool_analysis_id")
#' @export
#' @import tidyr
#' @import stringr
#' @import futile.logger
#' @import glue 
#' @import dplyr
aba.query <- function(analysis_ids, hgncid, ensemblid, rsid,
                      chrom, pos, pos_start, pos_end, p12_ge = 0.80, 
                      minpval1_le   = 1e-4, minpval2_le = 5e-8,
                      surround = 1e6, ncase_ge = 200, ncohort_ge = 200, 
                      protein_coding_only = FALSE, 
                      neale_only  = FALSE, 
                      gsk_only    = FALSE,
                      ttam_only   = FALSE,
                      db     = gtx::config_db(),
                      impala = getOption("gtx.impala", NULL)){
  # ---
  futile.logger::flog.debug("aba.query | validating input.")
  if(missing(hgncid) & 
     missing(ensemblid) & 
     missing(rsid) & 
     (missing(chrom) & ( (missing(pos_start) & missing(pos_end)) | missing(pos) ) ) & 
     missing(analysis_ids)){
    futile.logger::flog.error("aba.query | must specify either: analysis_ids, hgncid, ensemblid, or rsid. Skipping.")
    return()
  }
  # ---
  futile.logger::flog.debug("aba.query | validating impala connection.")
  impala <- validate_impala(impala = impala)
  
  # ---
  futile.logger::flog.debug("aba.query | establishing connection to database tables.")
  aba_tbl      <- dplyr::tbl(impala, glue::glue("{db}.coloc_results"))
  gwas_th_tbl  <- dplyr::tbl(impala, glue::glue("{db}.gwas_results_top_hits"))
  genes_tbl    <- dplyr::tbl(impala, glue::glue("{db}.genes"))
  analyses_tbl <- dplyr::tbl(impala, glue::glue("{db}.analyses"))
  sites_tbl    <- dplyr::tbl(impala, glue::glue("{db}.sites_ukb_500kv3"))
  
  # ---
  futile.logger::flog.debug("aba.query | set input")
  if(!missing(analysis_ids)){
    futile.logger::flog.debug("aba.query | processing analysis_ids input.")
    input <- dplyr::tibble(input = analysis_ids, analysis = analysis_ids)
  }
  else if(!missing(hgncid)){
    futile.logger::flog.debug("aba.query | processing hgncid input.")
    input <- dplyr::tibble(input = hgncid, hgncid = hgncid)
  }
  else if(!missing(ensemblid)){
    futile.logger::flog.debug("aba.query | processing ensemblid input.")
    input <- dplyr::tibble(input = ensemblid, ensemblid = ensemblid)
  }
  else if(!missing(rsid)){
    futile.logger::flog.debug("aba.query | processing rsid input.")
    input <- 
      dplyr::tibble(input = rsid) %>% 
      dplyr::mutate(rs = stringr::str_match(input, "\\d+") %>% as.numeric())
  }
  else if(!missing(chrom) & !missing(pos)){
    futile.logger::flog.debug("aba.query | processing chrom & pos input.")
    input <- 
      dplyr::tibble(chrom = chrom, in_start = pos - surround, in_end = pos + surround) %>% 
      dplyr::mutate(input = glue::glue("chr{chrom}:{in_start}-{in_end}") %>% as.character()) # have to use as.character for copy_to()
  }
  else if(!missing(chrom) & !missing(pos_start) & !missing(pos_end)){
    futile.logger::flog.debug("aba.query | processing chrom & pos input.")
    input <- 
      dplyr::tibble(chrom = chrom, in_start = pos_start, in_end = pos_end) %>% 
      dplyr::mutate(input = glue::glue("chr{chrom}:{in_start}-{in_end}") %>% as.character()) 
  }
  else {
    futile.logger::flog.error("aba.query | unable to properly handle input.")
    return()
  }
  
  # ---
  futile.logger::flog.debug("aba.query | copy input to RDIP")
  input_tbl <- impala_copy_to(df = input, dest = impala)
  
  # ---
  futile.logger::flog.debug("aba.query | harmonizing input")
  if(!missing(hgncid)){
    futile.logger::flog.debug("aba.query | harmonizing hgncid input.")
    input_tbl <- 
      dplyr::inner_join(
        input_tbl,
        genes_tbl %>% 
          dplyr::select(hgncid, ensemblid, chrom, in_pos = "pos_start"),
        by = "hgncid") %>% 
      dplyr::mutate(in_start = in_pos - surround, in_end = in_pos + surround) %>% 
      dplyr::select(-in_pos, -hgncid)
  }
  else if(!missing(ensemblid)){
    futile.logger::flog.debug("aba.query | harmonizing ensemblid input.")
    input_tbl <- 
      dplyr::inner_join(
        input_tbl,
        genes_tbl %>% 
          dplyr::select(ensemblid, chrom, in_pos = "pos_start"),
        by = "ensemblid") %>% 
      dplyr::mutate(in_start = in_pos - surround, in_end = in_pos + surround) %>% 
      dplyr::select(-in_pos)
  }
  else if(!missing(rsid)){
    futile.logger::flog.debug("aba.query | harmonizing rsid input.")
    input_tbl <- 
      dplyr::inner_join(
        input_tbl,
        sites_tbl,
        by = "rs") %>% 
      dplyr::select(input, chrom, in_pos = "pos") %>% 
      dplyr::mutate(in_start = in_pos - surround, in_end = in_pos + surround) %>% 
      dplyr::select(-in_pos)
  }
  else if(!missing(analysis_ids)){
    futile.logger::flog.debug("aba.query | harmonizing rsid input.")
    input_tbl <- 
      dplyr::inner_join(
        input_tbl,
        gwas_th_tbl,
        by = ("analysis")) %>% 
      dplyr::select(analysis, chrom, pos_index) %>% 
      dplyr::collect()
    
    input_tbl <- 
      input_tbl %>% 
      dplyr::rename(in_pos = "pos_index") %>% 
      dplyr::mutate(in_start = in_pos - surround, in_end = in_pos + surround) %>% 
      dplyr::mutate(input = glue::glue("{analysis}__chr{chrom}:{in_pos}") %>% as.character()) %>% 
      dplyr::select(input, chrom, in_start, in_end)
    
    input_tbl <- impala_copy_to(df = input_tbl, dest = impala)
  }
  
  # ---
  futile.logger::flog.debug("aba.query | prelim filter colocs")
  colocs_tbl <- 
    aba_tbl %>% 
    dplyr::filter(p12      >= p12_ge & 
           minpval1 <= minpval1_le &
           minpval2 <= minpval2_le) %>% 
    dplyr::inner_join(.,
               genes_tbl %>% 
                 dplyr::select(ensemblid, gene_start = pos_start, gene_end = pos_end, genetype, hgncid, chrom),
               by = c("entity" = "ensemblid"))
  
  if(protein_coding_only == TRUE){  
    colocs_tbl <- 
      colocs_tbl %>%
      dplyr::filter(genetype == "protein_coding")
  }
  
  # ---
  futile.logger::flog.debug("aba.query | append GWAS top hits & gene info")
  # Join GWAS top hits for each GWAS trait
  colocs_th_tbl <- 
    dplyr::inner_join(
      colocs_tbl,
      gwas_th_tbl %>% 
        dplyr::select(-signal, -num_variants, -rsq_index, -min_pval, 
                      -freq_index, -beta_index, -se_index),
      by = c("analysis2" = "analysis", "chrom")) %>% 
    dplyr::select(dplyr::everything(), 
                  th_pos   = pos_index, 
                  th_pval  = pval_index, 
                  th_start = pos_start, 
                  th_end   = pos_end,
                  th_ref   = ref_index,
                  th_alt   = alt_index) %>% 
    # Make sure the TH are in cis-windows of the coloc genes
    dplyr::filter((gene_start - 1e6 < th_start) & (gene_start + 1e6 > th_end)) %>%   
    # Join GWAS trait info - e.g. description & ncase
    dplyr::inner_join(.,
      analyses_tbl %>% dplyr::select(analysis, description, phenotype, ncase, ncohort),
      by = c("analysis2" = "analysis")) %>% 
    # Append RSID to each TH index
    dplyr::left_join(.,
      sites_tbl %>% dplyr::select("rs_chrom" = "chrom", "rs_pos" = "pos", "rs_ref" = "ref", "rs_alt" = "alt", "rs"),
      by = c("chrom" = "rs_chrom", "th_pos" = "rs_pos", "th_ref" = "rs_ref", "th_alt" = "rs_alt")) %>% 
  # we should then filter ncase
  dplyr::filter(ncase >= ncase_ge | ncohort >= ncohort_ge)
  
  # ---
  futile.logger::flog.debug("aba.query | filter colocs based on input")
  # If we gave a gene/region input filter based on that
  colocs_final_tbl <- 
    dplyr::inner_join(
      input_tbl,
      colocs_th_tbl,
      by = "chrom") %>% 
    # Make sure the TH are in the desired input window
    dplyr::filter((in_start < th_start) & (in_end > th_end)) 
    
  
  # ---
  futile.logger::flog.debug("aba.query | collect results")
  colocs_final <- 
    colocs_final_tbl %>% 
    dplyr::collect() %>% 
    dplyr::mutate(rsid = paste0("rs", rs)) %>% 
    dplyr::group_by(rsid) %>% 
    dplyr::mutate(genes_per_locus = dplyr::n_distinct(entity)) %>% 
    dplyr::ungroup()
  
  futile.logger::flog.debug("aba.query | verify collect results returned")
  if(nrow(colocs_final) == 0){
    futile.logger::flog.error("no results returned for the input.")
    return()
  }
  
  if(isTRUE(neale_only)){
    colocs_final <-
      colocs_final %>% 
      dplyr::filter(stringr::str_detect(analysis2, "neale"))
  }
  else if(isTRUE(gsk_only)){
    colocs_final <-
      colocs_final %>% 
      dplyr::filter(stringr::str_detect(analysis2, "GSK"))
  }
  else if(isTRUE(ttam_only)){
    colocs_final <-
      colocs_final %>% 
      dplyr::filter(stringr::str_detect(analysis2, "ttam_"))
  }
  
  # ---
  futile.logger::flog.debug("aba.query | clean up conn")
  close_int_conn(impala)
  
  # ---
  futile.logger::flog.debug("aba.query | complete")
  return(colocs_final)
}

#' \strong{aba.flatten() - Flatten full aba results}
#' 
#' The aba results and region queries will return colocs

#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param .data \code{\link{aba.query}} results object to filter
#' @return data.frame with input \code{\link{aba.query}} flattened to 1 gene / row. 
#' @examples 
#' Basic use:
#' colocs  <- aba.query(hgncid = "foo")
#' colocs_flat <- aba.flatten(colocs)
#' @export
#' @import tidyr
#' @import purrr
#' @import futile.logger
#' @import glue
#' @import dplyr
aba.flatten <- function(.data){
  input = .data
  # ---
  futile.logger::flog.debug("aba.flatten | verify input")
  mandatory_cols <- c("input", "rsid", "entity", "analysis1", "analysis2")
  if(any(purrr::map_lgl(mandatory_cols, ~ .x %in% names(input)) == FALSE)){
    futile.logger::flog.error("aba.flatten | missing input column.")
    return()
  }
  
  # ---
  futile.logger::flog.debug("aba.flatten | flatten data")
  input_filtered <-
    input %>%
    dplyr::group_by(input, rsid) %>% 
    dplyr::mutate(locus_n_genes     = dplyr::n_distinct(entity)) %>% 
    dplyr::mutate(locus_n_tissue    = dplyr::n_distinct(analysis1)) %>% 
    dplyr::group_by(input, entity) %>%
    dplyr::mutate(gene_n_top_hits   = dplyr::n_distinct(rsid)) %>%
    dplyr::mutate(gene_n_traits     = dplyr::n_distinct(analysis2)) %>% 
    dplyr::mutate(gene_n_tissue     = dplyr::n_distinct(analysis1)) %>% 
    dplyr::mutate(loci_min_n_genes  = min(locus_n_genes,  na.rm = TRUE) %>% as.integer()) %>% 
    dplyr::mutate(loci_max_n_genes  = max(locus_n_genes,  na.rm = TRUE) %>% as.integer()) %>% 
    dplyr::mutate(loci_min_n_tissue = min(locus_n_tissue, na.rm = TRUE) %>% as.integer()) %>%
    dplyr::mutate(loci_max_n_tissue = max(locus_n_tissue, na.rm = TRUE) %>% as.integer()) %>% 
    dplyr::ungroup() %>%
    dplyr::select(input, hgncid, entity,
           gene_n_tissue, gene_n_top_hits, gene_n_traits,
           loci_min_n_genes,  loci_max_n_genes, 
           loci_min_n_tissue, loci_max_n_tissue) %>% 
    dplyr::distinct() %>% 
    dplyr::arrange(gene_n_tissue, gene_n_top_hits, loci_max_n_genes, loci_min_n_genes)
  
  futile.logger::flog.debug("aba.flatten | return data")
  return(input_filtered)
}

#' \strong{aba.fill() - Fill in missing data to complete matrix of aba.query() results}
#' 
#' The \code{\link{aba.query}} will return only positive hits. Use 
#' this function to fill in the matrix of missing data based on positive hits.
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param .data 
#' @param db [Default = gtx::config_db()]
#' @return data.frame 
#' @examples 
#' Basic use:
#' colocs_pos  <- aba.query(hgncid = "foo")
#' colocs_full <- aba.fill(colocs_pos)
#' @export
#' @import tidyr
#' @import purrr
#' @import futile.logger
#' @import glue
#' @import dplyr
aba.fill <- function(.data, db = gtx::config_db(), impala = getOption("gtx.impala", NULL)){
  input = .data
  # ---
  futile.logger::flog.debug("aba.fill | verify input")
  mandatory_cols <- c("input", "rsid", "entity", "analysis1", "analysis2")
  if(any(purrr::map_lgl(mandatory_cols, ~ .x %in% names(input)) == FALSE)){
    futile.logger::flog.error("aba.fill | missing input column.")
    return()
  }
  # ---
  futile.logger::flog.debug("aba.fill | verify impala")
  impala <- validate_impala(impala = impala)
  # ---
  futile.logger::flog.debug("aba.fill | establishing connection to db.tables")
  aba_tbl <- dplyr::tbl(impala, glue::glue("{db}.coloc_results"))
  # ---
  futile.logger::flog.debug("aba.fill | expand matrix")
  data2pull <- 
    input %>% 
    dplyr::group_by(input) %>% 
    tidyr::expand(analysis2, analysis1, entity) %>% 
    dplyr::ungroup()
  # ---
  futile.logger::flog.debug("aba.fill | copy input to RDIP for join")
  data2pull_tbl <- impala_copy_to(df = data2pull, dest = impala)
  # ---
  futile.logger::flog.debug("aba.fill | query expanded matrix from aba results")
  expanded_dat <- 
    dplyr::inner_join(
      aba_tbl,
      data2pull_tbl,
      by = c("analysis2", "analysis1", "entity")) %>% 
    dplyr::collect() %>% 
    dplyr::inner_join(., 
      input %>% 
        dplyr::select(input, analysis2, rsid, th_start, th_end, th_pos, th_ref, th_alt, th_pval, 
                      description, phenotype, ncase, ncohort, in_start, in_end) %>% 
        dplyr::distinct(), 
      by = c("input", "analysis2")) %>% 
    # Make sure the TH are in the desired input window
    dplyr::filter((in_start < th_start) & (in_end > th_end)) %>% 
    dplyr::inner_join(.,
      input %>% 
        dplyr::select(input, entity, hgncid, chrom, gene_start, genetype) %>% 
        dplyr::distinct(), 
      by = c("input", "entity")) %>% 
    # Make sure the gene cis-windows are in the gwas TH window
    dplyr::filter((gene_start - 1e6 < th_start) & (gene_start + 1e6 > th_end)) 
  
  # ---
  futile.logger::flog.debug("aba.fill | clean up conn")
  close_int_conn(impala)
  
  # ---
  futile.logger::flog.debug("aba.fill | complete")
  return(expanded_dat)
  
}

#' aba.plot() - Plot aba results
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
#' @return ggplot2 object for viz and export.
#' @examples 
#' Basic use:
#' coloc_plot <- aba.plot(colocs)
#' coloc_plot <- aba.plot(colocs, max_dot_size = 3, title = "Cool plot")
#' 
#' View the returned object:
#' coloc_plot
#' 
#' Save the object to PDF:
#' ggsave(filename = "coloc_aba.pdf", 
#'        plot     = coloc_plot, 
#'        width    = 11,  
#'        height   = 8.5, 
#'        dpi      = 300)
#'        
#' Advanced use: query data, remove death related traits, filter, and then plot
#' colocs <- aba.wrapper(hgncid = "HMGCR")
#' colocs %>% dplyr::filter(input == "HMGCR") %>% pluck("figures", 1) + ggtitle("Best plot ever")
#' colocs %>% dplyr::filter(input == "HMGCR") %>% pluck("data", 1) %>% dplyr::filter(hgncid != "bad_gene") %>% aba.plot()
#' @export
#' @import dplyr
#' @import stringr
#' @import ggplot2
#' @import futile.logger
#' @import glue
aba.plot <- function(.data, ...){
  ## Make it accept a vector, return data frame of input <-> figures or just 1 figure
  ## Pull out plot function into internal fxn
  ## add option to return data nested. 
  # Verify input
  futile.logger::flog.debug("aba.plot | validating input")
  input <- .data
  required_cols <- c("analysis1", "analysis2", "description", "hgncid", "p12", "alpha21", "gene_start", "th_pos", "chrom", "rsid")
  if(!all(required_cols %in% (names(input)))){
    futile.logger::flog.error(paste0("aba.plot | input is missing required cols. Required cols include:", paste(required_cols, collapse = ", ")))
    stop();
  }
  
  # If we only have 1 "input" process as singular
  if(all((names(input) %in% "input")==FALSE)){
    futile.logger::flog.debug("aba.plot | process single input")
    ret = aba.int_coloc_plot(.data = input, ...) 
  }
  else {
    futile.logger::flog.debug("aba.plot | process multiple inputs")
    ret <- 
      input %>% 
      dplyr::group_by(input) %>% 
      tidyr::nest() %>% 
      dplyr::mutate(figures = purrr::map(data, aba.int_coloc_plot, ...))
  }
  futile.logger::flog.debug("aba.plot | complete")
  return(ret)
}

#' \strong{aba.int_coloc_plot() - Plot aba results}
#' 
#' This function is for internal aba plotting use. 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param p12_ge [Default = 0.80] Values below this will be colored as non-significant, while above have color for direction. 
#' @param max_dot_size [Default = 5] Adjust max size of dots.
#' @param title Set the plot title. NULL will remove title. 
#' @return ggplot2 object for viz and export.
#' @import dplyr
aba.int_coloc_plot <- function(.data, p12_ge = 0.80, max_dot_size = 5, title = NULL){
  futile.logger::flog.debug("aba.plot | validating input")
  input <- .data
  required_cols <- c("analysis1", "analysis2", "description", "hgncid", 
                     "p12", "alpha21", "gene_start", "th_pos", "chrom", "rsid")
  
  if(!all(required_cols %in% (names(input)))){
    futile.logger::flog.error(paste0("aba.plot | input is missing required cols. Required cols include:", paste(required_cols, collapse = ", ")))
    stop();
  }
  if(!is.numeric(max_dot_size)){
    futile.logger::flog.warn("aba.plot | max_dot_size parameter is not numeric.")
  }
  if(!is.null(title) & !is.character(title)){
    futile.logger::flog.warn("aba.plot | title is neither NULL or character string.")
  }
  # ---
  # Clean data (e.g. tissue & description), add simple direction of effect
  futile.logger::flog.debug("aba.plot | cleaing data for plotting")
  fig_dat <- 
    input %>% 
    # clean tissue names
    dplyr::mutate(analysis1 = stringr::str_replace_all(analysis1, "_", " ")) %>%
    dplyr::mutate(analysis1 = stringr::str_replace_all(analysis1, "gtex7", "")) %>%
    dplyr::mutate(analysis1 = stringr::str_to_title(analysis1)) %>% 
    dplyr::mutate(analysis1 = stringr::str_trim(analysis1)) %>% 
    # clean GWAS description names by removing "_" 
    dplyr::mutate(description = stringr::str_replace_all(description, "_", " ")) %>% 
    # Fill in missing hgncdids with ensemblids
    dplyr::mutate(hgncid = dplyr::case_when(!stringr::str_detect(hgncid, "\\w+") ~ entity,
                              is.na(hgncid)               ~ entity,
                              TRUE                        ~ hgncid)) %>% 
    # Remove Broad in description b/c only using broad data here
    dplyr::mutate(description = stringr::str_replace(description, "\\(UKB Broad\\)", "")) %>% 
    dplyr::mutate(description_wrap = stringr::str_wrap(description,  width = 60, exdent = 8)) %>%
    # Adjust alpha21 so that non-significant (H4 < 0.80) = NA 
    dplyr::mutate(direction = dplyr::case_when(p12 >= !! p12_ge & alpha21 >  0        ~ "Increased", 
                                 p12 >= !! p12_ge & alpha21 <  0        ~ "Decreased", 
                                 p12 >= !! p12_ge & alpha21 == 0        ~ "None",
                                 p12 <  !! p12_ge & is.numeric(alpha21) ~ "Not Sig.")) %>% 
    dplyr::mutate(rs_pos = glue::glue("{rsid}\n{chrom}:{th_pos}") %>% as.character())
  
  futile.logger::flog.debug("aba.plot | ordering results by chromosome - position coordinates")
  # Order genes (cols) and rows (rs ids) by position
  fig_dat <- within(fig_dat, {
    hgncid <- reorder(hgncid, gene_start)
    rs_pos <- reorder(rs_pos, th_pos) 
  })

  futile.logger::flog.debug("aba.plot | plotting data")
  fig <-
    fig_dat %>% 
    ggplot2::ggplot(ggplot2::aes(x = analysis1, y = description_wrap)) + 
    ggplot2::geom_point(ggplot2::aes(size = p12, fill = direction), shape = 21, color = "black", stroke = 0.2) +
    ggplot2::facet_grid(rs_pos ~ hgncid, scales = "free", space = "free", drop = TRUE, margins = FALSE) +
    ggplot2::scale_fill_manual(values = c("Decreased" = "blue", 
                                 "None"      = "white",
                                 "Increased" = "red",
                                 "Not Sig."  = "grey")) +
    ggplot2::scale_radius(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.00),
                 range  = c(1, max_dot_size)) + 
    ggplot2::guides(fill  = ggplot2::guide_legend(title = "Change in gene expression\nwith increased GWAS trait")) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(legend.position  = "bottom",
          legend.text      = ggplot2::element_text(size = 8),
          legend.title     = ggplot2::element_text(face = "bold"),
          axis.text.x      = ggplot2::element_text(size = 6, angle = 35, hjust = 1),
          axis.text.y      = ggplot2::element_text(size = 6, angle = 0 ),
          axis.line        = ggplot2::element_line(color = "black"),
          axis.title.x     = ggplot2::element_blank(),
          axis.title.y     = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_line(size = 0.25),
          panel.grid.minor = ggplot2::element_line(size = 0.25),
          strip.text.x     = ggplot2::element_text(face = "bold.italic"),
          strip.text.y     = ggplot2::element_text(size = 5, angle = 0),
          plot.title       = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::guides(size  = ggplot2::guide_legend(title = "Posterior probability\nof colocalization"))
  # ---
  # Add a title for T2 zoom
  if("tier2_zoom" %in% names(fig_dat)){
    t2_title <- 
      fig_dat %>% 
      dplyr::distinct(tier1_zoom, tier2_zoom)
    
    if(nrow(t2_title) > 1){
      futile.logger::flog.error("aba.plot | tier1/2 have too many variables.")
      stop()
    }
    if(pull(t2_title, tier1_zoom) == pull(t2_title, tier2_zoom)){
      fig <- fig + ggplot2::ggtitle(glue::glue("Tier 2 zoom on {glue::glue_collapse(t2_title$tier1_zoom, sep = ' & ')}"))
    }
    else {
      fig <- fig + ggplot2::ggtitle(glue::glue("Tier 1 zoom on {glue::glue_collapse(t2_title$tier1_zoom, sep = ' & ')}, and Tier 2 zoom on {glue::glue_collapse(t2_title$tier2_zoom, sep = ' & ')}"))  
    }
    
  }
  # Add a title for T1 zoom
  else if("tier1_zoom" %in% names(fig_dat)){
    t1_title <- 
      fig_dat %>% 
      dplyr::distinct(tier1_zoom)
    
    if(nrow(t1_title) > 1){
      futile.logger::flog.error("aba.plot | tier1 have too many variables.")
      stop()
    }
    fig <- fig + ggplot2::ggtitle(glue::glue("Tier 1 zoom on {t1_title}"))
  }
  # Add a title based on ARGS input
  else if(!is.null(title)){
    fig <- fig + ggplot2::ggtitle(title)
  }
  # ---
  # If we have >4 genes, flip gene names vertical
  if(fig_dat %>% dplyr::filter(p12 > 0.80) %>% dplyr::distinct(hgncid) %>% nrow() > 4){
    fig <- fig + ggplot2::theme(strip.text.x = ggplot2::element_text(angle = 90, size = 12))
  }
  # ---
  futile.logger::flog.debug("aba.plot | return plot")

  return(fig)
}

#' aba.wrapper - single function to query and plot the aba colocs
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @inheritParams aba.query
#' @return data.frame with the inputs used, all the data for each input, and default plots
#' @example 
#' colocs <- aba.wrapper(hgncid = "HMGCR")
#' @export
#' @import tidyr
#' @import stringr
#' @import futile.logger
#' @import glue 
#' @import dplyr
aba.wrapper <- 
  function(analysis_ids, hgncid, ensemblid, rsid,
           chrom, pos, pos_start, pos_end, 
           impala = getOption("gtx.impala", NULL), ...){
  # ---
  futile.logger::flog.debug("aba.wrapper | validating input.")
  if(missing(hgncid) & 
     missing(ensemblid) & 
     missing(rsid) & 
     (missing(chrom) & missing(pos_start) & missing(pos_end)) & 
     missing(analysis_ids)){
    futile.logger::flog.error("aba.wrapper | must specify either: analysis_ids, hgncid, ensemblid, or rsid. Skipping.")
    return()
  }
  # ---
  futile.logger::flog.debug("aba.wrapper | validating impala connection.")
  close_conn_later <- dplyr::case_when(is.null(impala) ~ TRUE, !is.null(impala) ~ FALSE)
  conn <- validate_impala(impala = impala)
  attr(conn, "internal_conn") <- FALSE
  
  # ---
  futile.logger::flog.debug("aba.wrapper | aba.query")
  if(!missing(analysis_ids)){
    futile.logger::flog.debug("aba.wrapper | processing analysis_ids input.")
    colocs <- aba.query(analysis_ids = analysis_ids, impala = conn, ...)
  }
  else if(!missing(hgncid)){
    futile.logger::flog.debug("aba.wrapper | processing hgncid input.")
    colocs <- aba.query(hgncid = hgncid, impala = conn, ...)
  }
  else if(!missing(ensemblid)){
    futile.logger::flog.debug("aba.wrapper | processing ensemblid input.")
    colocs <- aba.query(ensemblid = ensemblid, impala = conn, ...)
  }
  else if(!missing(rsid)){
    futile.logger::flog.debug("aba.wrapper | processing rsid input.")
    colocs <- aba.query(rsid = rsid, impala = conn, ...)
  }
  else if(!missing(chrom) & !missing(pos)){
    futile.logger::flog.debug("aba.wrapper | processing chrom & pos input.")
    colocs <- aba.query(chrom     = chrom, 
                        pos       = pos,
                        impala    = conn, ...)
  }
  else if(!missing(chrom) & (!missing(pos_start) | !missing(pos_end))){
    futile.logger::flog.debug("aba.wrapper | processing chrom & pos input.")
    colocs <- aba.query(chrom     = chrom, 
                        pos_start = pos_start, 
                        pos_end   = pos_end, 
                        impala    = conn, ...) 
  }
  else {
    futile.logger::flog.error("aba.wrapper | unable to properly handle input.")
    return()
  }
  # ---
  futile.logger::flog.debug("aba.wrapper | aba.fill")
  colocs <- aba.fill(colocs, impala = conn)
  
  # ---
  futile.logger::flog.debug("aba.wrapper | aba.plots")
  colocs <- aba.plot(colocs)
  
  # ---
  futile.logger::flog.debug("aba.wrapper | clean up conn")
  if(close_conn_later == TRUE){
    implyr::dbDisconnect(conn)
  }
  
  # ---
  futile.logger::flog.debug("aba.wrapper | complete")
  return(colocs)
  }

#' aba.save - save \code{aba.wrapper} figures
#' 
#' This function will iterate over all inputs in an aba.wrapper results 
#' and save each figures using the input to name the output figure. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param .data Results from aba.wrapper or aba.plot
#' @param path [Default = \code{getwd}] Specify the path to save the figures
#' @param suffix [Default = "aba-colocs"] File name = .data$input_$suffix.pdf
#' @param ... Parameters to pass to ggsave. e.g. dpi = 300
#' @examples 
#' Query aba colocs:
#' colocs <- aba.wrapper(hgncid = c("gene1", "gene2"))
#' 
#' The wraper returns 2 figures, aba.save() will save both figures
#' aba.save(colocs) 
#' @export
#' @import purrr
#' @import futile.logger
#' @import glue 
#' @import dplyr
#' @import ggplot2
aba.save <- function(.data, path = getwd(), suffix = "aba-colocs", ...){
  # ---
  futile.logger::flog.debug("aba.save | validate input")
  if(missing(.data)){futile.logger::flog.error("aba.save | missing input .data")}
  
  input = .data
  mandatory_cols <- c("input", "figures")
  if(any(purrr::map_lgl(mandatory_cols, ~ .x %in% names(input)) == FALSE)){
    futile.logger::flog.error("aba.save | missing input column.")
    return()
  }
  # ---
  futile.logger::flog.debug("aba.save | validate path")
  safe_system <- purrr::safely(base::system)
  
  exec <- safe_system(glue::glue("mkdir -p {path}"))
  if (!is.null(exec$error)){
    futile.logger::flog.error(glue::glue("aba.save | unable to create/validate {path}"))
    stop()
  }
  
  # ---
  futile.logger::flog.debug("aba.save | saving figures . . .")
  purrr::walk2(glue::glue("{path}/{input$input}_{suffix}.pdf"), input$figures, ggsave, ...)
  futile.logger::flog.debug("aba.save | saving figures complete")
}

#' aba.zoom - A function to help zoom in on specific data from aba.
#' 
#' This function is to enable users to zoom in on specific features
#' from the aba coloc data. The default aba plots can be very large and messy, 
#' and filtering these data can be difficult. To enable users to filter the data,
#' this function will help "zoom" in on specific features such as genes or GWAS. 
#' This function can be used 1-2 times depending on the desired level of zoom. 
#' 
#' Tier 1 zoom: All colocs for the input query, all genes with shared coloc traits (extended gene list), 
#' and all non-shared traits for the extended gene list.
#' 
#' Tier 2 zoom: All colocs for the input query, all genes with shared colocs. 
#' Excludes non-shared traits from the extended gene list. Note, Tier 2 zoom does not always 
#' further resolve data complexity.
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param hgncid HGNC symbol.
#' @param entity Typically ensemblid. 
#' @param analysis2 GWAS analysis ID
#' @param p12_ge [Default >= 0.80] This is the "H4" posterior probability cutoff 
#' @export
#' @import futile.logger
#' @import dplyr
#' @examples
#' Use aba.wrapper to run query and get data. 
#' colocs <- aba.wrapper(hgncid = "HMGCR", neale_only = TRUE, protein_coding_only = TRUE)
#' 
#' Tier 1 zoom
#' colocs %>% pluck("data", 1) %>% aba.zoom(hgncid = "HMGCR") %>% aba.plot()
#' 
#' Tier 2 zoom:
#' colocs %>% pluck("data", 1) %>% aba.zoom(hgncid = "HMGCR") %>% aba.zoom(hgncid = "HMGCR") %>% aba.plot()
aba.zoom <- function(.data, hgncid, entity, analysis2, p12_ge = 0.80){
  input <- .data
  # ---
  futile.logger::flog.debug("aba.zoom | validating args.")
  if(missing(hgncid) & missing(entity) & missing(analysis2)){
    futile.logger::flog.error("aba.zoom | must specify either:hgncid, entity, analysis2.")
    stop()
  }
  # ---
  futile.logger::flog.debug("aba.zoom | validating input.")
  required_cols <- c("analysis1", "analysis2", "hgncid", "entity", "p12")
  if(!all(required_cols %in% (names(input)))){
    futile.logger::flog.error(paste0("aba.zoom | input is missing required cols. Required cols include:", paste(required_cols, collapse = ", ")))
    stop();
  }
  if(!missing(hgncid)){
    if(nrow(filter(.data, hgncid == !!hgncid)) == 0){
      futile.logger::flog.error(glue::glue("aba.zoom | input data does not contain the input hgncid: {hgncid}"))
      stop()
    }
  }
  else if(!missing(entity)){
    if(nrow(filter(.data, entity == !!entity)) == 0){
      futile.logger::flog.error(glue::glue("aba.zoom | input data does not contain the input ensemblid: {ensemblid}"))
      stop()
    }
  }
  else if(!missing(analysis2)){
    if(nrow(filter(.data, analysis2 == !!analysis2)) == 0){
      futile.logger::flog.error(glue::glue("aba.zoom | input data does not contain the input analysis2: {analysis2}"))
      stop()
    }
  }
  # ---
  futile.logger::flog.debug("aba.zoom | Determining how to process data.")
  if("zoom2" %in% names(input)){
    futile.logger::flog.error("aba.zoom | previous zoom2 found, cannot further process.")
    stop()
  }
  else if("zoom1" %in% names(input)){
    futile.logger::flog.debug("aba.zoom | previous zoom1 found, tier 2 zoom processing")
    ret <- aba.int_zoom2(input, hgncid = hgncid, entity = entity, 
                         analysis2 = analysis2, p12_ge = p12_ge)
  } 
  else {
    futile.logger::flog.debug("aba.zoom | no previous zoom found, tier 1 zoom processing")
    ret <- aba.int_zoom1(input, hgncid = hgncid, entity = entity, 
                         analysis2 = analysis2, p12_ge = p12_ge)
  }
  # ---
  futile.logger::flog.debug("aba.zoom | processing complete.")
  return(ret)
}

#' aba.int_zoom1 - Internal fxn for \code{aba.zoom}
#' 
#' Tier 1 zoom
#' @import dplyr
aba.int_zoom1 <- function(.data, hgncid, entity, analysis2, p12_ge){
  # ---
  if(!missing(hgncid)){
    futile.logger::flog.debug(glue::glue("aba.zoom1 | Filtering using input hgncid: {hgncid}."))
    ret <- 
      .data %>% 
      dplyr::mutate(tier1_zoom = !!hgncid) %>% 
      dplyr::mutate(zoom1_pos_query = (hgncid == !!hgncid & p12 >= p12_ge)) %>% 
      dplyr::group_by(analysis2) %>% 
      dplyr::mutate(zoom1_pos_gwas = any(zoom1_pos_query))
  }
  else if(!missing(entity)){
    futile.logger::flog.debug(glue::glue("aba.zoom1 | Filtering using input entity: {entity}."))
    ret <- 
      .data %>% 
      dplyr::mutate(tier1_zoom = !!entity) %>% 
      dplyr::mutate(zoom1_pos_query = (entity == !!entity & p12 >= p12_ge)) %>% 
      dplyr::group_by(analysis2) %>% 
      dplyr::mutate(zoom1_pos_gwas = any(zoom1_pos_query))
  }
  else if(!missing(analysis2)){
    futile.logger::flog.debug(glue::glue("aba.zoom1 | Filtering using input analysis2: {analysis2}."))
    ret <- 
      .data %>% 
      dplyr::mutate(tier1_zoom = !!analysis2) %>% 
      dplyr::group_by(analysis2) %>%
      dplyr::mutate(zoom1_pos_gwas = (analysis2 == !!analysis2))
      
  }
  # ---
  ret <- 
    ret %>% 
    dplyr::group_by(hgncid) %>% 
    dplyr::mutate(zoom1_proxy_genes = any(zoom1_pos_gwas & p12 >= p12_ge)) %>% 
    dplyr::group_by(analysis2) %>% 
    dplyr::mutate(zoom1_proxy_gwas = any(zoom1_proxy_genes & p12 >= p12_ge)) %>% 
    dplyr::group_by(analysis1) %>% 
    dplyr::mutate(zoom1_proxy_tissues = any(zoom1_proxy_genes & zoom1_proxy_gwas & p12 >= p12_ge)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(zoom1 = (zoom1_proxy_genes & zoom1_proxy_gwas & zoom1_proxy_tissues)) %>% 
    dplyr::filter(zoom1) %>%
    dplyr::select(-dplyr::matches("zoom1_"))
  # ---
  
  return(ret)
}

#' aba.int_zoom2 - Internal fxn for \code{aba.zoom}
#' 
#' Tier 2 zoom
#' @import dplyr
aba.int_zoom2 <- function(.data, hgncid, entity, analysis2, p12_ge){
  # ---
  if(!missing(hgncid)){
    futile.logger::flog.debug(glue::glue("aba.zoom2 | Filtering using input hgncid: {hgncid}."))
    ret <- 
      .data %>% 
      dplyr::mutate(tier2_zoom = !!hgncid) %>% 
      dplyr::mutate(zoom2_pos_query = (hgncid == !!hgncid & p12 >= p12_ge)) %>% 
      dplyr::group_by(analysis2) %>% 
      dplyr::mutate(zoom2_pos_gwas = any(zoom2_pos_query))
      
  }
  else if(!missing(entity)){
    futile.logger::flog.debug(glue::glue("aba.zoom2 | Filtering using input entity: {entity}."))
    ret <- 
      .data %>% 
      dplyr::mutate(tier2_zoom = !!entity) %>% 
      dplyr::mutate(zoom2_pos_query = (entity == !!entity & p12 >= p12_ge)) %>% 
      dplyr::group_by(analysis2) %>% 
      dplyr::mutate(zoom2_pos_gwas = any(zoom2_pos_query))
  }
  else if(!missing(analysis2)){
    futile.logger::flog.debug(glue::glue("aba.zoom2 | Filtering using input analysis2: {analysis2}."))
    ret <- 
      .data %>% 
      dplyr::mutate(tier2_zoom = !!analysis2)
      dplyr::group_by(analysis2) %>%
      dplyr::mutate(zoom2_pos_gwas = (analysis2 == !!analysis2))
      
  }
  # ---
  ret <- 
    ret %>% 
    dplyr::group_by(hgncid) %>% 
    dplyr::mutate(zoom2_proxy_genes = any(zoom2_pos_gwas & p12 >= 0.80)) %>% 
    dplyr::group_by(analysis1) %>% 
    dplyr::mutate(zoom2_proxy_tissues = any(zoom2_proxy_genes & zoom2_pos_gwas & p12 >= 0.80)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(zoom2 = (zoom2_proxy_genes & zoom2_pos_gwas & zoom2_proxy_tissues)) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(zoom2) %>% 
    dplyr::select(-dplyr::matches("zoom2_"))
  # ---
  
  return(ret)
}
