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
#' @param protein_coding_only [Default = TRUE] Filter only for protein coding transcripts
#' @param neale_only [Default = FALSE] Filter onyl for Neale traits, reduces redundancy b/w GSK & Neale data.
#' @param gsk_only [Default = FALSE] Filter only for GSK traits, reduces redundancy b/w GSK & Neale data.
#' @param db [Defalt = "gene_gwas"] Database to pull data from. 
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
                      db     = "gene_gwas",
                      impala = getOption("gtx.impala", NULL)){
  ############################################
  flog.debug("aba.query | validating input.")
  if(missing(hgncid) & 
     missing(ensemblid) & 
     missing(rsid) & 
     (missing(chrom) & ( (missing(pos_start) & missing(pos_end)) | missing(pos) ) ) & 
     missing(analysis_ids)){
    flog.error("aba.query | must specify either: analysis_ids, hgncid, ensemblid, or rsid. Skipping.")
    return()
  }
  ############################################
  flog.debug("aba.query | validating impala connection.")
  impala <- validate_impala(impala = impala)
  
  ############################################
  flog.debug("aba.query | establishing connection to database tables.")
  aba_tbl      <- tbl(impala, glue("{db}.coloc_results"))
  gwas_th_tbl  <- tbl(impala, glue("{db}.gwas_results_top_hits"))
  genes_tbl    <- tbl(impala, glue("{db}.genes"))
  analyses_tbl <- tbl(impala, glue("{db}.analyses"))
  sites_tbl    <- tbl(impala, glue("{db}.sites_ukb_500kv3"))
  
  ############################################
  flog.debug("aba.query | set input")
  if(!missing(analysis_ids)){
    flog.debug("aba.query | processing analysis_ids input.")
    input <- tibble(input = analysis_ids, analysis = analysis_ids)
  }
  else if(!missing(hgncid)){
    flog.debug("aba.query | processing hgncid input.")
    input <- tibble(input = hgncid, hgncid = hgncid)
  }
  else if(!missing(ensemblid)){
    flog.debug("aba.query | processing ensemblid input.")
    input <- tibble(input = ensemblid, ensemblid = ensemblid)
  }
  else if(!missing(rsid)){
    flog.debug("aba.query | processing rsid input.")
    input <- 
      tibble(input = rsid) %>% 
      mutate(rs = str_match(input, "\\d+") %>% as.numeric())
  }
  else if(!missing(chrom) & !missing(pos)){
    flog.debug("aba.query | processing chrom & pos input.")
    input <- 
      tibble(chrom = chrom, in_start = pos - surround, in_end = pos + surround) %>% 
      mutate(input = glue("chr{chrom}:{in_start}-{in_end}") %>% as.character()) # have to use as.character for copy_to()
  }
  else if(!missing(chrom) & !missing(pos_start) & !missing(pos_end)){
    flog.debug("aba.query | processing chrom & pos input.")
    input <- 
      tibble(chrom = chrom, in_start = pos_start, in_end = pos_end) %>% 
      mutate(input = glue("chr{chrom}:{in_start}-{in_end}") %>% as.character()) 
  }
  else {
    flog.error("aba.query | unable to properly handle input.")
    return()
  }
  
  ############################################
  flog.debug("aba.query | copy input to RDIP")
  input_tbl <- impala_copy_to(df = input, dest = impala)
  
  ############################################
  flog.debug("aba.query | harmonizing input")
  if(!missing(hgncid)){
    flog.debug("aba.query | harmonizing hgncid input.")
    input_tbl <- 
      inner_join(
        input_tbl,
        genes_tbl %>% 
          select(hgncid, ensemblid, chrom, in_pos = "pos_start"),
        by = "hgncid") %>% 
      mutate(in_start = in_pos - surround, in_end = in_pos + surround) %>% 
      select(-in_pos, -hgncid)
  }
  else if(!missing(ensemblid)){
    flog.debug("aba.query | harmonizing ensemblid input.")
    input_tbl <- 
      inner_join(
        input_tbl,
        genes_tbl %>% 
          select(ensemblid, chrom, in_pos = "pos_start"),
        by = "ensemblid") %>% 
      mutate(in_start = in_pos - surround, in_end = in_pos + surround) %>% 
      select(-in_pos)
  }
  else if(!missing(rsid)){
    flog.debug("aba.query | harmonizing rsid input.")
    input_tbl <- 
      inner_join(
        input_tbl,
        sites_tbl,
        by = "rs") %>% 
      select(input, chrom, in_pos = "pos") %>% 
      mutate(in_start = in_pos - surround, in_end = in_pos + surround) %>% 
      select(-in_pos)
  }
  else if(!missing(analysis_ids)){
    flog.debug("aba.query | harmonizing rsid input.")
    input_tbl <- 
      inner_join(
        input_tbl,
        gwas_th_tbl,
        by = ("analysis")) %>% 
      select(analysis, chrom, pos_index) %>% 
      collect()
    
    input_tbl <- 
      input_tbl %>% 
      rename(in_pos = "pos_index") %>% 
      mutate(in_start = in_pos - surround, in_end = in_pos + surround) %>% 
      mutate(input = glue("{analysis}__chr{chrom}:{in_pos}") %>% as.character()) %>% 
      select(input, chrom, in_start, in_end)
    
    input_tbl <- impala_copy_to(df = input_tbl, dest = impala)
  }
  
  ############################################
  flog.debug("aba.query | prelim filter colocs")
  colocs_tbl <- 
    aba_tbl %>% 
    filter(p12      >= p12_ge & 
           minpval1 <= minpval1_le &
           minpval2 <= minpval2_le) %>% 
    inner_join(.,
               genes_tbl %>% 
                 select(ensemblid, gene_start = pos_start, gene_end = pos_end, genetype, hgncid, chrom),
               by = c("entity" = "ensemblid"))
  
  if(protein_coding_only == TRUE){  
    colocs_tbl <- 
      colocs_tbl %>%
      filter(genetype == "protein_coding")
  }
  
  ############################################
  flog.debug("aba.query | append GWAS top hits & gene info")
  # Join GWAS top hits for each GWAS trait
  colocs_th_tbl <- 
    inner_join(
      colocs_tbl,
      gwas_th_tbl %>% 
        select(-signal, -num_variants, -ref_index, 
               -alt_index, -rsq_index, -min_pval, 
               -freq_index, -beta_index, -se_index),
      by = c("analysis2" = "analysis", "chrom")) %>% 
    select(everything(), th_pos = pos_index, th_pval = pval_index, th_start = pos_start, th_end = pos_end) %>% 
    # Make sure the TH are in cis-windows of the coloc genes
    filter((gene_start - 1e6 < th_start) & (gene_start + 1e6 > th_end)) %>%   
    # Join GWAS trait info - e.g. description & ncase
    inner_join(
      .,
      analyses_tbl %>% select(analysis, description, phenotype, ncase, ncohort),
      by = c("analysis2" = "analysis")) %>% 
    # Append RSID to each TH index
    left_join(
      .,
      sites_tbl %>% select(rs_chrom = "chrom", rs_pos = "pos", rs),
      by = c("chrom" = "rs_chrom", "th_pos" = "rs_pos")) %>% 
  # we should then filter ncase
  filter(ncase >= ncase_ge | ncohort >= ncohort_ge)
  
  ############################################
  flog.debug("aba.query | filter colocs based on input")
  # If we gave a gene/region input filter based on that
  colocs_final_tbl <- 
    inner_join(
      input_tbl,
      colocs_th_tbl,
      by = "chrom") %>% 
    # Make sure the TH are in the desired input window
    filter((in_start < th_start) & (in_end > th_end)) 
    
  
  ############################################
  flog.debug("aba.query | collect results")
  colocs_final <- 
    colocs_final_tbl %>% 
    collect() %>% 
    mutate(rsid = paste0("rs", rs)) %>% 
    group_by(rsid) %>% 
    mutate(genes_per_locus = n_distinct(entity)) %>% 
    ungroup()
  
  flog.debug("aba.query | verify collect results returned")
  if(nrow(colocs_final) == 0){
    flog.error("no results returned for the input.")
    return()
  }
  
  if(isTRUE(neale_only)){
    colocs_final <-
      colocs_final %>% 
      filter(str_detect(analysis2, "neale"))
  }
  else if(isTRUE(gsk_only)){
    colocs_final <-
      colocs_final %>% 
      filter(str_detect(analysis2, "GSK"))
  }
  
  ############################################
  flog.debug("aba.query | clean up conn")
  close_int_conn(impala)
  
  ############################################
  flog.debug("aba.query | complete")
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
  ############################################
  flog.debug("aba.flatten | verify input")
  mandatory_cols <- c("input", "rsid", "entity", "analysis1", "analysis2")
  if(any(map_lgl(mandatory_cols, ~ .x %in% names(input)) == FALSE)){
    flog.error("aba.flatten | missing input column.")
    return()
  }
  
  ############################################
  flog.debug("aba.flatten | flatten data")
  input_filtered <-
    input %>%
    group_by(input, rsid) %>% 
    mutate(locus_n_genes     = n_distinct(entity)) %>% 
    mutate(locus_n_tissue    = n_distinct(analysis1)) %>% 
    group_by(input, entity) %>%
    mutate(gene_n_top_hits   = n_distinct(rsid)) %>%
    mutate(gene_n_traits     = n_distinct(analysis2)) %>% 
    mutate(gene_n_tissue     = n_distinct(analysis1)) %>% 
    mutate(loci_min_n_genes  = min(locus_n_genes,  na.rm = TRUE) %>% as.integer()) %>% 
    mutate(loci_max_n_genes  = max(locus_n_genes,  na.rm = TRUE) %>% as.integer()) %>% 
    mutate(loci_min_n_tissue = min(locus_n_tissue, na.rm = TRUE) %>% as.integer()) %>%
    mutate(loci_max_n_tissue = max(locus_n_tissue, na.rm = TRUE) %>% as.integer()) %>% 
    ungroup() %>%
    select(input, hgncid, entity,
           gene_n_tissue, gene_n_top_hits, gene_n_traits,
           loci_min_n_genes,  loci_max_n_genes, 
           loci_min_n_tissue, loci_max_n_tissue) %>% 
    distinct() %>% 
    arrange(gene_n_tissue, gene_n_top_hits, loci_max_n_genes, loci_min_n_genes)
  
  flog.debug("aba.flatten | return data")
  return(input_filtered)
}

#' \strong{aba.fill() - Fill in missing data to complete matrix of aba.query() results}
#' 
#' The \code{\link{aba.query}} will return only positive hits. Use 
#' this function to fill in the matrix of missing data based on positive hits.
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param .data 
#' @param db [Default = "gene_gwas"]
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
aba.fill <- function(.data, db = "gene_gwas", impala = getOption("gtx.impala", NULL)){
  input = .data
  ############################################
  flog.debug("aba.fill | verify input")
  mandatory_cols <- c("input", "rsid", "entity", "analysis1", "analysis2")
  if(any(map_lgl(mandatory_cols, ~ .x %in% names(input)) == FALSE)){
    flog.error("aba.fill | missing input column.")
    return()
  }
  ############################################
  flog.debug("aba.fill | verify impala")
  impala <- validate_impala(impala = impala)
  ############################################
  flog.debug("aba.fill | establishing connection to db.tables")
  aba_tbl <- tbl(impala, glue("{db}.coloc_results"))
  ############################################
  flog.debug("aba.fill | expand matrix")
  data2pull <- 
    input %>% 
    group_by(input) %>% 
    expand(analysis2, analysis1, entity) %>% 
    ungroup()
  ############################################
  flog.debug("aba.fill | copy input to RDIP for join")
  data2pull_tbl <- impala_copy_to(df = data2pull, dest = impala)
  ############################################
  flog.debug("aba.fill | query expanded matrix from aba results")
  expanded_dat <- 
    inner_join(
      aba_tbl,
      data2pull_tbl,
      by = c("analysis2", "analysis1", "entity")) %>% 
    collect() %>% 
    inner_join(
      ., 
      input %>% 
        select(input, analysis2, rsid, th_start, th_end, th_pos, th_pval, 
               description, phenotype, ncase, ncohort, in_start, in_end) %>% 
        distinct(), 
      by = c("input", "analysis2")) %>% 
    # Make sure the TH are in the desired input window
    filter((in_start < th_start) & (in_end > th_end)) %>% 
    inner_join(
      .,
      input %>% 
        select(input, entity, hgncid, chrom, gene_start, genetype) %>% 
        distinct(), 
      by = c("input", "entity")) %>% 
    # Make sure the gene cis-windows are in the gwas TH window
    filter((gene_start - 1e6 < th_start) & (gene_start + 1e6 > th_end)) 
  
  ############################################
  flog.debug("aba.fill | clean up conn")
  close_int_conn(impala)
  
  ############################################
  flog.debug("aba.fill | complete")
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
#' colocs %>% filter(input == "HMGCR") %>% pluck("figures", 1) + ggtitle("Best plot ever")
#' colocs %>% filter(input == "HMGCR") %>% pluck("data", 1) %>% filter(hgncid != "bad_gene") %>% aba.plot()
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
  flog.debug("aba.plot | validating input")
  input <- .data
  required_cols <- c("analysis1", "analysis2", "description", "hgncid", "p12", "alpha21", "gene_start", "th_pos", "chrom", "rsid")
  if(!all(required_cols %in% (names(input)))){
    flog.error(paste0("aba.plot | input is missing required cols. Required cols include:", paste(required_cols, collapse = ", ")))
    stop();
  }
  
  # If we only have 1 "input" process as singular
  if(all((names(input) %in% "input")==FALSE)){
    flog.debug("aba.plot | process single input")
    ret = aba.int_coloc_plot(.data = input, ...) 
  }
  else {
    flog.debug("aba.plot | process multiple inputs")
    ret <- 
      input %>% 
      group_by(input) %>% 
      nest() %>% 
      mutate(figures = map(data, aba.int_coloc_plot, ...))
  }
  flog.debug("aba.plot | complete")
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
  flog.debug("aba.plot | validating input")
  input <- .data
  required_cols <- c("analysis1", "analysis2", "description", "hgncid", "p12", "alpha21", "gene_start", "th_pos", "chrom", "rsid")
  if(!all(required_cols %in% (names(input)))){
    flog.error(paste0("aba.plot | input is missing required cols. Required cols include:", paste(required_cols, collapse = ", ")))
    stop();
  }
  if(!is.numeric(max_dot_size)){
    flog.warn("aba.plot | max_dot_size parameter is not numeric.")
  }
  if(!is.null(title) & !is.character(title)){
    flog.warn("aba.plot | title is neither NULL or character string.")
  }
  
  # Clean data (e.g. tissue & description), add simple direction of effect
  flog.debug("aba.plot | cleaing data for plotting")
  fig_dat <- 
    input %>% 
    # clean tissue names
    mutate(analysis1 = str_replace_all(analysis1, "_", " ")) %>%
    mutate(analysis1 = str_replace_all(analysis1, "gtex7", "")) %>%
    mutate(analysis1 = str_to_title(analysis1)) %>% 
    mutate(analysis1 = str_trim(analysis1)) %>% 
    # clean GWAS description names by removing "_" 
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
    mutate(rs_pos = glue("{rsid}\n{chrom}:{th_pos}") %>% as.character())
  
  flog.debug("aba.plot | ordering results by chromosome - position coordinates")
  # Order genes (cols) and rows (rs ids) by position
  fig_dat <- within(fig_dat, {
    hgncid <- reorder(hgncid, gene_start)
    rs_pos <- reorder(rs_pos, th_pos) 
  })
  
  flog.debug("aba.plot | plotting data")
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
          strip.text.x     = element_text(face = "bold.italic"),
          strip.text.y     = element_text(size = 5, angle = 0),
          plot.title       = element_text(hjust = 0.5)) +
    guides(size  = guide_legend(title = "Posterior probability\nof colocalization"))
  
  if(!is.null(title)){
    fig <- fig + ggtitle(title)
  }
  
  flog.debug("aba.plot | return plot")
  return(fig)
}

#' aba.wrapper - single function to query and plot the aba colocs
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @param analysis_ids ukbiobank analysis id(s), single string or vector.
#' @param hgncid HGNC symbol. Single string or vector. 
#' @param ensemblid Ensembl gene ID. Single string or vector.
#' @param rsid SNP rsid. Single string or vector.
#' @param surround [Default = 1e6] Distance from input gene/rsid to GWAS top hits. 
#' @param chrom chromosome - Used to define a specific region
#' @param pos_start start position - Used to define a specific region, overrides surround
#' @param pos_end end position - Used to define a specific region, overrides surround
#' @param p12_ge [Default >= 0.80] This is the "H4" posterior probability cutoff 
#' @param minpval1_le [Default <= 1e-4] Min pval seen in the eQTL                
#' @param minpval2_le [Default <= 5e-8] Min pval seen in the GWAS data           
#' @param ncase_ge [Default >= 200] Minimum ncases for traits.                   
#' @param ncohort_ge [Default >= 200] Minimum ncohort for traits.
#' @param protein_coding_only [Default = TRUE] Filter only for protein coding transcripts
#' @param neale_only [Default = FALSE] Filter onyl for Neale traits, reduces redundancy b/w GSK & Neale data.
#' @param gsk_only [Default = FALSE] Filter only for GSK traits, reduces redundancy b/w GSK & Neale data.
#' @return data.frame with the inputs used, all the data for each input, and default plots
#' @examples 
#' Query aba colocs:
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
  ############################################
  flog.debug("aba.wrapper | validating input.")
  if(missing(hgncid) & 
     missing(ensemblid) & 
     missing(rsid) & 
     (missing(chrom) & missing(pos_start) & missing(pos_end)) & 
     missing(analysis_ids)){
    flog.error("aba.wrapper | must specify either: analysis_ids, hgncid, ensemblid, or rsid. Skipping.")
    return()
  }
  ############################################
  flog.debug("aba.wrapper | validating impala connection.")
  close_conn_later <- case_when(is.null(impala) ~ TRUE, !is.null(impala) ~ FALSE)
  conn <- validate_impala(impala = impala)
  attr(conn, "internal_conn") <- FALSE
  
  ############################################
  flog.debug("aba.wrapper | aba.query")
  if(!missing(analysis_ids)){
    flog.debug("aba.wrapper | processing analysis_ids input.")
    colocs <- aba.query(analysis_ids = analysis_ids, impala = conn, ...)
  }
  else if(!missing(hgncid)){
    flog.debug("aba.wrapper | processing hgncid input.")
    colocs <- aba.query(hgncid = hgncid, impala = conn, ...)
  }
  else if(!missing(ensemblid)){
    flog.debug("aba.wrapper | processing ensemblid input.")
    colocs <- aba.query(ensemblid = ensemblid, impala = conn, ...)
  }
  else if(!missing(rsid)){
    flog.debug("aba.wrapper | processing rsid input.")
    colocs <- aba.query(rsid = rsid, impala = conn, ...)
  }
  else if(!missing(chrom) & !missing(pos)){
    flog.debug("aba.wrapper | processing chrom & pos input.")
    colocs <- aba.query(chrom     = chrom, 
                        pos       = pos,
                        impala    = conn, ...)
  }
  else if(!missing(chrom) & (!missing(pos_start) | !missing(pos_end))){
    flog.debug("aba.wrapper | processing chrom & pos input.")
    colocs <- aba.query(chrom     = chrom, 
                        pos_start = pos_start, 
                        pos_end   = pos_end, 
                        impala    = conn, ...)
  }
  else {
    flog.error("aba.wrapper | unable to properly handle input.")
    return()
  }
  ############################################
  flog.debug("aba.wrapper | aba.fill")
  colocs <- aba.fill(colocs, impala = conn)
  
  ############################################
  flog.debug("aba.wrapper | aba.plots")
  colocs <- aba.plot(colocs)
  
  ############################################
  flog.debug("aba.wrapper | clean up conn")
  if(close_conn_later == TRUE){
    implyr::dbDisconnect(conn)
  }
  
  ############################################
  flog.debug("aba.wrapper | complete")
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
  ############################################
  flog.debug("aba.save | validate input")
  if(missing(.data)){flog.error("aba.save | missing input .data")}
  
  input = .data
  mandatory_cols <- c("input", "figures")
  if(any(map_lgl(mandatory_cols, ~ .x %in% names(input)) == FALSE)){
    flog.error("aba.save | missing input column.")
    return()
  }
  ############################################
  flog.debug("aba.save | validate path")
  safe_system <- safely(system)
  
  exec <- safe_system(glue("mkdir -p {path}"))
  if (!is.null(exec$error)){
    flog.error(glue("aba.save | unable to create/validate {path}"))
    stop()
  }
  
  ############################################
  flog.debug("aba.save | saving figures . . .")
  walk2(glue("{path}/{input$input}_{suffix}.pdf"), input$figures, ggsave, ...)
  flog.debug("aba.save | saving figures complete")
}
