#' Query analysis, variant, gene etc. annotation
#'
#' Common queries of the annotation tables of the GWAS summary statistics.
#' 
#' All the main function arguments are optional, but at least one argument 
#' must be given in order to generate any results.
#'
#' @param analysis Analysis identifier
#' @param chrom Chromosome identifier (string)
#' @param pos Position
#' @param ref Reference allele
#' @param alt Alternate allele
#' @param rs dbSNP identifier (e.g. "rs123456")
#' @param hgnc HGNC gene symbol
#' @param ensg Ensembl gene identifier (e.g. "ENSG00000123456")
#' @param max_variants Maximum number of variants to return
#' @param max_genes Maximum number of genes to return
#' @param dbc Database connection (see \code{\link{gtxconnect}})
#'
#' @return
#'  Query results
#'
#' @author Toby Johnson \email{Toby.x.Johnson@gsk.com}
#' @export
annot <- function(analysis,
                  chrom, pos, ref, alt, rs,
                  hgnc, ensg,
                  max_variants = 30L, max_genes = 3L,
                  dbc = getOption("gtx.dbConnection", NULL)) {
  gtxdbcheck(dbc)
  if (!missing(analysis)) {
    return(gtxanalyses(analysis = analysis)) # short term fix, move gtxanalyses code here
  }
  if (!missing(chrom) || !missing(pos) || !missing(ref) || !missing(alt) || !missing(rs)) {
    w1 <- gtxwhere(chrom = chrom, pos = pos, ref = ref, alt = alt, rs = rs)
    v1 <- sqlWrapper(dbc,
                     glue('SELECT chrom, pos, ref, alt, rsid AS rs FROM sites WHERE {w1};'),
                     uniq = FALSE, zrok = TRUE) # default uniq = TRUE
    ## convert internal representations to userspace
    v1 <- within(v1, rs <- paste0('rs', rs))
    if (nrow(v1) == 0L) {
      gtx_warn('No variants match query in TABLE sites: {w1}')
      return(v1)
    }
    if (nrow(v1) > max_variants) {
      gtx_warn('>{max_variants} variants match query in TABLE sites: {w1}')
      gtx_warn('First {max_variants} results returned; to get more, refine query or increase max_variants')
      return(head(v1, max_variants))
    }
    return(v1)
  }
  if (!missing(hgnc) || !missing(ensg)) {
    w1 <- gtxwhere(hgncid = hgnc, ensemblid = ensg)
    g1 <- sqlWrapper(dbc,
                     glue('SELECT chrom, pos_start, pos_end, hgncid AS hgnc, ensemblid AS ensg, genetype FROM genes WHERE {w1};'),
                     uniq = FALSE, zrok = TRUE) # default uniq = TRUE
    if (nrow(g1) == 0L) {
      gtx_warn('No genes match query in TABLE genes: {w1}')
      return(g1)
    }
    if (nrow(g1) > max_genes) {
      gtx_warn('>{max_genes} genes match query in TABLE genes: {w1}')
      gtx_warn('First {max_genes} results returned; to get more, refine query or increase max_genes')
      return(head(g1, max_genes))
    }
    return(g1)
  }  
  gtx_warn('Incomplete query specified (use more arguments)')
  return(NULL)
}
  