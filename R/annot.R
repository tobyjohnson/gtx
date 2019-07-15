#' Query analysis, variant, gene etc. annotation
#'
#' Common queries of the annotation tables of the GWAS summary statistics.
#' 
#' All the main function arguments are optional, but at least one argument 
#' must be given in order to generate any results.
#'
#' @param x data frame, see Details
#' @param analysis Analysis identifier(s)
#' @param chrom Chromosome identifier(s) (string)
#' @param pos Position(s)
#' @param ref Reference allele(s)
#' @param alt Alternate allele(s)
#' @param rs dbSNP identifier (e.g. "rs123456")
#' @param hgnc HGNC gene symbol(s)
#' @param ensg Ensembl gene identifier(s) (e.g. "ENSG00000123456")
#' @param entity Entity identifier(s)
#' @param max_variants Maximum number of variants to return
#' @param max_genes Maximum number of genes to return
#' @param dbc Database connection (see \code{\link{gtxconnect}})
#'
#' @return
#'  Query results
#'
#' @author Toby Johnson \email{Toby.x.Johnson@gsk.com}
#' @export
annot <- function(x, 
                  analysis,
                  chrom, pos, ref, alt, rs,
                  hgnc, ensg,
                  entity, 
                  max_variants = 30L, max_genes = 3L,
                  dbc = getOption("gtx.dbConnection", NULL)) {
  gtxdbcheck(dbc)
  if (!missing(x)) {
    if (!is.data.frame(x)) {
      gtx_fatal_stop("annot(): x must be a dataframe")
    }
    ## If columns in x match legal annot() calls, consider recursive call
    ##   although if >1 argument, these should all be "vectorised" queries
    ## Note, no warnings or checks on number of rows returned if argument x was used, 
    ##   assume internal call or user knows what they are doing
    if ('analysis' %in% names(x)) {
      return(do.call(annot, x[ , 'analysis', drop = FALSE]))
    }
    if (any(c('chrom', 'pos', 'ref', 'alt', 'rs') %in% names(x))) {
      gtx_fatal_stop("annot(): vectorised sites queries not supported yet")
    }
    if (any(c('hgnc', 'ensg') %in% names(x))) {
      gtx_fatal_stop("annot(): vectorised genes queries not supported yet")
    }
    # if entity and entity_type present, use both since it may be more efficient
    if (all(c('entity', 'entity_type') %in% names(x))) {
      return(sqlWrapper(dbc, 
                        sprintf('SELECT * FROM entities WHERE %s', 
                                sql_where_x(x[ , c('entity', 'entity_type'), drop = FALSE])),
                        uniq = FALSE, zrok = TRUE))
    }
    # otherwise query entity only
    if ('entity' %in% names(x)) {
      return(sqlWrapper(dbc, 
                        sprintf('SELECT * FROM entities WHERE %s', 
                                sql_where_x(x[ , 'entity', drop = FALSE])),
                        uniq = FALSE, zrok = TRUE))
    }
    gtx_fatal_stop('annot(): could not determine query from names(x) = [{paste(names(x), collapse = ", ")}]')
  } else if (!missing(analysis)) {
    return(gtxanalyses(analysis = unique(analysis))) # short term fix, move gtxanalyses code here
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

## internal function, vectorized and general purpose version of gtxwhere, gtxwhat etc.
sql_where_x <- function(x) {
  if (!is.data.frame(x)) {
    gtx_fatal_stop("sql_where_x(): x must be a dataframe")
  }
  # define a list of legal column names and types
  col_types <- list(analysis = 'alphanum', 
                    entity = 'alphanum.-', 
                    entity_type = 'alphanum')
  # check all columns of are legal
  cols_illegal <- setdiff(names(x), names(col_types))
  if (length(cols_illegal) > 0L) {
    gtx_fatal_stop('sql_where_x(): x contains illegal columns [ {paste(cols_illegal, collapse = ", ")} ]')
  }
  # optimize query by removing duplicates
  x <- unique(x)
  # check no NAs
  rows_na <- apply(is.na(x), 1, any)
  if (any(rows_na)) {
    gtx_fatal_stop('sql_where_x(): x contains NAs')
  }
  for (col in names(col_types)) {
    if (col %in% names(x)) {
      x[[col]] <- paste0(col, '=\'', sanitize(x[[col]], type = col_types[[col]]), '\'')
    }
  }
  return(paste0('(', apply(x, 1, paste, collapse = ' AND '), ')', collapse = ' OR '))
}