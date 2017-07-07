gtxdbcheck <- function(dbc = getOption("gtx.dbConnection", NULL), verbose = FALSE) {
  if (! 'RODBC' %in% class(dbc)) stop('dbc does not appear to be a database connection (not of class RODBC)')
  tables <- sqlQuery(dbc, 'SHOW TABLES;')
  if (!is.data.frame(tables)) stop('dbc does not appear to be an open database connection (SHOW TABLES did not return a dataframe)')
  ## could add other checks for existence and schema of the tables present
  if (verbose) {
    if ('gwas_results' %in% tables$name) {
      message(sprintf('%s GWAS analyses', prettyNum(sqlQuery(dbc, 'SELECT count(1) AS n FROM analyses')$n, big.mark = ',', scientific = FALSE)))
      message(sprintf('%s GWAS results', prettyNum(sqlQuery(dbc, 'SELECT count(1) AS n FROM gwas_results')$n, big.mark = ',', scientific = FALSE)))
    } else {
      stop('dbc does not provide TABLE gwas_results')
    }
  }
  return(TRUE)
}

sanitize <- function(x, values, type) {
  ## function to sanitize x in preparation for constructing SQL queries
  if (!missing(values)) {
    return(values[na.omit(match(x, values))])
  } else if (!missing(type)) {
    if (identical(type, "int")) {
      return(suppressWarnings(na.omit(as.integer(x))))
    } else if (identical(type, "alphanum")) {
      ## note alphanum means starting with alphabetic then alpha or numeric, define other types e.g. broader
      return(grep("^[A-Za-z][A-Za-z0-9_]+$", x, value = TRUE))
    } else {
      stop("invalid type")
    }
  } else {
    stop("sanitize requires values or type")
  }
}

gtxwhere <- function(chrom, 
                     pos, pos_ge, pos_le, 
                     pos_end_ge, pos_start_le, 
                     pos_start_ge, pos_end_le, 
                     hgncid, ensemblid,
                     tablename) {
  ## function to construct a WHERE string for constructing SQL queries
  ## Notes:
  ##  if arguments have length > 1, WHERE string OR's over values, producing e.g.
  ##      "hgncid='ABCA1' OR hgncid='ABCA2'
  ##  if multiple arguments, WHERE string AND's over arguments, producing e.g.
  ##      "chrom='1' AND pos_start>=123456 AND pos_end<=123789"

  ## For a query region of interest, to finding segments (e.g. genes, recombination rate) that:
  ## wholly or partially overlap, use pos_end_ge=query_start, pos_start_le=query_end
  ## wholly overlap, use pos_start_ge=query_start, pos_end_le=query_end

  if (!missing(tablename)) {
    tablename <- sanitize(tablename, type = 'alphanum')
    if (identical(length(tablename), 1L) && all(!is.na(tablename)) || !is.null(tablename)) {
      tablename <- paste0(tablename, '.')
    } else {
      tablename <- ''
    }
  } else {
    tablename <- ''
  }
  
  ws1 <- list(
              if (missing(chrom)) NULL
              else sprintf("chrom='%s'", sanitize(chrom[1], values = c(as.character(1:22), "X", "Y"))), 
              
              if (missing(pos)) NULL
              else sprintf("pos=%i", sanitize(pos, type = "int")), 
              
              if (missing(pos_ge)) NULL
              else sprintf("pos>=%i", sanitize(pos_ge, type = "int")),
              
              if (missing(pos_le)) NULL
              else sprintf("pos<=%i", sanitize(pos_le, type = "int")),
              
              if (missing(pos_end_ge)) NULL
              else sprintf("pos_end>=%i", sanitize(pos_end_ge, type = "int")),
              
              if (missing(pos_start_le)) NULL
              else sprintf("pos_start<=%i", sanitize(pos_start_le, type = "int")),
              
              if (missing(pos_start_ge)) NULL
              else sprintf("pos_start>=%i", sanitize(pos_start_ge, type = "int")),
              
              if (missing(pos_end_le)) NULL
              else sprintf("pos_end<=%i", sanitize(pos_end_le, type = "int")),
              
              if (missing(hgncid)) NULL
              else sprintf("hgncid='%s'", sanitize(hgncid, type = "alphanum")),

              if (missing(ensemblid)) NULL
              else sprintf("ensemblid='%s'", sanitize(ensemblid, type = "alphanum"))
              )
  ws2 <- paste0("(", 
                unlist(sapply(ws1, function(x) if (is.null(x)) NULL else paste0(tablename, x, collapse = " OR "))), 
                ")", collapse = " AND ")
  return(ws2)
}
