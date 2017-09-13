## plan to move gtxdbcheck to separate file,
## and make this file all non-gtx-specific SQL sanitation functions

gtxdbcheck <- function(dbc = getOption("gtx.dbConnection", NULL), verbose = FALSE) {
  if (! 'RODBC' %in% class(dbc)) stop('dbc does not appear to be a database connection (not of class RODBC)')
  dbs <- sqlQuery(dbc, 'SHOW DATABASES;')
  if (!is.data.frame(dbs)) stop('dbc does not appear to be an open database connection (SHOW DATABASES did not return a dataframe)')
  tables <- sqlQuery(dbc, 'SHOW TABLES;')
  if (!is.data.frame(tables)) stop('dbc does not appear to be an open database connection (SHOW TABLES did not return a dataframe)')
  ## could add other checks for existence and schema of the tables present
  if (verbose) {
    if ('analyses' %in% tables$name) {
      res <- sqlQuery(dbc, 'SELECT count(1) AS n FROM analyses')
      if (is.data.frame(res)) {
        message('TABLE analyses: ', prettyNum(res$n, big.mark = ',', scientific = FALSE), ' records')
      } else {
        stop('SQL error:\n', as.character(res))
      }
    } else {
      stop('dbc does not provide TABLE analyses')
    }
    res <- sqlQuery(dbc, 'SELECT DISTINCT results_db FROM ANALYSES')
    if (is.data.frame(res)) {
      res1 <- sanitize(res$results_db, type = 'alphanum')
      res1a <- res1 %in% dbs$name # identify whether user has access
      message('TABLE analyses: gwas_results in databases: ', 
              paste(res1, ifelse(res1a, '[OK]', '[no access]'), collapse = ', '))
      for (results_db in res1[res1a]) { 
        res2 <- sqlQuery(dbc, sprintf('SELECT count(1) AS n FROM %s.gwas_results', results_db))
        if (is.data.frame(res2)) {
          message('TABLE ', results_db, '.gwas_results: ', prettyNum(res2$n, big.mark = ',', scientific = FALSE), ' records')
        } else {
          stop('SQL error:\n', as.character(res2))
        }
      }
    } else {
      stop('SQL error:\n', as.character(res))
    }
  }
  return(TRUE)
}
      
## Two uses for this function:
##  resolve = TRUE, return the results_db for a single analysis, throw informative
##               error if user does not have access
##  resolve = FALSE, return a table of analysis,results_db,has_access
##
## included in sanitation because it's primarily to make access
## control generate comprehensible error messages
gtxanalysisdb <- function(analysis,
                          resolve = TRUE,
                          dbc = getOption("gtx.dbConnection", NULL)) {
    gtxdbcheck(dbc)
    dbs <- sqlQuery(dbc, 'SHOW DATABASES;')
    if (resolve) {
        res <- sqlWrapper(dbc,
                          sprintf('SELECT results_db FROM analyses WHERE %s;',
                                  gtxwhat(analysis = sanitize1(analysis, type = 'alphanum')))) # default uniq = TRUE 
        if (res$results_db %in% dbs$name) {
            return(res$results_db)
        } else {
            stop('analysis [ ', analysis, ' ] no access to required database [ ', res$resultsdb, ' ]')
        }
    } else {
        res <- sqlWrapper(dbc,
                          sprintf('SELECT analysis, results_db FROM analyses WHERE %s;',
                                  gtxwhat(analysis = sanitize(analysis, type = 'alphanum'))),
                          uniq = FALSE)
        res$has_access <- res$results_db %in% dbs$name
        return(res)
    }
}  

##
## sanitize x for use in SQL
## specialized for some types relevant to genetic data
## 
sanitize <- function(x, values, type) {
    ## function to sanitize x for use in SQL queries
    if (!missing(values)) {
        x <- as.character(na.omit(x)) # converts to character vector for (almost?) all possible inputs
        values <- as.character(na.omit(values))
        xm <- match(x, values)
        if (any(is.na(xm))) {
            stop('SQL input [ ', paste(x[is.na(xm)], collapse = ', '),
                 ' ] not in [ ', paste(values, collapse = ', '), ' ]')
        }
        return(values[xm])
    } else if (!missing(type)) {
        if (identical(type, "int")) {
            x <- na.omit(x)
            xi <- suppressWarnings(as.integer(x))
            if (any(is.na(xi))) {
                stop('SQL input [ ', paste(x[is.na(xi)], collapse = ', '),
                     ' ] not integer')
            }
            return(as.character(xi))
        } else if (identical(type, "double")) {
            x <- na.omit(x)
            xd <- suppressWarnings(as.double(x))
            if (any(is.na(xd))) {
                stop('SQL input [ ', paste(x[is.na(xd)], collapse = ', '),
                     ' ] not double')
            }
            return(as.character(xd))
        } else if (identical(type, "alphanum")) {
            ## note that here, "alphanum" means starting with alphabetic then alpha or numeric or underscore
            x <- as.character(na.omit(x))
            xa <- grepl("^[A-Za-z][A-Za-z0-9_]+$", x)
            if (any(!xa)) {
                stop('SQL input [ ', paste(x[!xa], collapse = ', '),
                     ' ] not alphanumeric')
            }
            return(x)
        } else if (identical(type, "alphanum-") || identical(type, "hgnc")) {
            ## "alphanum-" allows hyphens after first character
            ## hgnc is a synonym, for gene names e.g. HLA-DQB
            ## break up this case if broader pattern set needed for hgnc
            x <- as.character(na.omit(x))
            xa <- grepl("^[A-Za-z][A-Za-z0-9_-]+$", x)
            if (any(!xa)) {
                stop('SQL input [ ', paste(x[!xa], collapse = ', '),
                     ' ] not alphanumeric (hyphen allowed after first character)')
            }
            return(x)
        } else if (identical(type, "alphanum.")) {
            ## "alphanum." allows periods after first character, not sure what use cases are?
            x <- as.character(na.omit(x))
            xa <- grepl("^[A-Za-z][A-Za-z0-9_.]+$", x)
            if (any(!xa)) {
                stop('SQL input [ ', paste(x[!xa], collapse = ', '),
                     ' ] not alphanumeric (period allowed after first character)')
            }
            return(x)
        } else if (identical(type, "alphanum.-")) {
            ## "alphanum.-" allows periods and hyphens after first character, not sure what use cases are?
            x <- as.character(na.omit(x))
            xa <- grepl("^[A-Za-z][A-Za-z0-9_.-]+$", x)
            if (any(!xa)) {
                stop('SQL input [ ', paste(x[!xa], collapse = ', '),
                     ' ] not alphanumeric (period and hyphen allowed after first character)')
            }
            return(x)
        } else if (identical(type, "rs")) {
            x <- tolower(as.character(na.omit(x)))
            xrs <- grepl("^rs[1-9][0-9]*$", x)
            if (any(!xrs)) {
                stop('SQL input [ ', paste(x[!xrs], collapse = ', '),
                     ' ] not rs identifier(s)')
            }
            return(substr(x, 3, nchar(x))) # note rs parts of identifiers are stripped
        } else if (identical(type, "ENSG")) {
            x <- toupper(as.character(na.omit(x)))
            xens <- grepl("^ENSG[0-9]+$", x)
            if (any(!xens)) {
                stop('SQL input [ ', paste(x[!xens], collapse = ', '),
                     ' ] not ENSG identifier(s)')
            }
            return(x)
            # return(substr(x, 5, nchar(x))) # note in future ENSG parts of identifiers are stripped
        } else if (identical(type, "ACGT+")) {
            x <- as.character(na.omit(x))
            xa <- grepl("^[ACGT]+$", x)
            if (any(!xa)) {
                stop('SQL input [ ', paste(x[!xa], collapse = ', '),
                     ' ] not an ACGT sequence')
            }
            return(x)
        } else {
            stop("invalid type")
        }
    }
    stop("sanitize requires values or type")
}

##
## wrapper to throw error if length>1
## 
sanitize1 <- function(x, values, type) {
    x <- sanitize(x, values, type)
    if (identical(length(x), 1L)) return(x)
    if (identical(length(x), 0L)) {
        stop('SQL empty input where one is needed')
    }
    stop('SQL multiple inputs [ ', paste(x, collapse = ', '), ' ] where one is needed')
}

##
## wrapper for sqlQuery()
## checks return value is data frame with exactly
## [or at least] one row (if uniq is TRUE [or FALSE])
##
sqlWrapper <- function(dbc, sql, uniq = TRUE) {
    ## Note this function is for generic SQL RODBC usage
    ## and therefore does NOT take dbc from options('gtx.dbConnection')
    if (! 'RODBC' %in% class(dbc)) stop('dbc does not appear to be a database connection (not of class RODBC)')
    res <- sqlQuery(dbc, sql, as.is = TRUE)
    if (is.data.frame(res)) {
        if (identical(nrow(res), 0L)) {
            stop('SQL [ ', sql, ' ] returned 0 rows, expected ', if (uniq) '1 row' else '>=1 rows')
        } else if (identical(nrow(res), 1L)) {
            return(res)
        } else {
            if (uniq) {
                stop('SQL [ ', sql, ' ] returned ', nrow(res), 'rows, expected 1 row')
            } else {
                return(res)
            }
        }
    } else {
        stop('SQL [ ', sql, ' ] returned error:\n', as.character(res))
    }
    stop('internal error in sqlWrapper()') # should never happen
}

