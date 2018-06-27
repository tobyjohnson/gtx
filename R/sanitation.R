# plan to move gtxdbcheck to separate file,
## and make this file all non-gtx-specific SQL sanitation functions

gtxdbcheck <- function(dbc = getOption("gtx.dbConnection", NULL), verbose = FALSE) {
  if (! 'Impala' %in% class(dbc)) stop('dbc does not appear to be an Impala connection (not of class Impala)')
  dbs <- dbGetQuery(dbc, 'SHOW DATABASES;')
  if (!is.data.frame(dbs)) stop('dbc does not appear to be an open database connection (SHOW DATABASES did not return a dataframe)')
  tables <- dbGetQuery(dbc, 'SHOW TABLES;')
  if (!is.data.frame(tables)) stop('dbc does not appear to be an open database connection (SHOW TABLES did not return a dataframe)')
  ## could add other checks for existence and schema of the tables present
  if (verbose) {
    if ('analyses' %in% tables$name) {
      res <- dbGetQuery(dbc, 'SELECT count(1) AS n FROM analyses')
      if (is.data.frame(res)) {
        message('TABLE analyses: ', prettyNum(res$n, big.mark = ',', scientific = FALSE), ' records')
      } else {
        stop('SQL error:\n', as.character(res))
      }
    } else {
      stop('dbc does not provide TABLE analyses')
    }
    res <- dbGetQuery(dbc, 'SELECT DISTINCT results_db FROM ANALYSES')
    if (is.data.frame(res)) {
      res1 <- sanitize(res$results_db, type = 'alphanum')
      res1a <- res1 %in% dbs$name # identify whether user has access
      message('TABLE analyses: gwas_results in databases: ', 
              paste(res1, ifelse(res1a, '[OK]', '[no access]'), collapse = ', '))
      for (results_db in res1[res1a]) { 
        res2 <- dbGetQuery(dbc, sprintf('SELECT count(1) AS n FROM %s.gwas_results', results_db))
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
    dbs <- sqlWrapper(dbc, 'SHOW DATABASES;', uniq = FALSE)
    if (resolve) {
        res <- sqlWrapper(dbc,
                          sprintf('SELECT results_db FROM analyses WHERE %s;',
                                  gtxwhat(analysis1 = analysis))) # sanitation by gtxwhat; default uniq = TRUE 
        if (res$results_db %in% dbs$name) {
            ## should not require sanitation because of %in% check, but sanitize anyway to be safe...
            return(sanitize(res$results_db, type = 'alphanum'))
        } else {
            stop('analysis [ ', analysis, ' ] no access to required database [ ', res$resultsdb, ' ]')
        }
    } else {
        res <- sqlWrapper(dbc,
                          sprintf('SELECT analysis, results_db FROM analyses WHERE %s;',
                                  gtxwhat(analysis = analysis)), # sanitation by gtxwhat
                          uniq = FALSE)
        res$has_access <- res$results_db %in% dbs$name
        return(res)
    }
}  

##
## sanitize x for use in SQL
## specialized for some types relevant to genetic data
##
## "alphanum" and friends are intended for matching symbols, e.g. db table names, gene names
## "text" is intended for matching free text
sanitize <- function(x, values, type) {
    ## function to sanitize x for use in SQL queries
    ## FIXME for more general usage, the error message 'SQL input' may not always be appropriate
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
            x <- na.omit(x) ## silently drop missing values
            xi <- suppressWarnings(as.integer(x))
            if (any(is.na(xi))) {
                stop('SQL input [ ', paste(x[is.na(xi)], collapse = ', '),
                     ' ] not integer')
            }
            return(as.character(xi)) ## FIXME this could potentially render some ints as 1.234e7 for example, depending on behaviour of as.character() for input int type
        } else if (identical(type, "double")) {
            x <- na.omit(x)
            xd <- suppressWarnings(as.double(x))
            if (any(is.na(xd))) {
                stop('SQL input [ ', paste(x[is.na(xd)], collapse = ', '),
                     ' ] not double')
            }
            return(as.character(xd))
        } else if (identical(type, "count")) {
            ## "count" means a counting integer starting at 0, more strict about scientific notation than "int" above
            x <- as.character(na.omit(x))
            xa <- grepl("^[0-9]*$", x)
            if (any(!xa)) {
                stop('SQL input [ ', paste(x[!xa], collapse = ', '),
                     ' ] not counting integer')
            }
            return(x)
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
        } else if (identical(type, "text")) {
            ## intended for matching free text query terms
            x <- as.character(na.omit(x))
            xa <- grepl("^[A-Za-z0-9 ]*$", x)
            if (any(!xa)) {
                stop('SQL input [ ', paste(x[!xa], collapse = ', '),
                     ' ] not text (strict definition)')
            }
            return(x)
        } else if (identical(type, "rs")) {
            x <- tolower(as.character(na.omit(x)))
            xrs <- grepl("^rs[1-9][0-9]*$", x)
            if (any(!xrs)) {
                stop('SQL input [ ', paste(x[!xrs], collapse = ', '),
                     ' ] not rs identifier(s)')
            }
            return(substr(x, 3, nchar(x))) # note rs parts of identifiers are stripped AND STRINGS ARE RETURNED
        } else if (identical(type, "ENSG")) {
            x <- toupper(as.character(na.omit(x)))
            xens <- grepl("^ENSG[0-9]+$", x)
            if (any(!xens)) {
                stop('SQL input [ ', paste(x[!xens], collapse = ', '),
                     ' ] not ENSG identifier(s)')
            }
            return(x)
            # return(substr(x, 5, nchar(x))) # note in future ENSG parts of identifiers WILL BE stripped
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
## -- allow zero rows is zrok=TRUE
##
sqlWrapper <- function(dbc, sql, uniq = TRUE, zrok = FALSE) {
    ## Note this function is for generic SQL usage
    ## and therefore does NOT take dbc from options('gtx.dbConnection')
    if (! 'Impala' %in% class(dbc) && ! 'SQLiteConnection' %in% class(dbc)) {
        stop('dbc does not appear to be an Impala connection (not of class Impala or SQLiteConnection)')
    }
    res <- dbGetQuery(dbc, sql) # !!! removed as.is=TRUE when switched from RODBC to odbc
    if (is.data.frame(res)) {
        if (identical(nrow(res), 0L)) {
            if (zrok) return(res)
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

