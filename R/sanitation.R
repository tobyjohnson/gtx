#' @import DBI

##' @export
## FIXME, change to required(?) argument odbc=...
gtxconnect <- function(dbc = dbConnect(odbc::odbc(), dsn = 'impaladsn'), 
                       use_database,
                       do_stop = TRUE, # deprecate this, it's cleaner for caller to use a tryCatch()
                       cache = TRUE, cache_analyses = TRUE, cache_genes = TRUE) {
  if (missing(use_database)) {
    stop('gtxconnect() requires use_database argument')
    ## FIXME instead don't execute USE statement
  }
  tmp <- try(eval(dbc))
  if (identical('try-error', class(tmp))) {
    if (do_stop) {
      stop(tmp) # verbatim error message
    } else {
      return(list(check = FALSE,
                  status = 'Error creating db connection'))
    }
  }
  options(gtx.dbConnection = tmp)
  tmp <- gtxdbcheck(check_databases = use_database, do_stop = do_stop) # overwrite previous value of tmp
  if (!do_stop) gtx_debug(tmp$status) # FIXME this was a temporary workaround
  if (do_stop || tmp$check) { # if do_stop, would have stopped previously instead of setting tmp$check=FALSE
    invisible(DBI::dbExecute(getOption('gtx.dbConnection'), 
                        sprintf('USE %s;', sanitize1(use_database, type = 'alphanum'))))
    tmp <- gtxdbcheck(do_stop = do_stop) # overwrite previous value of tmp
    if (!do_stop) gtx_debug(tmp$status) # FIXME this was a temporary workaround
  }
  
  if (cache) gtxcache(cache_analyses = cache_analyses, cache_genes = cache_genes)
  
  # Document DB connection & GTX version.
  futile.logger::flog.info(glue::glue("{gtx:::gtxversion()}"))
  futile.logger::flog.info(glue::glue("GTX connection established to: {use_database}"))
  
  return(tmp)
}

## gtxdbcheck(dbc), called immediately within every gtx function that uses a database connection
##     should be the most lightweight check possible (minimize actual queries against dbc)
##     should throw an error with stop() in the most common failure modes
##       especially those that would otherwise throw unintelligible error messages
##       e.g. dbc not initialized, server unresponsive/down/no authentication 
##
## supports several different uses:
## do_stop=TRUE, default, will throw error, for use in interactive code/notebooks
## do_stop=FALSE, will return a text string, for use in applications that will handle the error themselves
##
## check_databases is optional as this can be expensive
##
## tables_required has a sensible default for the main dbc, but may need to be reduced for
## checking a connection to a cache
##                                       
## added verbose option otherwise more useful logging messages get swamped

#' @import DBI
#' @import RSQLite

##' @export
gtxdbcheck <- function(dbc = getOption("gtx.dbConnection", NULL),
                       do_stop = TRUE,
                       check_databases,
                       tables_required = c('analyses', 'gwas_results', 'genes'), 
                       verbose = FALSE) {

  ## FIXME: make this work for other client/server db connections especially MySQL
  ## Noting several blocks below behave slightly differently depending on class(dbc)

  ##
  ## Check dbc is of a class we recognize
  ## Construct a string describing current database connection  
  if ('Impala' %in% class(dbc)) {
      ## FIXME for Kerberized environment, consider SELECT effective_user() as more relevant than clientside USER@HOST 
      currdb <- paste0(Sys.getenv('USER'), '@', Sys.getenv('HOSTNAME'), ' <Impala>')
  } else if ('SQLiteConnection' %in% class(dbc)) {
      currdb <- '<SQLite>'
  } else {
      if (do_stop) {
          stop('dbc class [ ', paste(class(dbc), collapse = ', '), ' ] not recognized')
      } else {
          return(list(check = FALSE,
                      status = paste0('dbc class [ ', paste(class(dbc), collapse = ', '), ' ] not recognized')))
      }
  }
      
  ##
  ## If check_databases requested and we know how, check we can access them
  ## (i.e. check whether a later USE statement should be successful)
  if (!missing(check_databases)) {
      if ('Impala' %in% class(dbc)) {
          if (verbose) gtx_debug('Trying query: SHOW DATABASES;')
          dbs <- try(DBI::dbGetQuery(dbc, 'SHOW DATABASES;'))
	  if (identical('try-error', class(dbs))) {
	      # FIXME grepl the content of dbs for 'Ticket expired' to help user
              if (do_stop) {
                  stop(currdb, ' is not an open database connection (SHOW DATABASES returned an error)')
              } else {
                  return(list(check = FALSE,
                              status = paste(currdb, 'error listing available databases')))
              }
	  }
          if (!is.data.frame(dbs)) {
              if (do_stop) {
                  stop(currdb, ' does not appear to be an open database connection (SHOW DATABASES did not return a dataframe)')
              } else {
                  return(list(check = FALSE,
                              status = paste(currdb, 'cannot list available databases')))
              }
          }
          if (!all(check_databases %in% dbs$name)) {
              if (do_stop) {
                  stop(currdb, ' cannot access database(s) [ ', paste(setdiff(check_databases, dbs$name), collapse = ', '), ' ]')
              } else {
                  return(list(check = FALSE,
                              status = paste(currdb, 'cannot access database(s) [ ', paste(setdiff(check_databases, dbs$name), collapse = ', '), ' ]')))
              }
          }
      } else if ('SQLiteConnection' %in% class(dbc)) {
          stop('check_databases not possible for SQLiteConnection')    
      } else {
          stop('gtxdbcheck internal error')
      }
      if (do_stop) {
          return(TRUE)
      } else {
          return(list(check = TRUE,
                      status = paste0(currdb, ' can access database(s) [ ', check_databases, ' ]')))
      }
  }
  
  ## Update string describing current database connection to include db name or SQLite path
  if ('Impala' %in% class(dbc)) {
      if (verbose) gtx_debug('Trying query: SELECT current_database();')
      scdb <- try(DBI::dbGetQuery(dbc, 'SELECT current_database();'))
      if (identical('try-error', class(scdb))) {
	# FIXME grepl the content of scdb for 'Ticket expired' to help user
        if (do_stop) {
          stop(currdb, ' is not an open database connection (SELECT current_database() returned an error)')
        } else {
          return(list(check = FALSE,
                      status = paste(currdb, 'error determining current database')))
        }
      }
      currdb <- paste0(Sys.getenv('USER'), '@', Sys.getenv('HOSTNAME'), ' <Impala:', 
                       paste(scdb[ , 'current_database()'], collapse = ', '), '>')    
  } else if ('SQLiteConnection' %in% class(dbc)) {
      currdb <- paste0('<SQLite:', dbc@dbname, '>')
  } else {
      stop('gtxdbcheck internal error')
  }  
    
  ##
  ## Always check tables since this is cheapest way we know to check
  ## db connection is working and pointing to the expected content
  if ('Impala' %in% class(dbc)) {
      if (verbose) gtx_debug('Trying query: SHOW TABLES;')
      tables <- try(DBI::dbGetQuery(dbc, 'SHOW TABLES;'))
      if (identical('try-error', class(tables))) {
	# FIXME grepl the content of tables for 'Ticket expired' to help user
        if (do_stop) {
          stop(currdb, ' is not an open database connection (SHOW TABLES returned an error)')
        } else {
          return(list(check = FALSE,
                      status = paste(currdb, 'error listing available tables')))
        }
      }
      if (!is.data.frame(tables)) {
          if (do_stop) {
              stop(currdb, ' does not appear to be an open database connection (SHOW TABLES did not return a dataframe)')
          } else {
              return(list(check = FALSE,
                          status = paste(currdb, 'cannot list available database tables')))
          }
      }
      tables <- tables$name
  } else if ('SQLiteConnection' %in% class(dbc)) {
      tables <- RSQLite::dbListTables(dbc)
  }

  if (!is.null(tables_required) && !all(tables_required %in% tables)) {
      if (do_stop) {
          stop(currdb, ' cannot access table(s) [ ', paste(setdiff(tables_required, tables), collapse = ', '), ' ]')
      } else {
          return(list(check = FALSE,
                 status = paste(currdb,'cannot access table(s) [ ', paste(setdiff(tables_required, tables), collapse = ', '), ' ]')))
      }
  }
  ##
  ## Checks completed okay, return
  if (do_stop) {
      return(TRUE)
  } else {
      return(list(check = TRUE,
                  status = paste(currdb, 'OK')))
  }

  ## could add other checks for existence and schema of the tables present
  ## could add options to attempt to fix common problems
  ## e.g. Impala connections where SELECT * LIMIT 10 fails, can try INVALIDATE METADATA
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

    ## Determine databases available
    if ('Impala' %in% class(dbc)) {
        ## Avoid frequently issuing SHOW DATABASES; commands
        dbs <- getOption('gtx.dbConnection_SHOW_DATABASES', NULL)
        if (is.null(dbs)) {
            gtxdbcheck(dbc)
            dbs <- sqlWrapper(dbc, 'SHOW DATABASES;', uniq = FALSE)
            options(gtx.dbConnection_SHOW_DATABASES = dbs)
        }
    } else if ('SQLiteConnection' %in% class(dbc)) {
        ## When dbc is an SQLite handle, there is no concept of a database
        dbs <- NULL
    } else {
        stop('dbc class [ ', paste(class(dbc), collapse = ', '), ' ] not recognized')
    }

    ## Check the connection we will use for TABLE analyses (not necessarily dbc)
    gtxdbcheck(getOption('gtx.dbConnection_cache_analyses', dbc), tables_required = 'analyses')

    if (resolve) {
        ## Query results_db for this analysis
        res <- getDataFromDB(connectionType = 'SQL',
                             connectionArguments = list(getOption('gtx.dbConnection_cache_analyses', dbc),
                                                        sprintf('SELECT results_db FROM analyses WHERE %s;',
                                                                gtxwhat(analysis1 = analysis)))) # sanitation by gtxwhat; default uniq = TRUE
        
        ## First, if null or empty string, return empty string (hence caller will use unqualified table name)
        ## (note, database nulls are returned as R NAs; following is safe because length(res$results_db)==1)
        if (is.na(res$results_db) || res$results_db == '') {
            return('')
        }
        ## Next, if results_db is specified but dbs is NULL, throw error
        ## (this could be undesirable if this function is called with a different dbc to what will be used
        ##  to actually query the results)
        if ('SQLiteConnection' %in% class(dbc)) {
            stop('analysis [ ', analysis, ' ] has results_db [ ', res$results_db, ' ] but this connection does not support databases')
        }
        ## Otherwise, check whether database is accessible
        if (res$results_db %in% dbs$name) {
            ## Return results_db with period appended
            ## (should not require sanitation because of %in% check, but sanitize anyway to be safe...)
            return(paste0(sanitize(res$results_db, type = 'alphanum'), '.'))
        } else {
            stop('analysis [ ', analysis, ' ] has results_db [ ', res$results_db, ' ] which is not a database this connectionhas access to')
        }
    } else {
        ## Query results_db for all analyses requested
        res <- getDataFromDB(connectionType = 'SQL',
                             connectionArguments = list(getOption('gtx.dbConnection_cache_analyses', dbc),
                                                        sprintf('SELECT analysis, results_db FROM analyses WHERE %s;',
                                                                gtxwhat(analysis = analysis)), # sanitation by gtxwhat
                                                        uniq = FALSE))
        
        ## If accessible databases known, add 'has_access' column
        if (!is.null(dbs)) res$has_access <- res$results_db %in% dbs$name
        return(res)
    }
}  

##
## sanitize x for use in SQL
## specialized for some types relevant to genetic data
##
## "alphanum" and friends are intended for matching symbols, e.g. db table names, gene names
## "text" is intended for matching free text
#' @export
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
#' @export
sanitize1 <- function(x, values, type) {
    x <- sanitize(x, values, type)
    if (identical(length(x), 1L)) return(x)
    if (identical(length(x), 0L)) {
        stop('SQL empty input where one is needed')
    }
    stop('SQL multiple inputs [ ', paste(x, collapse = ', '), ' ] where one is needed')
}

##
## wrapper for dbGetQuery()
## checks return value is data frame with exactly
## [or at least] one row (if uniq is TRUE [or FALSE])
## -- allow zero rows is zrok=TRUE
##
#' @import DBI
#' 
#' @export
sqlWrapper <- function(dbc, sql, uniq = TRUE, zrok = FALSE) {
    ## Note this function is for generic SQL usage
    ## and therefore does NOT take dbc from options('gtx.dbConnection')
    if (! 'Impala' %in% class(dbc) && ! 'SQLiteConnection' %in% class(dbc)) {
        stop('dbc is not an odbc database connection of class Impala or SQLiteConnection')
    }
    flog.debug(paste0(class(dbc), ' query: ', sql))
    t0 <- as.double(Sys.time())
    res <- DBI::dbGetQuery(dbc, sql) # !!! removed as.is=TRUE when switched from RODBC to odbc
    t1 <- as.double(Sys.time())
    flog.debug(paste0(class(dbc), ' query returned ', nrow(res), ' rows in ', round(t1 - t0, 3), 's.'))
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


# retrive data from a DB
# this function wraps around sqlWrapper and allows
# to implement different types of connections if needed
# connectionType (character): what kind of connection is required
# connectionArguments(list): the arguments required to query the connected database
# returns queried results from the database in the format specifeid by the appropriate call
getDataFromDB <- function(connectionType = 'SQL', connectionArguments){
  isValidType <- checkType(connectionType)
  if (!is.null(isValidType)){
    stop(isValidType)
  }
  connectionResults <- switch(connectionType,
                              'SQL' = do.call(sqlWrapper, connectionArguments),
                              stop('Unknown connection type'))
  return(connectionResults)
}


# check the type of connection passed to getDataFromDB()
checkType <- function(typeDB){
  isValidType <- NULL
  if (is.null(typeDB)){isValidType <- 'Connection type cannot be NULL'}
  if (!is.character(typeDB)){isValidType <- 'Connection type cannot must be a character vector'}
  if (length(typeDB) != 1){isValidType <- 'Please specify only one connection type'}
  
  return(isValidType)
}
