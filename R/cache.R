#' Speed up search of genes and analyses
#' 
#' Create a in-memory SQLite db storing analyses metadata and genes coordinates, to speed up searches
#'
#' @param	cache_analysis	Cache the Analyses table
#' @param	cache_genes	Cache Gene coordinates
#' @param	disconnect	Disconnects the db. E.g. run gtxcache(disconnect=T) to free up memory
#' @param	dbc_from	DB to be cached (e.g. gwases)
#' @param 	dbc_to		If there is already a gtxcache database available, use this instead of creating a new one
#' 
#' @export
gtxcache <- function(cache_analyses = TRUE, 
                     cache_genes = TRUE,  
                     disconnect = FALSE, 
                     dbc_from = getOption("gtx.dbConnection", NULL), 
                     dbc_to = getOption('gtx.dbConnection_cache', NULL)) {

    ## create connection to in-memory SQLite db where the cache will be kept,
    ## if it hasn't already been created
    if (is.null(dbc_to)) {
        futile.logger::flog.info('Creating new in-memory SQLite db')
        dbc_to <- dbConnect(RSQLite::SQLite(), ':memory:')
        options(gtx.dbConnection_cache = dbc_to)
    } else {
        futile.logger::flog.info('Using existing in-memory SQLite db')
    }
    ## we should check this is a valid dbc

    if (disconnect) {
        futile.logger::flog.info('Disconnecting from in-memory SQLite db')
        dbDisconnect(getOption('gtx.dbConnection_cache'))
        options(gtx.dbConnection_cache = NULL, 
                gtx.dbConnection_cache_analyses = NULL, 
                gtx.dbConnection_cache_genes = NULL)
        return(invisible(NULL))
    }
 
    gtxdbcheck(dbc_from) # this is the *source* database that we are caching content from
   
    if (cache_analyses) {
        futile.logger::flog.info('Caching TABLE analyses ...')
        t0 <- as.double(Sys.time())
        DBI::dbWriteTable(getOption('gtx.dbConnection_cache'),
                     'analyses',
                     sqlWrapper(dbc_from, 'SELECT * FROM analyses;', uniq = FALSE),
                         overwrite = TRUE)
        DBI::dbExecute(getOption('gtx.dbConnection_cache'),
                  'CREATE UNIQUE INDEX analyses_analysis ON analyses (analysis)') # should be unique
        t1 <- as.double(Sys.time())
        futile.logger::flog.info(paste0('Cached and INDEXed TABLE analyses with ',
                    sqlWrapper(getOption('gtx.dbConnection_cache'),
                            'SELECT count(1) AS n FROM analyses')$n,
                    ' rows in ', round(t1 - t0, 3), 's.'))
        futile.logger::flog.info('Caching TABLE analyses_tags ...')
        t0 <- as.double(Sys.time())
        DBI::dbWriteTable(getOption('gtx.dbConnection_cache'),
                          'analyses_tags',
                          sqlWrapper(dbc_from, 'SELECT * FROM analyses_tags;', uniq = FALSE),
                          overwrite = TRUE)
        DBI::dbExecute(getOption('gtx.dbConnection_cache'),
                       'CREATE INDEX analyses_tags_analysis ON analyses_tags (analysis)') # may not be unique
        t1 <- as.double(Sys.time())
        futile.logger::flog.info(paste0('Cached and INDEXed TABLE analyses_tags with ',
                  sqlWrapper(getOption('gtx.dbConnection_cache'),
                             'SELECT count(1) AS n FROM analyses_tags')$n,
                  ' rows in ', round(t1 - t0, 3), 's.'))
        options(gtx.dbConnection_cache_analyses = getOption('gtx.dbConnection_cache'))
    }
    
    if (cache_genes) {
        futile.logger::flog.info('Caching TABLE genes ...')
        t0 <- as.double(Sys.time())
        DBI::dbWriteTable(getOption('gtx.dbConnection_cache'),
                     'genes',
                     sqlWrapper(dbc_from, 'SELECT * FROM genes;', uniq = FALSE),
                         overwrite = TRUE)
        DBI::dbExecute(getOption('gtx.dbConnection_cache'),
                  'CREATE INDEX genes_chrom_pos_start_pos_end ON genes (chrom, pos_start, pos_end)')
        t1 <- as.double(Sys.time())
        futile.logger::flog.info(paste0('Cached and INDEXed TABLE genes with ',
                    sqlWrapper(getOption('gtx.dbConnection_cache'),
                            'SELECT count(1) AS n FROM genes')$n,
                    ' rows in ', round(t1 - t0, 3), 's.'))
        options(gtx.dbConnection_cache_genes = getOption('gtx.dbConnection_cache'))
    }
    
    return(invisible(NULL))
}

