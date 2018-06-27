## for debugging, we have a function to print messages depending on a global option gtx.debug
gtxlog <- function(...) {
    if (getOption('gtx.debug', FALSE)) return(message(...))
    return(invisible(NULL))
}

##
## convenience function to construct WHERE
## part of SQL for genomic data tables
##
gtxwhere <- function(chrom, 
                     pos, pos_ge, pos_le, 
                     pos_end_ge, pos_start_le, 
                     pos_start_ge, pos_end_le, 
                     rs, hgncid, ensemblid,
                     entity, 
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
        tablename <- paste0(sanitize1(tablename, type = 'alphanum'), '.')
    } else {
        tablename <- ''
    }
    
    ws1 <- list(
        if (missing(chrom)) NULL
        else sprintf("chrom='%s'", sanitize1(chrom, values = c(as.character(1:22), "X", "Y"))),
        
        if (missing(pos)) NULL
        else sprintf("pos=%s", sanitize(pos, type = "int")), 
        
        if (missing(pos_ge)) NULL
        else sprintf("pos>=%s", sanitize(pos_ge, type = "int")),
        
        if (missing(pos_le)) NULL
        else sprintf("pos<=%s", sanitize(pos_le, type = "int")),
        
        if (missing(pos_end_ge)) NULL
        else sprintf("pos_end>=%s", sanitize(pos_end_ge, type = "int")),
        
        if (missing(pos_start_le)) NULL
        else sprintf("pos_start<=%s", sanitize(pos_start_le, type = "int")),
        
        if (missing(pos_start_ge)) NULL
        else sprintf("pos_start>=%s", sanitize(pos_start_ge, type = "int")),
        
        if (missing(pos_end_le)) NULL
        else sprintf("pos_end<=%s", sanitize(pos_end_le, type = "int")),
        
        if (missing(rs)) NULL
        else sprintf("rsid=%s", sanitize(rs, type = "rs")), # note sanitize(type="rs") strips the rs parts hence returns integers
        
        if (missing(hgncid)) NULL
        else sprintf("hgncid='%s'", sanitize(hgncid, type = "alphanum-")), # thousands of HGNC ids contain hyphens, e.g. HLA-A
        
        if (missing(ensemblid)) NULL
        else sprintf("ensemblid='%s'", sanitize(ensemblid, type = "ENSG")),

        if (missing(entity)) NULL
        else sprintf("feature='%s'", sanitize(entity, type = "alphanum")) # will be entity=INT FIXME
    )
    ws2 <- paste0("(", 
                  unlist(sapply(ws1, function(x) if (is.null(x)) NULL else paste0(tablename, x, collapse = " OR "))), 
                  ")", collapse = " AND ")
    return(ws2)
}

##
## convenience function to construct WHERE
## part of SQL for analyses table
## Note behaviour for most arguments here is OR/OR, different to gtxwhere()
##
gtxwhat <- function(analysis1,
                    analysis,
                    analysis_not, 
                    description_contains,
                    phenotype_contains,
                    has_tag, 
                    ncase_ge,
                    ncohort_ge,
                    tablename) {
  ## function to construct a WHERE string for constructing SQL queries
  ## Notes:
  ##  if arguments have length > 1, WHERE string OR's over values
  ##  if multiple arguments, WHERE string OR's over arguments

    if (!missing(tablename)) {
        tablename <- paste0(sanitize1(tablename, type = 'alphanum'), '.')
    } else {
        tablename <- ''
    }

    ## use of analysis1 is a special case where exactly one analysis key is being provided
    ## ignore all other arguments
    if (!missing(analysis1)) {
        if (getOption('gtx.analysisIsString', TRUE)) { # Future release will change default to FALSE
            return(sprintf("(%sanalysis='%s')", tablename, sanitize1(analysis1, type = "alphanum")))
        } else {
            return(sprintf("(%sanalysis=%s)", tablename, sanitize1(analysis1, type = "count"))) # note no quotes
        }
    }
    
    ##
    ## analysis, description_contains, phenotype_contains etc. are OR'ed within and between arguments
    ws1 <- list(
        if (missing(analysis)) NULL
        else {
            if (getOption('gtx.analysisIsString', TRUE)) { # Future release will change default to FALSE
                sprintf("analysis='%s'", sanitize(analysis, type = "alphanum"))
            } else {
                sprintf("analysis=%s", sanitize1(analysis, type = "count")) # note no quotes
            }
        },
        
        if (missing(description_contains)) NULL
        else sprintf("description ILIKE '%%%s%%'", sanitize(description_contains, type = "text")), # Sanitation may be too restrictive, should do something intelligent with whitespace

        if (missing(phenotype_contains)) NULL
        else sprintf("phenotype ILIKE '%%%s%%'", sanitize(phenotype_contains, type = "text")), # Sanitation may be too restrictive

        if (missing(has_tag)) NULL
        else sprintf("tag='%s'", sanitize(has_tag, type = "alphanum")) # Sanitation may be too restrictive
    )
    ## format
    ws1f <- paste0("(", 
                  unlist(sapply(ws1, function(x) if (is.null(x)) NULL else paste0(tablename, x, collapse = " OR "))), 
                  ")", collapse = " OR ")

    ##
    ## analysis_not, ncase_ge, ncohort_ge etc. are AND'ed within and between arguments
    ws2 <- list(
        if (missing(analysis_not)) NULL
        else {
            if (getOption('gtx.analysisIsString', TRUE)) { # Future release will change default to FALSE
                sprintf("analysis!='%s'", sanitize(analysis_not, type = "alphanum"))
            } else {
                sprintf("analysis!=%s", sanitize1(analysis_not, type = "count")) # note no quotes
            }
        },
        
        if (missing(ncase_ge)) NULL
        else sprintf("ncase >= %s", sanitize(ncase_ge, type = "int")),
        
        if (missing(ncohort_ge)) NULL
        else sprintf("ncohort >= %s", sanitize(ncohort_ge, type = "int"))
    )
    ws2f <- paste0("(", 
                  unlist(sapply(ws2, function(x) if (is.null(x)) NULL else paste0(tablename, x, collapse = " AND "))), 
                  ")", collapse = " AND ")

    ## FIXME bad hack make this better
    if (ws1f == '()' && ws2f == '()') {
        stop('gtxwhat() produced no conditions')
    } else if (ws1f != '()' && ws2f == '()') {
        return(paste0('(', ws1f, ')'))
    } else if (ws1f == '()' && ws2f != '()') {
        return(paste0('(', ws2f, ')'))
    } else if (ws1f != '()' && ws2f != '()') {
        return(paste0('((', ws1f, ') AND (', ws2f, '))'))
    } else {
        stop('gtxwhat() logical error')
    }
}

gtxfilter <- function(pval_le, pval_gt,
                      maf_ge, maf_lt,
                      rsq_ge, rsq_lt, 
                      analysis,
                      tablename,
                      dbc = getOption("gtx.dbConnection", NULL)) {
    ## function to construct a WHERE string for constructing SQL queries
    ## logical AND within and between arguments
    ## queries TABLE analyses to determine results_db table name and
    ## to efficiently handle whether freq and rsq are NULL

    ## query analysis metadata, using cache if it exists
    amd <- sqlWrapper(getOption('gtx.dbConnection_cache_analyses', dbc),
                      sprintf('SELECT results_db, has_freq, has_rsq
                                    FROM analyses
                                    WHERE %s;',
                            gtxwhat(analysis1 = analysis)))
    
    if (!missing(tablename)) {
        tablename <- paste0(sanitize1(tablename, type = 'alphanum'), '.')
    } else {
        tablename <- ''
    }

    ## set flags for whether to print warning messages
    ## can be set at multiple places in code, but we only want to print messages once
    warn_freq <- FALSE
    warn_rsq <- FALSE

    ## unlike similar functions, gtxwhat(), gtxwhere(), that add the tablename at the
    ## end, for this function the specific OR logic needed for maf_lt means we have to
    ## paste the tablename within every clause (grrr)

    ws1 <- list(
        if (missing(pval_le)) NULL
        else sprintf('%spval<=%s', tablename, sanitize1(pval_le, type = 'double')),

        if (missing(pval_gt)) NULL
        else sprintf('%spval>%s', tablename, sanitize1(pval_gt, type = 'double')),

        if (missing(maf_ge)) NULL # Two list elements for freq>= and freq<= if maf_ge argument used
        else {
            if (amd$has_freq == 1) { # type BOOLEAN in Hive/Impala is returned to R as 0/1
                sprintf('(%sfreq>=%s AND %sfreq<=%s)',
                        tablename, sanitize1(maf_ge, type = 'double'),
                        tablename, sanitize1(1 - maf_ge, type = 'double'))
            } else {
                warn_freq <- TRUE
                NULL
            }
        },

        ## Note, maf_lt needs OR logic, where everything else uses AND
        if (missing(maf_lt)) NULL # Two list elements for freq< and freq> if maf_lt argument used
        else {
            if (amd$has_freq == 1) { # type BOOLEAN in Hive/Impala is returned to R as 0/1
                sprintf('(%sfreq<%s OR %sfreq>%s)',
                        tablename, sanitize1(maf_lt, type = 'double'),
                        tablename, sanitize1(1 - maf_lt, type = 'double'))
            } else {
                warn_freq <- TRUE
                NULL
            }
        },

        if (missing(rsq_ge)) NULL
        else {
            if (amd$has_rsq == 1) { # type BOOLEAN in Hive/Impala is returned to R as 0/1
                sprintf('%srsq>=%s', tablename, sanitize1(rsq_ge, type = 'double'))
            } else {
                warn_rsq <- TRUE
                NULL
            }
        },

        if (missing(rsq_lt)) NULL
        else {
            if (amd$has_rsq == 1) { # type BOOLEAN in Hive/Impala is returned to R as 0/1
                sprintf('%srsq<%s', tablename, sanitize1(rsq_lt, type = 'double'))
            } else {
                warn_rsq <- TRUE
                NULL
            }
        }
    )
    if (warn_freq) warning('gtxfilter MAF filtering skipped because has_freq=False for analysis [ ', analysis, ' ]')
    if (warn_rsq) warning('gtxfilter Rsq filtering skipped because has_rsq=False for analysis [ ', analysis, ' ]')

    ## format, can simplify the inner paste0 since all list elements have length 1
    ws1f <- paste0("(", 
                   unlist(ws1),
#                   unlist(sapply(ws1, function(x) if (is.null(x)) NULL else paste0(tablename, x))), 
                   ")", collapse = " AND ")
    ## handle case where no arguments were non-missing but still need to generate
    ## valid SQL clause for substitution into '... WHERE %s' or '... WHERE %s AND %s'
    if (ws1f == '()') ws1f <- '(True)'
    return(ws1f)
}

gtxfilter_label <- function(maf_ge, maf_lt,
                            rsq_ge, rsq_lt, 
                            analysis,
                            dbc = getOption("gtx.dbConnection", NULL)) {
    ## function to construct a nice label for plot

    ## query analysis metadata, note this is inefficient since gtxfilter() may
    ## be making the same/similar query (potentially multiple times)
    ## added, using cache if it exists
    amd <- sqlWrapper(getOption('gtx.dbConnection_cache_analyses', dbc),
                      sprintf('SELECT has_freq, has_rsq
                                    FROM analyses
                                    WHERE %s;',
                            gtxwhat(analysis1 = analysis)))
    
    ## set flags for whether to print warning messages
    ## can be set at multiple places in code, but we only want to print messages once
    warn_freq <- FALSE
    warn_rsq <- FALSE

    ws1 <- list(
        if (missing(maf_ge)) NULL
        else {
            if (amd$has_freq == 1) { # type BOOLEAN in Hive/Impala is returned to R as 0/1
                sprintf('MAF>=%s', sanitize1(maf_ge, type = 'double'))
            } else {
                warn_freq <- TRUE
                NULL
            }
        },

        if (missing(maf_lt)) NULL
        else {
            if (amd$has_freq == 1) { # type BOOLEAN in Hive/Impala is returned to R as 0/1
                sprintf('MAF<%s', sanitize1(maf_lt, type = 'double'))
            } else {
                warn_freq <- TRUE
                NULL
            }
        },

        if (missing(rsq_ge)) NULL
        else {
            if (amd$has_rsq == 1) { # type BOOLEAN in Hive/Impala is returned to R as 0/1
                sprintf('rsq>=%s', sanitize1(rsq_ge, type = 'double'))
            } else {
                warn_rsq <- TRUE
                NULL
            }
        },

        if (missing(rsq_lt)) NULL
        else {
            if (amd$has_rsq == 1) { # type BOOLEAN in Hive/Impala is returned to R as 0/1
                sprintf('rsq<%s', sanitize1(rsq_lt, type = 'double'))
            } else {
                warn_rsq <- TRUE
                NULL
            }
        }
    )

    ## FIXME make nicer formatting like 0.01>MAF>=0.001
    
    ws1f <- paste0(unlist(ws1), collapse = " and ")
    if (ws1f == '') {
        ws1f <- 'All variants'
    } else {
        ws1f <- paste('Variants with', ws1f)
    }
    return(ws1f)

    ## FIXME add text if filtering requested but could not be applied
    ##    if (warn_freq) warning('gtxfilter MAF filtering skipped because has_freq=False for analysis [ ', analysis, ' ]')
    ##    if (warn_rsq) warning('gtxfilter Rsq filtering skipped because has_rsq=False for analysis [ ', analysis, ' ]')


}

## function to return pretty printing label for an analysis
## analysis is the analysis id
## entity is the result of a call to gtxentity (i.e. a list with elements entity, entity_label)
gtxanalysis_label <- function(analysis, entity, nlabel = TRUE,
                              dbc = getOption("gtx.dbConnection", NULL)) {
    ares <- sqlWrapper(dbc,
                       sprintf('SELECT label, ncase, ncontrol, ncohort FROM analyses WHERE %s;',
                               gtxwhat(analysis1 = analysis)))
    if (nlabel) {
        if (!is.na(ares$ncase) && !is.na(ares$ncontrol)) {
            alabel <- sprintf('%s, n=%i vs %i', ares$label, ares$ncase, ares$ncontrol)
        } else if (!is.na(ares$ncohort)) {
            alabel <- sprintf('%s, n=%i', ares$label, ares$ncohort)
        } else {
            alabel <- sprintf('%s, n=?', ares$label)
        }
    } else {
        alabel <- ares$label
    }
    if (!is.null(entity)) alabel <- paste(entity$entity_label, alabel)
    return(alabel)
}

gtxanalyses <- function(analysis, analysis_not, 
                        phenotype_contains,
                        description_contains,
                        has_tag,
                        ncase_ge,
                        ncohort_ge,
                        ## if extra filters are added, be sure to update definition of all_analyses below
                        analysis_fields = c('description', 'phenotype', 'covariates', 'cohort', 'unit',
                                    'ncase', 'ncontrol', 'ncohort'),
                        tag_is, with_tags = FALSE,
                        has_access_only = FALSE, 
                        dbc = getOption("gtx.dbConnection", NULL)) {
    gtxdbcheck(dbc)
    dbs <- sqlWrapper(dbc, 'SHOW DATABASES;', uniq = FALSE)

    ## sanitize and check desired analysis_fields, create sanitized string for SQL
    ## added, using cache if it exists
    acols <- names(sqlWrapper(getOption('gtx.dbConnection_cache_analyses', dbc),
                              'SELECT * FROM analyses LIMIT 0', zrok = TRUE)) # columns present in TABLE analyses
    analysis_fields <- sanitize(union(c('analysis', 'entity_type', 'results_db'),
                              analysis_fields),
                        type = 'alphanum')
    analysis_fields_bad <- !(analysis_fields %in% acols)
    if (any(analysis_fields_bad)) {
        error('analysis_fields [ ', paste(analysis_fields[analysis_fields_bad], collapse = ', '), ' ] not found in TABLE analyses')
    }
    rm(acols) # clean up
    analysis_fields <- paste0('analyses.', analysis_fields, collapse = ', ')

    if (!missing(has_tag)) with_tags <- TRUE
    
    ## sanitize and check desired tag_is logicals
    if (!missing(tag_is)) {
        acols <- names(sqlWrapper(dbc, 'SELECT * FROM analyses_tags LIMIT 0', zrok = TRUE)) # columns present in TABLE analyses
        tag_is <- sanitize(tag_is, type = 'alphanum')
        tag_is_bad <- !(tag_is %in% acols)
        if (any(tag_is_bad)) {
            error('tag_is logicals [ ', paste(tag_is[tag_is_bad], collapse = ', '), ' ] not found in TABLE analyses_tags')
        }
        rm(acols)
        paste0('(analysis_tags.', tag_is, ')', collapse = ' AND ')
        with_tags <- TRUE
    } else {
        tag_is <- '(True)'
    }
    
    all_analyses <- (missing(analysis) && missing(analysis_not) && missing(phenotype_contains) &&
                     missing(description_contains) && missing(ncase_ge) && missing(ncohort_ge) &&
                     missing(has_tag))
    if (all_analyses) {
        ## This query cannot use cache unless JOIN with analyses_tags, and subquery, are supported
        res <- sqlWrapper(dbc,
                          sprintf('SELECT %s %s FROM analyses %s',
                                  analysis_fields,
                                  if (with_tags) ', t1.tag' else '',
                                  if (with_tags) sprintf('LEFT JOIN (SELECT analysis, tag FROM analyses_tags WHERE %s) AS t1 USING (analysis)', tag_is) else ''),
                                  ## This LEFT JOIN does not work as expected when tag_is is used
                                  ## because, logically, the left join happens BEFORE the WHERE clause is applied
                                  ## need the tag_is as a subquery because ... EXPLAIN
                          uniq = FALSE)
    } else {
        res <- sqlWrapper(dbc,
                          sprintf('SELECT %s %s FROM analyses %s WHERE %s AND %s',
                                  analysis_fields,
                                  if (with_tags) ', analyses_tags.tag' else '',
                                  if (with_tags) 'LEFT JOIN analyses_tags USING (analysis)' else '',
                                  gtxwhat(analysis = analysis, analysis_not = analysis_not, 
                                          description_contains = description_contains,
                                          phenotype_contains = phenotype_contains,
                                          has_tag = has_tag, 
                                          ncase_ge = ncase_ge, ncohort_ge = ncohort_ge),
                                  tag_is),
                          uniq = FALSE)
    }
    res$has_access <- res$results_db %in% dbs$name
    if (has_access_only) {
        res <- subset(res, has_access)
    }
    return(res)
}

gtxregion <- function(chrom, pos_start, pos_end, 
                      hgncid, ensemblid, pos, rs, surround = 500000, 
                      dbc = getOption("gtx.dbConnection", NULL)) {
  gtxdbcheck(dbc)
  if (!missing(chrom) && !missing(pos_start) && !missing(pos_end)) {
    stopifnot(identical(length(chrom), 1L))
    stopifnot(identical(length(pos_start), 1L))
    stopifnot(identical(length(pos_end), 1L))
    stopifnot(pos_end > pos_start)
  } else if (!missing(hgncid)) {
    gp <- sqlWrapper(dbc,
                     sprintf('SELECT chrom, min(pos_start) as pos_start, max(pos_end) as pos_end FROM genes WHERE %s GROUP BY chrom',
                             gtxwhere(hgncid = hgncid)),
                     uniq = FALSE) # require unique result but catch below to provide more informative error
    if (!identical(nrow(gp), 1L)) stop('hgncid(s) [ ', paste(hgncid, collapse = ', '), ' ] do not map to unique chromosome')
    chrom <- gp$chrom[1] # do we need the [1] here?
    pos_start <- gp$pos_start[1] - surround
    pos_end <- gp$pos_end[1] + surround
  } else if (!missing(ensemblid)) {
    gp <- sqlWrapper(dbc,
                     sprintf('SELECT chrom, min(pos_start) as pos_start, max(pos_end) as pos_end FROM genes WHERE %s GROUP BY chrom',
                             gtxwhere(ensemblid = ensemblid)),
                     uniq = FALSE) # require unique result but catch below to provide more informative error
    if (!identical(nrow(gp), 1L)) stop('ensemblid(s) [ ', paste(ensemblid, collapse = ', '), ' ] do not map to unique chromosome')
    chrom <- gp$chrom[1]
    pos_start <- gp$pos_start[1] - surround
    pos_end <- gp$pos_end[1] + surround
    ## Note pos and rs low in priority order to allow secondary use in regionplot()
  } else if (!missing(chrom) && !missing(pos)) {
    stopifnot(identical(length(chrom), 1L))
    stopifnot(identical(length(pos), 1L))
    stopifnot(identical(length(surround), 1L))
    pos_start <- pos - surround
    pos_end <- pos + surround
  } else if (!missing(rs)) {
    gp <- sqlWrapper(dbc,
                     sprintf('SELECT chrom, pos FROM sites WHERE %s;',
                             gtxwhere(rs = rs))) # default uniq = TRUE
    chrom <- gp$chrom
    pos_start <- gp$pos - surround
    pos_end <- gp$pos + surround
  } else {
    stop('gtxregion() failed due to inadequate arguments')
  }

  return(list(chrom = chrom, pos_start = pos_start, pos_end = pos_end,
              label =   sprintf('chr%s:%s-%s', chrom, pos_start, pos_end)))
}


## infer the entity id according to the type required for the analysis
gtxentity <- function(analysis, entity, hgncid, ensemblid, 
                      dbc = getOption("gtx.dbConnection", NULL)) {
    gtxdbcheck(dbc)

    entity_type <- sqlWrapper(dbc,
                              sprintf('SELECT entity_type FROM analyses WHERE %s;',
                                      gtxwhat(analysis1 = analysis))
                              )$entity_type
    ## FIXME this logging message could be improved when entity_type is null or empty string
    gtxlog('analysis [ ', analysis, ' ] requires entity type [ ', entity_type, ' ]')
    
    if (identical(entity_type, 'ensg')) {
        ## analysis requires an entity that is Ensembl gene id
        if (!missing(entity)) {
            ## entity argument supplied, infer whether Ensembl gene id or HGNC id
            entity <- sanitize1(entity, type = 'alphanum')
            if (grepl('^ENSG[0-9]+$', entity)) {
                gtxlog('Recognized entity [ ', entity, ' ] as Ensembl gene id')
                ## do nothing
            } else {
                ## attempt to convert assumed HGNC id to Ensembl gene id
                gtxlog('Looking for Ensembl gene id for entity [ ', entity, ' ]')
                entity <- sqlWrapper(dbc,
                                     sprintf('SELECT ensemblid FROM genes WHERE %s;', 
                                             gtxwhere(hgncid = entity))
                                     )$ensemblid
                gtxlog('Converted entity to Ensembl gene id [ ', entity, ' ]')
            }
        } else {
            ## infer entity from other arguments if supplied
            if (!missing(hgncid)) {
                entity <- sqlWrapper(dbc,
                                     sprintf('SELECT ensemblid FROM genes WHERE %s;', 
                                             gtxwhere(hgncid = hgncid))
                                     )$ensemblid
                gtxlog('Inferred entity [ ', entity, ' ] from HGNC id [ ', hgncid, ' ]')
            } else if (!missing(ensemblid)) {
                entity <- sanitize1(ensemblid, type = 'ENSG')
                gtxlog('Inferred entity [ ', entity, ' ] from Ensembl gene id [ ', ensemblid, ' ]')
            } else {
                stop('Analysis [ ', analysis, ' ] requires an Ensembl gene id entity but none could be inferred')
            }
        }
        entity_label <- sqlWrapper(dbc,
                                   sprintf('SELECT hgncid FROM genes WHERE %s;',
                                           gtxwhere(ensemblid = entity))
                                   )$hgncid
        if (is.na(entity_label) || entity_label == '') entity_label <- entity # may have to paste 'ENSG' back on when switch to ints
        return(list(entity = entity, entity_label = entity_label))
    } else {
        return(NULL)
        # entity_type is NULL or unrecognized type
    }
    stop('internal error in gtxentity')
}

# for certain types, sanitize and gtxlabel are (almost) inverses
gtxlabel <- function(x, type) {
    if (identical(type, 'ensg') || identical(type, 'rs')) { ## could nest up other ENS[PGT] types...
        x <- na.omit(x) ## silently drop missing values
        xi <- suppressWarnings(as.integer(x))
        if (any(is.na(xi) | xi < 0)) {
            stop('gtxlabel invalid input [ ', paste(x[is.na(xi) | xi < 0], collapse = ', '),
                 ' ] for type [ ', type, ' ]')
        }
        if (identical(type, 'ensg')) {
            return(sprintf('ENSG%012i', x))
        } else if (identical(type, 'rs')) {
            return(sprintf('rs%i', x))
        } else {
            stop('internal error in gtxlabel()')
        }
    } else {
        stop('invalid type [ ', type, ' ] in gtxlabel()')
    }
}

