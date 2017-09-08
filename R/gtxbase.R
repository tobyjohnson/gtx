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
        else sprintf("chrom='%s'", sanitize(chrom[1], values = c(as.character(1:22), "X", "Y"))),
        ## FIXME: Should this be sanitize1() ?
        
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
        else sprintf("rs='%s'", sanitize(rs, type = "rs")),
        
        if (missing(hgncid)) NULL
        else sprintf("hgncid='%s'", sanitize(hgncid, type = "alphanum-")), # thousands of HGNC ids contain hyphens, e.g. HLA-A
        
        if (missing(ensemblid)) NULL
        else sprintf("ensemblid='%s'", sanitize(ensemblid, type = "ENSG"))
    )
    ws2 <- paste0("(", 
                  unlist(sapply(ws1, function(x) if (is.null(x)) NULL else paste0(tablename, x, collapse = " OR "))), 
                  ")", collapse = " AND ")
    return(ws2)
}

##
## convenience function to construct WHERE
## part of SQL for analyses table
## Note behaviour here is OR/OR, different to gtxwhere()
##
gtxwhat <- function(analysis,
                    description_like,
                    phenotype_like,
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
    
    ws1 <- list(
        if (missing(analysis)) NULL
        else sprintf("analysis='%s'", sanitize(analysis, type = "alphanum"))
        
        if (missing(description_like)) NULL
        else sprintf("description ILIKE '%s'", sanitize(description_like, type = "alphanum")), # Too restrictive

        if (missing(phenotype_like)) NULL
        else sprintf("phenotype ILIKE '%s'", sanitize(phenotype_like, type = "alphanum")) # Too restrictive
    )
    ws2 <- paste0("(", 
                  unlist(sapply(ws1, function(x) if (is.null(x)) NULL else paste0(tablename, x, collapse = " OR "))), 
                  ")", collapse = " AND ")
    return(ws2)
}

gtxanalyses <- function(analysis,
                        phenotype_like,
                        description_like,
                        dbc = getOption("gtx.dbConnection", NULL)) {
    gtxdbcheck(dbc)

    return(sqlWrapper(dbc, sprintf('SELECT analysis, description, phenotype, ncase, ncontrol, ncohort FROM analyses %s',
                                   gtxwhat(analysis = analysis,
                                           description_like = description_like,
                                           phenotype_like = phenotype_like)),
                      uniq = FALSE))
}

gtxregion <- function(chrom, pos_start, pos_end,
                           hgncid, ensemblid, surround = 500000, 
                           dbc = getOption("gtx.dbConnection", NULL)) {
  gtxdbcheck(dbc)
  if (!missing(chrom) & !missing(pos_start) & !missing(pos_end)) {
    stopifnot(identical(length(chrom), 1L))
    stopifnot(identical(length(pos_start), 1L))
    stopifnot(identical(length(pos_end), 1L))
    stopifnot(pos_end > pos_start)
  } else if (!missing(hgncid)) {
    ## verify dbc is a database connection
    gp <- sqlQuery(dbc, sprintf('SELECT chrom, min(pos_start) as pos_start, max(pos_end) as pos_end FROM genes WHERE %s GROUP BY chrom',
                                gtxwhere(hgncid = hgncid)))
    if (!identical(nrow(gp), 1L)) stop('hgncid(s) [ ', paste(hgncid, collapse = ', '), ' ] do not map to unique chromosome')
    chrom <- gp$chrom[1]
    pos_start <- gp$pos_start[1] - surround
    pos_end <- gp$pos_end[1] + surround
  } else if (!missing(ensemblid)) {
    ## verify dbc is a database connection
    gp <- sqlQuery(dbc, sprintf('SELECT chrom, min(pos_start) as pos_start, max(pos_end) as pos_end FROM genes WHERE %s GROUP BY chrom',
                                gtxwhere(ensemblid = ensemblid)))
    if (!identical(nrow(gp), 1L)) stop('ensemblid(s) [ ', paste(ensemblid, collapse = ', '), ' ] do not map to unique chromosome')
    chrom <- gp$chrom[1]
    pos_start <- gp$pos_start[1] - surround
    pos_end <- gp$pos_end[1] + surround
  }
  return(list(chrom = chrom, pos_start = pos_start, pos_end = pos_end))
}


## infer the entity id according to the type required for the analysis
gtxentity <- function(analysis, entity, hgncid, ensemblid, 
                      dbc = getOption("gtx.dbConnection", NULL)) {
    gtxdbcheck(dbc)

    entity_type <- sqlWrapper(dbc,
                              sprintf('SELECT entity_type FROM analyses WHERE analysis=\'%s\';',
                                      sanitize1(analysis, type = 'alphanum'))
                              )$entity_type
    gtxlog('analysis [ ', analysis, ' ] requires entity type [ ', entity_type, ' ]')
    
    if (identical(entity_type, factor('ensg'))) {
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
                                     sprintf('SELECT ensemblid FROM genes WHERE hgncid=\'%s\';', 
                                             entity)
                                     )$ensemblid
                gtxlog('Converted entity to Ensembl gene id [ ', entity, ' ]')
            }
        } else {
            ## infer entity from other arguments if supplied
            if (!missing(hgncid)) {
                entity <- sqlWrapper(dbc,
                                     sprintf('SELECT ensemblid FROM genes WHERE hgncid=\'%s\';', 
                                             sanitize(hgncid, type = 'hgnc'))
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
                                   sprintf('SELECT hgncid FROM genes WHERE ensemblid=\'%s\';',
                                           sanitize1(entity, type = 'ENSG'))
                                   )$hgncid
        if (is.na(entity_label) || entity_label == '') entity_label <- entity # may have to paste 'ENSG' back on when switch to ints
        return(list(entity = entity, entity_label = entity_label))
    } else {
        return(NULL)
        # entity_type is NULL or unrecognized type
    }
    stop('internal error in gtxentity')
}




