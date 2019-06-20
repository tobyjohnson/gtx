#define s4 class for gwas
setClass(
  'regionplot',
  representation(v1Data = 'data.table',
                 phewasArguments = 'list'),
  prototype(v1Data = NULL,
            phewasArguments = NULL)
)

#validate phewas user defined inputs
checkPhewasInput <- function(chrom,
                             pos,
                             ref,
                             alt,
                             rs,
                             analysis,
                             analysis_not,
                             description_contains,
                             phenotype_contains,
                             has_tag,
                             ncase_ge,
                             ncohort_ge,
                             analysis_fields,
                             tag_is,
                             plot_ymax,
                             tags_to_right,
                             dbc,
                             connectionType,
                             phewasObj = NULL,
                             callFromAPI) {
  #dbc
  gtxdbcheck(dbc)
  
  #connectionType
  checkType(connectionType)
  
  #phewasObj
  if (!is.null(phewasObj)) {
    if (class(phewasObj) != 'phewas') {
      stop('phewasObj must be an object of class phewas')
    }
  }
  
  #callFromAPI
  if (!is.logical(callFromAPI)) {
    stop('callFraomAPI must be logical')
  }
  
}


getPhewasDataArgs <- function(phewasArgs) {
  ###backwards way to call phewas.data given that some arguments might be missing
  #get names of args in phewas.data
  phewasDataArgsNames <- names(as.list(args(phewas.data)))
  phewasDataArgsNames <-
    phewasDataArgsNames[phewasDataArgsNames != '']
  phewasDataArgs <-
    phewasArgs[names(phewasArgs) %in% phewasDataArgsNames]
  phewasDataArgs$with_tags <- TRUE
  
  return(phewasDataArgs)
}


#' Phenome Wide Assocation Study
#'
#' @description  Look up and plot results for a Phenome Wide Association Study.
#'
#' @details Look up association summary statistics for a single variant of interest, specified
#' either by chromosome and position, or by dbSNP rs identifier.  (FIXME need to add optional
#' ref/alt arguments to disambiguate when there are multiple variants at the same position,
#' or multiple variants with the same rs.)
#'
#' By default, summary statistics are returned for all analyses where the variant of interest
#' was tested.  If results are required for only a specific set of analyses, these can be
#' specified using a combination of arguments that are passed through to \code{gtxanalyses()}.
#'
#' @param chrom: UNDEFINED
#' @param pos: Position of variant of interest
#' @param ref: Reference allele of variant of interest
#' @param alt: Alternate allele of variant of interest
#' @param rs: dbSNP rs identifier of variant of interest
#' @param analysis: Identifiers of analyses to include
#' @param analysis_not: Identifiers of analyses to exclude
#' @param description_contains: Search term for analyses to include
#' @param phenotype_contains: UNDEFINED
#' @param has_tag: Tag value for analyses to include
#' @param ncase_ge: Threshold case sample size greater-or-equal for inclusion
#' @param ncohort_ge: Threshold cohort sample size greater-or-equal for inclusion
#' @param analysis_fields: Fields from analysis metadata to return
#' @param tag_is: Tag flag to use for plot grouping
#' @param nearby: Additionally return minimum pvalue at position+-nearby
#' @param plot_ymax: Y-axis value maximum
#' @param tags_to_right: Tag values to force to right of plot, next to "No tag"
#' @param dbc: Database connection
#' @param connectionType (character) The kind of connection used to retrive
#' data. Only 'SQL' allowed
#' @param callFromAPI (logical) Deterimines if the function return an S4 object
#' (used by API), or data in legacy format (i.e. dataframe)
#'
#' @usage
#' phewas(chrom, pos, ref, alt, rs, analysis, analysis_not,
#' description_contains, phenotype_contains, has_tag, ncase_ge, ncohort_ge,
#' analysis_fields = c("description", "label", "unit", "ncase",
#' "ncontrol", "ncohort"), tag_is = "phewas_group", nearby,
#' plot_ymax = 30, tags_to_right = NULL,
#' dbc = getOption("gtx.dbConnection", NULL))
#'
#' @author Toby Johnson \email{Toby.x.Johnson@gsk.com}
#' @export
# Could inherit param documentation from gtxanalyses
# Should merge phewas.qq into main phewas as an alternative plot style
phewas <- function(chrom,
                   pos,
                   ref,
                   alt,
                   rs,
                   analysis,
                   analysis_not,
                   description_contains,
                   phenotype_contains,
                   has_tag,
                   ncase_ge,
                   ncohort_ge,
                   ## if extra filters are added, be sure to update definition of all_analyses below
                   analysis_fields = c('description', 'label', 'unit',
                                       'ncase', 'ncontrol', 'ncohort'),
                   tag_is = 'phewas_group',
                   nearby,
                   plot_ymax = 30,
                   tags_to_right = NULL,
                   dbc = getOption("gtx.dbConnection", NULL),
                   connectionType = 'SQL',
                   callFromAPI = FALSE) {
  #deal with optional arguments
  ## validate nearby argument
  if (!missing(nearby)) {
    nearby <- as.integer(nearby)
    if (nearby <= 0)
      nearby <- NULL
  }
  
  phewasObj <- new('phewas')
  
  #get list of calling args
  phewasArgs <- as.list(match.call())
  phewasArgs[[which(names(phewasArgs) == '')]] <- NULL
  phewasObj@phewasArguments <- phewasArgs
  
  #checkInputs
  do.Call(checkPhewasInput, phewasArgs)
  
  phewasDataArgs <- getPhewasDataArgs(phewasArgs)
  phewasObj@phewasDataArgs <- phewasDataArgs
  #we add the object now to avoid storing multiple copis of the objetc within itself
  phewasDataArgs$phewasObj = phewasObj
  
  
  
  #get the data
  p1 <- do.call(phewas.data, phewasDataArgs)
  
  if (nrow(p1) == 0L) {
    gtx_warn('PheWAS had zero results, skipping plotting')
    return(invisible(p1))
  }
  
  p1orig <- p1
  
  ymax <- max(10, ceiling(-log10(min(p1$pval))))
  if (ymax > plot_ymax) {
    # Hard coded threshold makes sense for control of visual display
    ymax <- plot_ymax
    gtx_warn('Truncating P-values at 1e-{ymax}')
    truncp <- TRUE
  } else {
    truncp <- FALSE
  }
  
  p1 <-
    p1[order(p1$tag, p1$pval), ] # with this sort, tag==NA will be last
  p1$tag[is.na(p1$tag)] <- 'No tag'
  p1$x1 <-
    1:nrow(p1) # initial attempttag_is at x axis position, to be refined below
  x1 <- with(
    aggregate(
      p1$x1,
      by = list(tag = p1$tag),
      FUN = function(x)
        return(c(min = min(x), max = max(x)))
    ),
    cbind(data.frame(tag), as.data.frame(x))
  )
  # reorder because 'No tag' will be sorted within x1:
  x1 <-
    x1[order(match(x1$tag, c(tags_to_right, 'No tag')), x1$tag, na.last = FALSE), ]
  group_interspace <-
    50 # this should be a user configurable parameter
  x1$offset <-
    c(0, cumsum(x1$max - x1$min + group_interspace))[1:nrow(x1)] - x1$min + group_interspace
  x1$midpt <- 0.5 * (x1$max + x1$min) + x1$offset
  if (x1$tag[nrow(x1)] == 'No tag') {
    tmp_col <-
      col2rgb(c(rainbow(nrow(x1) - 1), 'grey')) # rainbow for all tags, grey for last
  } else {
    tmp_col <- col2rgb(rainbow(nrow(x1))) # rainbow for all tags
  }
  x1$col_red <- tmp_col['red', , drop = TRUE]
  x1$col_green <- tmp_col['green', , drop = TRUE]
  x1$col_blue <- tmp_col['blue', , drop = TRUE]
  p1m <- match(p1$tag, x1$tag)
  p1$x <- p1$x1 + x1$offset[p1m]
  if (!missing(nearby) && !is.null(nearby)) {
    # colour scaling so that:
    # strength = 1      when pval == pval_nearby
    # strength = 0.61   when pval == pval_nearby * 10 (one log10 less significant)
    # strength = 0.37   when pval == pval_nearby * 100
    # strength = 0.22   when pval == pval_nearby * 1000
    # etc.
    p1$strength <- (p1$pval_nearby / p1$pval) ^ (0.5 / log(10))
    w0 <- which(p1$pval == 0.)
    if (length(w0) > 0) {
      gtx_warn('PheWAS results include {length(w0)} analyses with pval=0.; assuming strength=1')
      p1$strength[w0] <- 1.
    }
    wna <-
      which(!is.finite(p1$strength) |
              p1$strength < 0. | p1$strength > 1.)
    if (length(wna) > 0) {
      gtx_warn(
        'PheWAS results include {length(wna)} analyses with invalid strength; assuming strength=0'
      )
      p1$strength[wna] <- 0.
    }
  } else {
    p1$strength <- 1.
  }
  # for point background, interpolate RGB space linearly,
  # (strength) of the way from main (point outline) colour to white
  p1$col <-
    rgb(x1$col_red[p1m], x1$col_green[p1m], x1$col_blue[p1m], maxColorValue = 255)
  p1$col_bg <-
    rgb((p1$strength * x1$col_red[p1m] + (1. - p1$strength) * 255),
        (p1$strength * x1$col_green[p1m] + (1. - p1$strength) * 255),
        (p1$strength * x1$col_blue[p1m] + (1. - p1$strength) * 255),
        maxColorValue = 255
    )
  # note still sorted by tag then pval
  p1$y <- pmin(-log10(p1$pval), ymax) # guarantee all y <= ymax
  p1$yo <-
    c(Inf, ifelse(p1$tag[-1] != p1$tag[-nrow(p1)], Inf, p1$y[-nrow(p1)] - p1$y[-1])) # yaxis space within group relative to previous point
  
  p1 <-
    p1[order(p1$pval, decreasing = TRUE), ] # reorder to make sure most significant points plotted over less significant ones
  
  plot.new()
  plot.window(range(p1$x), c(0, 1.1 * ymax))
  abline(h = 0, col = "grey")
  
  if (!missing(nearby) && !is.null(nearby)) {
    ### Ad hoc legend for colour scheme
    nearby_legend <-
      c(
        paste0('association within ', round(nearby * 1.e-3), 'kb'),
        'equal',
        '10x',
        '100x',
        '1000x more significant'
      )
    xoff <-
      cumsum(
        strwidth(nearby_legend, units = 'user', cex = 0.5) +
          3. * strwidth('M', units = 'user', cex = 0.5)
      )
    # FIXME dynamically shrink cex if bigger than par('usr')[1:2]
    yoff <- 0.5 * (ymax + par('usr')[4])
    points(
      xoff[1:4],
      rep(yoff, 4),
      pch = 22,
      cex = 1,
      col = rgb(1, 0, 0),
      bg = rgb(1, 1 - exp(-c(
        0.0, 0.5, 1.0, 1.5
      )), 1 - exp(-c(
        0.0, 0.5, 1.0, 1.5
      )))
    )
    text(
      c(0, xoff[1:4]),
      rep(yoff, 4),
      pos = 4,
      cex = 0.5,
      labels = nearby_legend
    )
  }
  max_y = max(p1$y)
  ## plot text then overplot with points
  with(subset(p1, y > min(6, max_y * .96) &
                yo > strheight('M', cex = 0.5)),
       # should be user controllable
       text(x, y, as.character(label), pos = 4, cex = 0.5))
  points(
    p1$x,
    p1$y,
    pch = ifelse(p1$beta > 0., 24, 25),
    col = p1$col,
    bg = p1$col_bg,
    cex = 0.75
  )
  if (truncp) {
    ## would be nice to more cleanly overwrite y axis label
    axis(2,
         at = ymax,
         labels = substitute({
           
         } >= ymax, list(ymax = ymax)),
         las = 1)
    # could e.g. do setdiff(pretty(...), ymax)
  }
  ypretty <- pretty(c(0, ymax))
  axis(2, at = ypretty[ypretty <= ymax], las = 1)
  mtext(
    x1$tag,
    side = 1,
    line = 0.5,
    las = 2,
    at = x1$midpt,
    col = x1$col,
    las = 2,
    cex = 0.5
  )
  mtext(expression(-log[10](paste(
    italic(P), "-value"
  ))), 2, 3)
  mtext.fit(main = paste('PheWAS for', attr(p1, 'variant')))
  box()
  
  return(invisible(p1orig))
}




#' @description phewas Look up results without plotting
#'
#' @inheritParams phewas
#' @param phewasObj (phewas) If not null, it will be used to get arguments
#' values from the original phewas() call. Users will normally not need to
#' concern with this argument, as it will be used by the call from the API
#'
#' @usage
#' #' phewas.data(chrom, pos, ref, alt, rs, analysis, analysis_not,
#' description_contains, phenotype_contains, has_tag, ncase_ge, ncohort_ge,
#' analysis_fields = c("description", "label", "unit", "ncase",
#' "ncontrol", "ncohort"), tag_is, with_tags = FALSE, nearby,
#' dbc = getOption("gtx.dbConnection", NULL))
#'
#' @return Data frame of summary statistics
#' @export
phewas.data <- function(chrom,
                        pos,
                        ref,
                        alt,
                        rs,
                        analysis,
                        analysis_not,
                        description_contains,
                        phenotype_contains,
                        has_tag,
                        ncase_ge,
                        ncohort_ge,
                        ## if extra filters are added, be sure to update
                        #definition of all_analyses below
                        analysis_fields = c('description', 'label', 'unit',
                                            'ncase', 'ncontrol', 'ncohort'),
                        tag_is,
                        with_tags = FALSE,
                        # not clear why missing with_tags
                        #cannot be passed through gtxanalyses()
                        nearby = NULL,
                        dbc = getOption("gtx.dbConnection", NULL),
                        connectionType = 'SQL',
                        phewasObj = NULL) {
  ## validate nearby argument
  if (!is.null(nearby)) {
    nearby <- as.integer(nearby)
    if (nearby <= 0) {
      nearby <- NULL
    }
  }
  
  #get calling arguments
  #get list of calling args
  phewasDataArgs <- as.list(match.call())
  phewasDataArgs[[which(names(phewasDataArgs) == '')]] <- NULL
  
  if (is.null(phewasObj)) {
    phewasObj <- new('phewas')
    phewasObj@phewasDataArgs <- phewasDataArgs
  }
  
  #checkInputs
  do.Call(checkPhewasInput, phewasDataArgs)
  
  ## Look up variant
  v1 <- switch(connectionType,
               'SQL' = phewasGetV1(chrom, pos, ref, alt, rs, dbc,
                                   uniq = FALSE, zrok = TRUE,
                                   connectionType = 'SQL'))
  
  phewasObj@v1Data <- v1
  
  
  phewasObj <- phewasDataProcessing(phewasObj)
  
  if (callFromAPI) {
    #return object
  } else {
    #return only df res
  }
  
  #was returning res
}


checkedV1 <- function(v1) {
  v1_label <- label_variant(v1$chrom, v1$pos, v1$ref, v1$alt, v1$rs)
  if (nrow(v1) == 0L) {
    gtx_warn('PheWAS query did not match any variant (check TABLE sites)')
    v1 <-
      NULL # go through rest of code to have single point where no results are handled
  } else if (nrow(v1) > 1L) {
    gtx_warn('PheWAS query matches variant to {v1_label}') # normally a debug message but raised to warning in this situation
    gtx_warn(
      'PheWAS query matched multiple variants, cannot provide results without disambiguating ref/alt arguments'
    )
    v1 <-
      NULL # go through rest of code to have single point where no results are handled
  } else {
    gtx_debug('PheWAS query matches variant to {v1_label}')
    
    return(v1)
  }
}

getAnalysisMetadata(phewasObj) {
  #get the appropriate arguments
  phewasDataArgs <- phewasObj@phewasDataArgs
  
  gtxanalysesArgsNames <- names(as.list(args(gtxanalyses)))
  gtxanalysesArgsNames <-
    gtxanalysesArgsNames[gtxanalysesArgsNames != '']
  gtxanalysesArgs <-
    phewasDataArgs[names(phewasDataArgs) %in% gtxanalysesArgsNames]
  gtxanalysesArgs$has_access_only <- TRUE
  
  a1 <-
    do.Call(gtxanalyses, gtxanalysesArgs) # will work fine if all
  #filtering arguments are missing, as it internally sets all_analyses<-TRUE
  gtx_debug('Queried metadata for {nrow(a1)} analyses')
  gtx_debug(
    '{sum(is.na(a1$results_db) | a1$results_db == "", na.rm = TRUE)} analyses in current database, {sum(a1$results_db != "", na.rm = TRUE)} analyses in other accessible databases'
  )
  
  phewasObj@a1 <- a1
  return(phewasObj)
}

phewasDataProcessing <- function(phewasObj) {
  ## Look up analysis metadata
  phewasObj <- getAnalysisMetadata(phewasObj)
  
  a1 <- phewasObj@a1
  
  
  
  ## Handle when results_db is NULL in database (returned as NA to R), since not using gtxanalysisdb()
  ## Add period if results_db is a database name, otherwise empty string (use pattern %sgwas_results in sprintf's below)
  ## (FIXME it would be better if this did go through gtxanalysisdb() for future maintenance)
  loop_clean_db <- sapply(unique(a1$results_db), function(x) {
    return(if (is.na(x) ||
               x == '')
      ''
      else
        sanitize1(paste0(x, '.'), type = 'alphanum.'))
  })
  
  phewasObj@loop_clean_db <- loop_clean_db
  
  all_analyses <-
    (
      missing(analysis) &&
        missing(analysis_not) && missing(phenotype_contains) &&
        missing(description_contains) &&
        missing(has_tag) &&
        missing(ncase_ge) && missing(ncohort_ge)
    )
  
  phewasObj@all_analyses <- all_analyses
  
  phewasObj <- getResPhewas(phewasObj)
  
  getResPhewas <- function(phewasObj) {
    # note, following do.call(rbind, lapply(...)) constructs generate R NULL when a1 has zero rows
    # okay but need to check below
    
    #extract data from object
    v1 <- phewasObj@v1Data
    loop_clean_db <- phewasObj@loop_clean_db
    connectionType <- phewasObj@conectionType
    all_analyses <- phewasObj@all_analyses
    phewasObj <- phewasObj@analysis
    phewasObj <- phewasObj@analysis_not
    description_contains <- phewasObj@description_contains
    phenotype_contains <- phewasObj@phenotype_contains
    ncase_ge <- phewasObj@ncase_ge
    ncohort_ge <- phewasObj@ncohort_ge
    if ('nearby' %in% names(phewasObj@phewasDataArgs)) {
      nearby <- phewasObj@phewasDataArgs$nearby
    } else {
      nearby <- NULL
    }
    
    
    if (is.null(v1)) {
      # falling through when unique variant not specified
      res <- NULL
    } else if (nrow(a1) == 0L) {
      # no analyses identified using TABLE analyses, don't run real query
      res <- NULL
    } else if (all_analyses) {
      ## Optimize for the case where all analyses are desired, to avoid having a
      ## very long SQL string with thousands of 'OR analysis=' clauses
      
      res <-
        do.call(rbind, lapply(loop_clean_db, function(results_db) {
          outputLoop <- switch(connectionType,
                               'SQL' = getRes(results_db, v1, uniq = FALSE,
                                              zrok = TRUE, connectionType = 'SQL'),
                               stop('Unkown connectionType'))
        }))
      
      
      if (!is.null(nearby)) {
        ## FIXME check nearby is an integer
        res_nearby <-
          do.call(rbind, lapply(loop_clean_db, function(results_db) {
            outputLoop <- switch(connectionType,
                                 'SQL' = getResLoop(results_db, v1, nearby,
                                                    uniq = FALSE, zrok = TRUE,
                                                    connectionType = 'SQL'),
                                 stop('Unkown connectionType'))
            }))
        
        res <- merge(res, res_nearby, all.x = TRUE, all.y = FALSE)
      }
    } else {
      if (nrow(a1) > 10) {
        # FIXME this should be a tuning parameter, unclear what near-optimal value is
        # Optimize for large number of analyses expected
        # evaluate WHERE clause because missing(analysis) etc do not work correctly inside anonymous function body
        w1 <-
          gtxwhat(
            analysis = analysis,
            analysis_not = analysis_not,
            description_contains = description_contains,
            phenotype_contains = phenotype_contains,
            ncase_ge = ncase_ge,
            ncohort_ge = ncohort_ge,
            tablename = 'phewas_results'
          )
        # FIXME, temporary workaround, analyses_tags should (?) be in VIEW definition for phewas_results
        # FIXME current LEFT JOIN analyses_tags seems to result in inefficient query plans
        if (!missing(has_tag)) {
          w1 <-
            paste(w1,
                  'AND',
                  gtxwhat(has_tag = has_tag, tablename = 'analyses_tags'))
        }
      } else {
        # Optimize for small number of analyses expected
        w1 <-
          gtxwhat(analysis = a1$analysis, tablename = 'phewas_results')
      }
      gtx_debug('PheWAS using: {w1}')
      res <-
        do.call(rbind, lapply(loop_clean_db, function(results_db) {
          outputLoop <- switch(connectionType,
                               'SQL' = getResCase2(results_db, v1, w1, 
                                                   uniq = FALSE, zrok = TRUE,
                                                   connectionType),
                               stop('Unknown connectionType'))
        }))
      if (!is.null(nearby)) {
        ## FIXME check nearby is an integer
        res_nearby <-
          do.call(rbind, lapply(loop_clean_db, function(results_db) {
            outputLoop <- switch(connectionType,
                                 'SQL' = getResCase3(results_db, v1, nearby, w1,
                                                     uniq = FALSE, zrok = TRUE,
                                                     connectionType),
                                 stop('Unkown connectionType'))
          }))
        res <- merge(res, res_nearby, all.x = TRUE, all.y = FALSE)
      }
    }
    if (is.null(res)) {
      res <- data.frame(
        analysis = character(0),
        entity = character(0),
        beta = double(0),
        se = double(0),
        pval = double(0),
        rsq = double(0),
        freq = double(0)
      )
    }
    ## Merge results of phewas query with metadata about analyses
    ## Note analyses missing either results or metadata are dropped
    ## [in future we may wish to have a facility to include in forest plots,
    ##  empty lines for phenotypes with missing data]
    res <- within(merge(a1, res, by = 'analysis'), {
      results_db <- NULL
      has_access <- NULL
    })
    
    # add chrom/pos/ref/alt columns to make easier to pass through to other functions, or to help interpretation if saved to a file
    # (note, not added if v1<-NULL above)
    res$chrom <- v1$chrom
    res$pos <- v1$pos
    res$ref <- v1$ref
    res$alt <- v1$alt
    
    # Fix labels for entities, finding labels for the NA ones then pasting on
    tmpl <-
      within(annot(na.omit(res[, c('entity', 'entity_type')])), {
        entity_label <- paste0(label, ' ')
        label <- NULL
      })
    res <- within(merge(res,
                        tmpl[, c('entity', 'entity_type', 'entity_label')],
                        all.x = TRUE, all.y = FALSE),
                  {
                    label <-
                      paste0(ifelse(is.na(entity_label), '', entity_label), label)
                  })
    
    # need to sort AFTER merge with labels
    res <- res[!is.na(res$pval), ]
    res <- res[order(res$pval), ]
    
    # add attribute, used for plot labelling etc., FIXME should we do this when multiple variants matched hence no query was run?
    attr(res, 'variant') <- v1_label
    return(res)
  }
  
  
  return(phewasObj)
  
  
}

#' @describe phewas Draw simple QQ plot of PheWAS results
#'
#' @inheritParams phewas
#'
#' @usage
#' phewas.qq(chrom, pos, ref, alt, rs, analysis, analysis_not,
#' description_contains, phenotype_contains, ncase_ge, ncohort_ge,
#' dbc = getOption("gtx.dbConnection", NULL))
#'
#' @export
phewas.qq <- function(chrom,
                      pos,
                      ref,
                      alt,
                      rs,
                      analysis,
                      analysis_not,
                      description_contains,
                      phenotype_contains,
                      ncase_ge,
                      ncohort_ge,
                      dbc = getOption("gtx.dbConnection", NULL)) {
  gtxdbcheck(dbc)
  
  res <-
    phewas.data(
      chrom = chrom,
      pos = pos,
      ref = ref,
      alt = alt,
      rs = rs,
      analysis = analysis,
      analysis_not = analysis_not,
      description_contains = description_contains,
      phenotype_contains = phenotype_contains,
      ncase_ge = ncase_ge,
      ncohort_ge = ncohort_ge,
      dbc = dbc
    )
  
  with(res, qq10(res$pval))
  
  return(invisible(res))
}
