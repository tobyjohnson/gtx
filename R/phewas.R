#' Phenome Wide Assocation Study
#' 
#' Look up and plot results for a Phenome Wide Association Study.
#' 
#' Look up association summary statistics for a single variant of interest, specified 
#' either by chromosome and position, or by dbSNP rs identifier.  (FIXME need to add optional 
#' ref/alt arguments to disambiguate when there are multiple variants at the same position,
#' or multiple variants with the same rs.)
#' 
#' By default, summary statistics are returned for all analyses where the variant of interest 
#' was tested.  If results are required for only a specific set of analyses, these can be 
#' specified using a combination of arguments that are passed through to \code{gtxanalyses()}.
#' 
#' @param plot_ymax Y-axis value maximum 
#' 
#' @author Toby Johnson \email{Toby.x.Johnson@gsk.com}
#' @export
# Could inherit param documentation from gtxanalyses
# Should merge phewas.qq into main phewas as an alternative plot style
phewas <- function(chrom, pos, ref, alt, rs,
                   analysis, analysis_not, 
                   description_contains,
                   phenotype_contains,
                   has_tag, 
                   ncase_ge,
                   ncohort_ge,
                   ## if extra filters are added, be sure to update definition of all_analyses below                   
                   analysis_fields = c('description', 'label', 'unit',
                                       'ncase', 'ncontrol', 'ncohort'),
                   tag_is = 'phewas_group', 
                   plot_ymax = 30,
                   dbc = getOption("gtx.dbConnection", NULL)) {

    gtxdbcheck(dbc)
    p1 <- phewas.data(chrom = chrom, pos = pos, ref = ref, alt = alt, rs = rs,
                      analysis = analysis, analysis_not = analysis_not,
                      description_contains = description_contains,
                      phenotype_contains = phenotype_contains,
                      has_tag = has_tag, 
                      ncase_ge = ncase_ge, ncohort_ge = ncohort_ge,
                      analysis_fields = analysis_fields,
                      tag_is = tag_is, with_tags = TRUE, # force with_tags = TRUE
                      dbc = dbc)
    p1orig <- p1

    ymax <- max(10, ceiling(-log10(min(p1$pval))))
    if (ymax > plot_ymax) { # Hard coded threshold makes sense for control of visual display
        ymax <- plot_ymax
        warning('Truncating P-values at 1e-', ymax)
        truncp <- TRUE
    } else {
        truncp <- FALSE
    }

    p1 <- p1[order(p1$tag, p1$pval), ] # with this sort, tag==NA will be last
    p1$tag[is.na(p1$tag)] <- 'No tag'
    p1$x1 <- 1:nrow(p1) # initial attempt at x axis position, to be refined below
    x1 <- with(aggregate(p1$x1, 
                         by = list(tag = p1$tag), 
                         FUN = function(x) return(c(min = min(x), max = max(x)))),
               cbind(data.frame(tag), as.data.frame(x)))
    group_interspace <- 50 # this should be a user configurable parameter
    x1$offset <- c(0, cumsum(x1$max - x1$min + group_interspace))[1:nrow(x1)] - x1$min + group_interspace
    x1$midpt <- 0.5*(x1$max + x1$min) + x1$offset
    if (x1$tag[nrow(x1)] == 'No tag') {
        x1$col <- c(rainbow(nrow(x1) - 1), 'grey') # rainbow for all tags, grey for last
    } else {
        x1$col <- rainbow(nrow(x1)) # rainbow for all tags
    }
    p1m <- match(p1$tag, x1$tag)
    p1$x <- p1$x1 + x1$offset[p1m]
    p1$col <- x1$col[p1m]
    # note still sorted by tag then pval
    p1$y <- pmin(-log10(p1$pval), ymax)
    p1$yo <- c(Inf, ifelse(p1$tag[-1] != p1$tag[-nrow(p1)], Inf, p1$y[-nrow(p1)] - p1$y[-1])) # yaxis space within group relative to previous point

    p1 <- p1[order(p1$pval, decreasing = TRUE), ] # reorder to make sure most significant points plotted over less significant ones
    plot.new()
    plot.window(range(p1$x), c(0, max(p1$y)))

    max_y = max(p1$y)
    ## plot text then overplot with points
    with(subset(p1, y > min(6, max_y*.96) & yo > strheight('M', cex = 0.5)), # should be user controllable
         text(x, y, as.character(label), pos = 4, cex = 0.5))
    points(p1$x, p1$y, pch = 19, col = p1$col, cex = 0.5)
    if (truncp) {
        ## would be nice to more cleanly overwrite y axis label
        axis(2, at = ymax, labels = substitute({}>=ymax, list(ymax = ymax)), las = 1)
        # could e.g. do setdiff(pretty(...), ymax)
    }
    axis(2, las = 1)
    mtext(x1$tag, side = 1, line = 0.5, las = 2, at = x1$midpt, col = x1$col, las = 2, cex = 0.5)
    box()
    
    return(invisible(p1orig))
}

#' @describeIn phewas Look up results without plotting
#' 
#' @param chrom Chromosome of variant of interest
#' @param pos Position of variant of interest
#' @param ref Reference allele of variant of interest
#' @param alt Alternate allele of variant of interest
#' @param rs dbSNP rs identifier of variant of interest
#' @param analysis Identifiers of analyses to include
#' @param analysis_not Identifiers of analyses to exclude
#' @param description_contains Search term for analyses to include
#' @param has_tag Tag value for analyses to include
#' @param ncase_ge Threshold case sample size greater-or-equal for inclusion 
#' @param ncohort_ge Threshold cohort sample size greater-or-equal for inclusion
#' @param analysis_fields Fields from analysis metadata to return
#' @param tag_is Tag flag to use for plot grouping
#' @param nearby Additionally return minimum pvalue at position+-nearby
#' @param dbc Database connection
#'
#' @return Data frame of summary statistics
#' @export
phewas.data <- function(chrom, pos, ref, alt, rs,
                        analysis, analysis_not, 
                        description_contains,
                        phenotype_contains,
                        has_tag, 
                        ncase_ge,
                        ncohort_ge,
                        ## if extra filters are added, be sure to update definition of all_analyses below
                        analysis_fields = c('description', 'label', 'unit',
                                    'ncase', 'ncontrol', 'ncohort'),
                        tag_is, with_tags = FALSE, # not clear why missing with_tags cannot be passed through gtxanalyses()
                        nearby, 
                        dbc = getOption("gtx.dbConnection", NULL)) {
    gtxdbcheck(dbc)

    ## Look up variant
    v1 <- sqlWrapper(dbc,
                     sprintf('SELECT chrom, pos, ref, alt FROM sites WHERE %s;',
                             gtxwhere(chrom = chrom, pos = pos, ref = ref, alt = alt, rs = rs))) # default uniq = TRUE
    flog.debug(paste0('Resolved PheWAS varint to chr', v1$chrom, ':', v1$pos, ':', v1$ref, '>', v1$alt))

    ## Look up analysis metadata
    a1 <- gtxanalyses(analysis = analysis, analysis_not = analysis_not, 
                      description_contains = description_contains,
                      phenotype_contains = phenotype_contains,
                      has_tag = has_tag, 
                      ncase_ge = ncase_ge, ncohort_ge = ncohort_ge,
                      analysis_fields = analysis_fields, 
                      tag_is = tag_is, with_tags = with_tags, 
                      has_access_only = TRUE, 
                      dbc = dbc) # will work fine if all filtering arguments are missing, as it internally sets all_analyses<-TRUE
    flog.debug(paste0('Queried metadata for ', nrow(a1), ' analyses with access to'))

    ## Handle when results_db is NULL in database (returned as NA to R), since not using gtxanalysisdb()
    ## Add period if results_db is a database name, otherwise empty string (use pattern %sgwas_results in sprintf's below)
    ## (FIXME it would be better if this did go through gtxanalysisdb() for future maintenance)
    loop_clean_db <- sapply(unique(a1$results_db), function(x) {
      return(if (is.na(x) || x == '') '' else sanitize1(paste0(x, ''), type = 'alphanum.'))
    })

    all_analyses <- (missing(analysis) && missing(analysis_not) && missing(phenotype_contains) &&
                     missing(description_contains) && missing(has_tag) && 
                     missing(ncase_ge) && missing(ncohort_ge))

    # note, following do.call(rbind, lapply(...)) constructs generate R NULL when a1 has zero rows
    # okay but need to check below

    if (all_analyses) {
        ## Optimize for the case where all analyses are desired, to avoid having a
        ## very long SQL string with thousands of 'OR analysis=' clauses
        res <- do.call(rbind, lapply(loop_clean_db, function(results_db) {
            flog.debug(paste0('PheWAS all_analyses query in db ', results_db, ' ...'))
            sqlWrapper(getOption('gtx.dbConnection'),
                       sprintf('SELECT analysis, entity, beta, se, pval, rsq, freq FROM %sgwas_results WHERE chrom=\'%s\' AND pos=%s AND ref=\'%s\' AND alt=\'%s\' AND pval IS NOT NULL;',
                               results_db,
                               sanitize1(v1$chrom, values = c(as.character(1:22), "X", "Y")),
                               sanitize1(v1$pos, type = "int"),
                               sanitize1(v1$ref, type = "ACGT+"),
                               sanitize1(v1$alt, type = "ACGT+")),
                       uniq = FALSE, zrok = TRUE)
        }))
        if (!missing(nearby)) {
            ## FIXME check nearby is an integer
            res_nearby <- do.call(rbind, lapply(loop_clean_db, function(results_db) {
                sqlWrapper(getOption('gtx.dbConnection'),
                           sprintf('SELECT analysis, entity, min(pval) AS pval_nearby FROM %sgwas_results WHERE %s AND pval IS NOT NULL GROUP BY analysis, entity;',
                                   results_db,
                                   gtxwhere(chrom = v1$chrom, pos_ge = v1$pos - nearby, pos_le = v1$pos + nearby)),
                           uniq = FALSE, zrok = TRUE)
            }))
            res <- merge(res, res_nearby, all.x = TRUE, all.y = FALSE)
        }
    } else {
        res <- do.call(rbind, lapply(loop_clean_db, function(results_db) {
            sqlWrapper(getOption('gtx.dbConnection'),
                       sprintf('SELECT analysis, entity, beta, se, pval, rsq, freq FROM %sgwas_results WHERE %s AND chrom=\'%s\' AND pos=%s AND ref=\'%s\' AND alt=\'%s\' AND pval IS NOT NULL;',
                               results_db,
                               gtxwhat(analysis = a1$analysis),
                               sanitize1(v1$chrom, values = c(as.character(1:22), "X", "Y")),
                               sanitize1(v1$pos, type = "int"),
                               sanitize1(v1$ref, type = "ACGT+"),
                               sanitize1(v1$alt, type = "ACGT+")),
                       uniq = FALSE, zrok = TRUE)
        }))
        if (!missing(nearby)) {
            ## FIXME check nearby is an integer
            res_nearby <- do.call(rbind, lapply(loop_clean_db, function(results_db) {
                sqlWrapper(getOption('gtx.dbConnection'),
                           sprintf('SELECT analysis, entity, min(pval) AS pval_nearby FROM %sgwas_results WHERE %s AND %s AND pval IS NOT NULL GROUP BY analysis, entity;',
                                   results_db,
                                   gtxwhat(analysis = a1$analysis),
                                   gtxwhere(chrom = v1$chrom, pos_ge = v1$pos - nearby, pos_le = v1$pos + nearby)),
                           uniq = FALSE, zrok = TRUE)
            }))
            res <- merge(res, res_nearby, all.x = TRUE, all.y = FALSE)
        }
    }
    if (is.null(res)) {
        res <- data.frame(analysis = character(0), entity = character(0), 
                          beta = double(0), se = double(0), pval = double(0), rsq = double(0), freq = double(0))
    }
    ## Merge results of phewas query with metadata about analyses
    ## Note analyses missing either results or metadata are dropped
    ## [in future we may wish to have a facility to include in forest plots,
    ##  empty lines for phenotypes with missing data]
    res <- within(merge(a1, res, by = 'analysis'), {
        results_db <- NULL
        has_access <- NULL
    })
    
    res <- res[!is.na(res$pval),]
    res <- res[order(res$pval), ]
    return(res)
}

#' @describeIn phewas Draw simple QQ plot of PheWAS results

#' @export 
phewas.qq <- function(chrom, pos, ref, alt, rs,
                      analysis, analysis_not, 
                      description_contains,
                      phenotype_contains,
                      ncase_ge,
                      ncohort_ge,
                      dbc = getOption("gtx.dbConnection", NULL)) {
    gtxdbcheck(dbc)

    res <- phewas.data(chrom = chrom, pos = pos, ref = ref, alt = alt, rs = rs,
                       analysis = analysis, analysis_not = analysis_not, 
                       description_contains = description_contains,
                       phenotype_contains = phenotype_contains,
                       ncase_ge = ncase_ge, ncohort_ge = ncohort_ge,
                       dbc = dbc)
    
    with(res, qq10(res$pval))

    return(invisible(res))
}
