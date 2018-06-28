## Internal functions for finemapping (fm) calculations

## The user interface to finemapping functionality is regionplot() and regionplot.data()

## Detailed multisignal results are exposed via cleo()

## Trying to use the following notation
## ppc = posterior probability of causality
## ccs = credibly causal set
## _signal = computed uasing single signal assumption
## _cleo = computed using CLEO analyses

## Note, rough exploration with AD meta-analysis of IGAP and UKB results
## suggests that, for lnodds units, and with priorsd=1, priorc=1e-05 gives
## sensible results for nullppc, as a function of lead variant P-value,
## across a range of highly-borderline-non-significant loci
## FIXME make a small notebook documenting this
## Broader discussion with StatGen, change the default priors for coloc?

fm_signal <- function(data,
                      priorsd = 1, priorc = 1e-5, ccs_size = 0.95, ccs_only = FALSE) {
    stopifnot(all(c('beta', 'se') %in% names(data)))
    data$lnabf <- with(data, abf.Wakefield(beta, se, priorsd, log = TRUE))
    ##abf <- norm1(c(0, data$lnabf), log = TRUE) # general decision not to return this intermediate quantity
    ppc <- norm1(c(0, log(priorc) + data$lnabf), log = TRUE)
    ##abf[is.na(abf)] <- 0.
    ppc[is.na(ppc)] <- 0. # should be option to leave as NA or set to zero
    ##data$abf_signal <- abf[-1]
    data$ppc_signal <- ppc[-1]
    data$ccs_signal <- credset(ppc, cred = ccs_size)[-1]
    if (ccs_only) data <- data[which(data$ccs_signal), ]
    attr(data, 'params_signal') <- list(priorsd = priorsd, priorc = priorc, ccs_size = ccs_size)
    ##attr(data, 'nullabf_signal') <- abf[1]
    attr(data, 'nullppc_signal') <- ppc[1]
    return(data)
}

## currently a single function to query all conditional signals that
## (partly or fully) overlap a query region, and compute causal probs
## currently, selection by entity is not supported
fm_cleo.data <- function(analysis,
                         chrom, pos_start, pos_end, pos, 
                         hgncid, ensemblid, rs, surround = 500000,
                         dbc = getOption("gtx.dbConnection", NULL)) {
    ## check database connection
    gtxdbcheck(dbc)

    xregion <- gtxregion(chrom = chrom, pos_start = pos_start, pos_end = pos_end, pos = pos, 
                         hgncid = hgncid, ensemblid = ensemblid, rs = rs, surround = surround,
                         dbc = dbc)
    
    gtxlog('Querying distinct signals from CLEO analyses')
    t0 <- as.double(Sys.time())
    ds <- sqlWrapper(dbc,
                     sprintf('SELECT DISTINCT signal FROM %s.gwas_results_cond WHERE %s AND %s;', 
                             gtxanalysisdb(analysis), 
                             gtxwhat(analysis1 = analysis),
                             gtxwhere(chrom = xregion$chrom, pos_ge = xregion$pos_start, pos_le = xregion$pos_end)),
                     uniq = FALSE, zrok = TRUE)$signal
    t1 <- as.double(Sys.time())
    gtxlog('Query returned ', length(ds), ' signal(s) overlapping query region ', xregion$label, ' in ', round(t1 - t0, 3), 's.')
    if (length(ds) > 0) {
        gtxlog('Querying summary statistics from CLEO analyses')
        t0 <- as.double(Sys.time())
        ss <- sqlWrapper(dbc,
                         sprintf('SELECT chrom,pos,ref,alt,signal,beta_cond,se_cond FROM %s.gwas_results_cond WHERE %s AND %s AND (%s);', 
                                 gtxanalysisdb(analysis), 
                                 gtxwhat(analysis1 = analysis),
                                 gtxwhere(chrom = xregion$chrom),
                                 paste0('signal=', sanitize(ds, type = 'int'), collapse = ' OR ')),
                           uniq = FALSE, zrok = TRUE)
        t1 <- as.double(Sys.time())
        gtxlog('Query returned ', nrow(ss), ' summary statistics in ', round(t1 - t0, 3), 's.')
        ## modify xregion, extend to include the range of positions actually spanned by these signals
        xregion$pos_start <- min(ss$pos)
        xregion$pos_end <- max(ss$pos)
    } else {
        ss <- sqlWrapper(dbc,
                         sprintf('SELECT chrom,pos,ref,alt,signal,beta_cond,se_cond FROM %s.gwas_results_cond LIMIT 0;', 
                                 gtxanalysisdb(analysis)), 
                         uniq = FALSE, zrok = TRUE)
    }
    attr(ss, 'analysis') <- analysis
    attr(ss, 'region') <- xregion
    return(ss)
}

## currently, selection by entity is not supported
fm_cleo <- function(analysis,
                    chrom, pos_start, pos_end, pos, 
                    hgncid, ensemblid, rs, surround = 500000,
                    priorsd = 1., priorc = 1e-5, ccs_size = 0.95, ccs_only = FALSE,
                    dbc = getOption("gtx.dbConnection", NULL)) {
    ## database connection is checked within fm_cleo.data() call

    data <- fm_cleo.data(analysis = analysis,
                         chrom = chrom, pos_start = pos_start, pos_end = pos_end, pos = pos,
                         hgncid = hgncid, ensemblid = ensemblid, rs = rs, surround = surround, 
                         dbc = dbc)

    if (nrow(data) > 0) {
        ## compute lnABF+ln(priorc) for each variant for each signal,
        ## relative to lnabf==lnabfp==0. for each signal null
        data$lnabf <- with(data, abf.Wakefield(beta_cond, se_cond, priorsd = priorsd, log = TRUE))
        data$lnabfp <- data$lnabf + log(priorc)
        ## rescale so max is 0, over variants *and the null*
        ## including null should not be needed *here*, if CLEO has selected significant signals, unless priorc is *very* small
        ## FIXME this is a more stable scaling to use in other fm calculations
        maxlnabfp <- aggregate(data[ , 'lnabfp', drop = FALSE], by = list(signal = data$signal), FUN = max, na.rm = TRUE)
        maxlnabfp$lnabfp <- pmax(maxlnabfp$lnabfp, 0.)
        ## compute unnormalized ppc (ppcu)
        data$ppcu <- exp(data$lnabfp - maxlnabfp$lnabfp[match(data$signal, maxlnabfp$signal)])
        ## compute sum over variants and nulls within signals, normalize
        sumppcu <- aggregate(data[ , 'ppcu', drop = FALSE], by = list(signal = data$signal), FUN = sum, na.rm = TRUE)
        sumppcu$nullppcu <- exp(0. - maxlnabfp$lnabfp)[match(sumppcu$signal, maxlnabfp$signal)]
        sumppcu$sumppcu <- sumppcu$ppcu + sumppcu$nullppcu
        sumppcu$nullppc <- sumppcu$nullppcu/sumppcu$sumppcu
        data$ppc_cleo <- data$ppcu/sumppcu$sumppcu[match(data$signal, sumppcu$signal)]
        ## clean up data
        data$lnabfp <- NULL
        data$ppcu <- NULL
        ## extract null ppcs
        nullppc <- sumppcu$nullppc 
        names(nullppc) <- sumppcu$signal
        ## compute credibly causal sets
        data$ccs_cleo <- NA
        for (sig in sumppcu$signal) {
            ww <- which(data$signal == sig)
            data$ccs_cleo[ww] <- gtx:::credset(c(nullppc[match(sig, names(nullppc))], data$ppc[ww]), cred = ccs_size)[-1]
        }
        ## add as attr()ibutes
        attr(data, 'params_cleo') <- list(priorsd = priorsd, priorc = priorc, ccs_size = ccs_size)
        attr(data, 'nullppc_cleo') <- nullppc
        if (ccs_only) data <- subset(data, ccs_cleo)
        return(data)
    } else {
        ## data have no rows, FIXME add expected columns and attr()ibutes
        return(data)
    }
}

cleo <- function(analysis,
                 chrom, pos_start, pos_end, pos, 
                 hgncid, ensemblid, rs, surround = 500000,
                 priorsd = 1., priorc = 1e-5, ccs_size = 0.95, 
                 dbc = getOption("gtx.dbConnection", NULL)) {
    ## database connection is checked within fm_cleo.data() call

    ## If in future, priorsd and priorc are columns in analyses
    ## *here* is the place we would query and use those values

    ## call fm_cleo with ccs_only = FALSE
    data <- fm_cleo(analysis = analysis,
                    chrom = chrom, pos_start = pos_start, pos_end = pos_end, pos = pos,
                    hgncid = hgncid, ensemblid = ensemblid, rs = rs, surround = surround,
                    priorsd = priorsd, priorc = priorc, ccs_size = ccs_size, ccs_only = FALSE,
                    dbc = dbc)
    xregion <- attr(data, 'region')

    ## annotate data by signals
    gtxlog('Merging with protein coding annotation')
    t0 <- as.double(Sys.time())
    vep <- sqlWrapper(dbc,
                      sprintf('SELECT chrom, pos, ref, alt, impact FROM vep WHERE %s;',
                              gtxwhere(chrom = xregion$chrom, pos_ge = xregion$pos_start, pos_le = xregion$pos_end)),
                      uniq = FALSE, zrok = TRUE)    
    t1 <- as.double(Sys.time())
    gtxlog('Query returned ', nrow(vep), ' annotations within query region ', xregion$label, ' in ', round(t1 - t0, 3), 's.')
    t0 <- as.double(Sys.time())
    ## note, merge does not preserve attr()ibutes
    datas <- merge(data, vep, by = c('chrom', 'pos', 'ref', 'alt'), all.x = TRUE, all.y = FALSE) 
    ## FIXME there's an interesting possibility to use all.y=TRUE to capture impactful variants that are not in the FM analysis
    ## and potentially flag this to the user
    t1 <- as.double(Sys.time())
    gtxlog('Merged with CLEO results in ', round(t1 - t0, 3), 's.')

    ## PCI is probability of causal impact
    datai <- subset(datas, !is.na(impact), drop = FALSE) # FIXME when content of impact changes
    if (nrow(datai) > 0) {
        pci1 <- aggregate(datai[, 'ppc_cleo', drop = FALSE], 
                         by = list(signal = datai$signal),
                         FUN = sum, na.rm = TRUE)
        pci1 <- pci1[order(pci1$ppc_cleo, decreasing = TRUE), ]
        pci <- pci1$ppc_cleo
        names(pci) <- pci1$signal
        ## FIXME do setdiff with names(attr(data, 'nullppc_cleo')) to add true zeros
    } else {
        pci <- NULL # FIXME should be all zero with signals from names(attr(data, 'nullppc_cleo'))
    }

    ## aggregate over signals
    gtxlog('Aggregating over signals')
    t0 <- as.double(Sys.time())
    dsumppc <- aggregate(data[ , 'ppc_cleo', drop = FALSE], 
                         by = with(data, list(chrom = chrom, pos = pos, ref = ref, alt = alt)),
                         FUN = sum, na.rm = TRUE)
    dorccs <- aggregate(data[ , 'ccs_cleo', drop = FALSE], 
                        by = with(data, list(chrom = chrom, pos = pos, ref = ref, alt = alt)),
                        FUN = any)
    dataa <- merge(merge(dsumppc, dorccs, by = c('chrom', 'pos', 'ref', 'alt'), all = TRUE),
                   vep, by = c('chrom', 'pos', 'ref', 'alt'), all.x = TRUE, all.y = FALSE)
    dataa <- dataa[with(dataa, order(chrom, pos, ref, alt)), ]
    t1 <- as.double(Sys.time())
    gtxlog('Aggregated in ', round(t1 - t0, 3), 's.')

    datas = subset(datas, ccs_cleo) # drop to signalwise CCS
    
    ##Preserve original attributes, note this list should be introspected not hard coded (FIXME)
    for (a in c('analysis', 'region', 'params_cleo', 'nullppc_cleo')) {
        attr(datai, a) <- attr(data, a)
        attr(dataa, a) <- attr(data, a)
    }
    
    return(list(nsignals = length(attr(data, 'nullppc_cleo')),
                pci = pci, 
                signals = datas,
                aggregate = dataa))
}
