## Internal functions for finemapping (fm) calculations

## The user interface to finemapping functionality is regionplot() and regionplot.data()

## Detailed multisignal results are exposed via cleo()

## Trying to use the following notation in user-facing objects:
## pp = posterior probability (of causality, for one variant)
## nullpp = null/emptry model posterior probability
## cs = credible set (of models, both with variants and empty/null)
## _signal = computed using single signal assumption
## _cleo = computed using CLEO analyses
##
## lnabf = ln(ABF), for model with one variant, versus null/empty model
## Note that lnabf>0 indicates support for model with one variant, unlike original Wakefield definition

## Consistent notation for internal calculations:
## lnabf = ln(ABF) as above
## lnabfp = ln(ABF) + ln(prior), with relative prior for variant added, relative to prior=1 for null/empty model
##    hence lnabfp==0 for the null/empty model
## ppu = exp(lnabfp - K), posterior probability (unnormalized), where K is a constant within signals
## maxlnabfp = usual method for choosing K, max over variants *and null*
## sumppu = denominator for calculating pp, sum over variants *and null*

## Note, rough exploration with AD meta-analysis of IGAP and UKB results
## suggests that, for lnodds units, and with priorsd=1, priorc=1e-05 gives
## sensible results for nullpp, as a function of lead variant P-value,
## across a range of highly-borderline-non-significant loci
## FIXME make a small notebook documenting this
## Broader discussion with StatGen, change the default priors for coloc?

##' @export
fm_signal <- function(data,
                      priorsd = 1, priorc = 1e-5, cs_size = 0.95, cs_only = FALSE) {
  stopifnot(all(c('beta', 'se') %in% names(data)))
    
  ## Compute lnABF for each variant versus empty model
  data$lnabf <- with(data, abf.Wakefield(beta, se, priorsd, log = TRUE))
  ## Make a vector of (ABFs and) posteriors *with empty model as first element*
  ## (*) Thus the i-th variant has position (i+1) in this vector
  ## Normalization sum-to-one is thus over the empty model and all one-causal-variant models
  ##abf <- norm1(c(0, data$lnabf), log = TRUE) # general decision not to return this intermediate quantity
  pp <- norm1(c(0, data$lnabf + log(priorc)), log = TRUE)
  ##abf[is.na(abf)] <- 0.
  pp[is.na(pp)] <- 0. # should be option to leave as NA or set to zero
  ## Add normalised probabilities for variants back to dataframe;
  ## drop first element see note (*) above
  ##data$abf_signal <- abf[-1]
  data$pp_signal <- pp[-1]
  ## Compute credible set on whole vector including empty model, then
  ## add TRUE/FALSE indicators for variants back to the dataframe
  ## drop first element see note (*) above
  data$cs_signal <- credset(pp, cred = cs_size)[-1]
  if (cs_only) data <- data[which(data$cs_signal), ]
  attr(data, 'params_signal') <- list(priorsd = priorsd, priorc = priorc, cs_size = cs_size)
  ##attr(data, 'nullabf_signal') <- abf[1]
  attr(data, 'nullpp_signal') <- pp[1]
  return(data)
}

## internal function to query all conditional signals that (partly or
## fully) overlap a query region
## FIXME: currently, selection by entity is not supported
fm_cleo.data <- function(analysis,
                         chrom, pos_start, pos_end, pos, 
                         hgncid, ensemblid, rs, surround = 500000,
                         dbc = getOption("gtx.dbConnection", NULL)) {
    ## check database connection
    gtxdbcheck(dbc)

    xregion <- gtxregion(chrom = chrom, pos_start = pos_start, pos_end = pos_end, pos = pos, 
                         hgncid = hgncid, ensemblid = ensemblid, rs = rs, surround = surround,
                         dbc = dbc)
    
    futile.logger::flog.info('Querying distinct signals from CLEO analyses')
    t0 <- as.double(Sys.time())
    ds <- sqlWrapper(dbc,
                     sprintf('SELECT DISTINCT signal FROM %sgwas_results_cond WHERE %s AND %s;', 
                             gtxanalysisdb(analysis), 
                             gtxwhat(analysis1 = analysis),
                             gtxwhere(chrom = xregion$chrom, pos_ge = xregion$pos_start, pos_le = xregion$pos_end)),
                     uniq = FALSE, zrok = TRUE)$signal
    t1 <- as.double(Sys.time())
    futile.logger::flog.info(paste0('Query returned ', length(ds), ' signal(s) overlapping query region ', xregion$label, ' in ', round(t1 - t0, 3), 's.'))
    if (length(ds) > 0) {
        futile.logger::flog.info('Querying summary statistics from CLEO analyses')
        # First get summary statistics for index variants (ssi)
        t0 <- as.double(Sys.time())
        ssi <- sqlWrapper(dbc,
                       sprintf('SELECT t1.chrom, t1.pos, t1.ref, t1.alt, signal, beta_joint, se_joint, beta, se, pval FROM %sgwas_results AS t1 JOIN %sgwas_results_joint AS t2 USING (chrom, pos, ref, alt) WHERE %s AND %s;', 
                               gtxanalysisdb(analysis),
                               gtxanalysisdb(analysis), 
                               where_from(analysisu = analysis, chrom = xregion$chrom, tablename = 't1'),
                               where_from(analysisu = analysis, chrom = xregion$chrom, signal = ds, tablename = 't2')),
                       uniq = FALSE, zrok = TRUE)
        t1 <- as.double(Sys.time())
        futile.logger::flog.info(paste0('Index variant query returned ', nrow(ssi), ' summary statistics in ', round(t1 - t0, 3), 's.'))
        t0 <- as.double(Sys.time())
        ss <- sqlWrapper(dbc,
                         sprintf('SELECT chrom,pos,ref,alt,signal,beta_cond,se_cond FROM %sgwas_results_cond WHERE %s;', 
                                 gtxanalysisdb(analysis), 
                                 where_from(analysisu = analysis, chrom = xregion$chrom, signal = ds)),
                           uniq = FALSE, zrok = TRUE)
        t1 <- as.double(Sys.time())
        futile.logger::flog.info(paste0('CLEO query returned ', nrow(ss), ' summary statistics in ', round(t1 - t0, 3), 's.'))
        ## modify xregion, extend to include the range of positions actually spanned by these signals
        xregion$pos_start <- min(ss$pos)
        xregion$pos_end <- max(ss$pos)
    } else {
        ## FIXME db query with LIMIT 0 is an expensive way to generate a df with 0 rows
        ss <- sqlWrapper(dbc,
                         sprintf('SELECT chrom,pos,ref,alt,signal,beta_cond,se_cond FROM %sgwas_results_cond LIMIT 0;', 
                                 gtxanalysisdb(analysis)), 
                         uniq = FALSE, zrok = TRUE)
    }
    attr(ss, 'analysis') <- analysis
    attr(ss, 'region') <- xregion
    attr(ss, 'index_cleo') <- ssi
    return(ss)
}

## internal function to calculate posterior probabilities and credible
## sets from data returned by fm_cleo.data()
## FIXME: currently, selection by entity is not supported
fm_cleo <- function(analysis,
                    chrom, pos_start, pos_end, pos, 
                    hgncid, ensemblid, rs, surround = 500000,
                    priorsd = 1., priorc = 1e-5, cs_size = 0.95, cs_only = FALSE,
                    dbc = getOption("gtx.dbConnection", NULL)) {
    
    ## no database connection check needed, because it happens within fm_cleo.data()

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
        ## FIXME: this is a more stable scaling, propagate its use in other fm calculations
        maxlnabfp <- aggregate(data[ , 'lnabfp', drop = FALSE], by = list(signal = data$signal), FUN = max, na.rm = TRUE)
        maxlnabfp$lnabfp <- pmax(maxlnabfp$lnabfp, 0.)
        ## compute unnormalized pp (ppu)
        data$ppu <- exp(data$lnabfp - maxlnabfp$lnabfp[match(data$signal, maxlnabfp$signal)])
        ## compute sum over variants and nulls within signals, normalize
        sumppu <- aggregate(data[ , 'ppu', drop = FALSE], by = list(signal = data$signal), FUN = sum, na.rm = TRUE)
        sumppu$nullppu <- exp(0. - maxlnabfp$lnabfp)[match(sumppu$signal, maxlnabfp$signal)]
        sumppu$sumppu <- sumppu$ppu + sumppu$nullppu
        sumppu$nullpp <- sumppu$nullppu/sumppu$sumppu
        data$pp_cleo <- data$ppu/sumppu$sumppu[match(data$signal, sumppu$signal)]
        ## remove unwanted columns from data
        data$lnabfp <- NULL
        data$ppu <- NULL
        ## extract nullpp for each signal
        nullpp <- sumppu$nullpp 
        names(nullpp) <- sumppu$signal
        ## compute credible sets for each signal
        ## FIXME: Optimize to remove loop over repeated which (i.e. groupby).  Slightly complicated by
        ## dependence on nullpp
        ## Consider adding dummy rows with chrom=pos=ref=alt=beta=se=NA to represent the empty model for each signal
        ## This would make the normalizing easier
        data$cs_cleo <- NA
        for (sig in sumppu$signal) {
            ww <- which(data$signal == sig)
            data$cs_cleo[ww] <- gtx:::credset(c(nullpp[match(sig, names(nullpp))], data$pp[ww]), cred = cs_size)[-1]
        }
        ## add attr()ibutes
        attr(data, 'params_cleo') <- list(priorsd = priorsd, priorc = priorc, cs_size = cs_size)
        attr(data, 'nullpp_cleo') <- nullpp

        ## return either the cs only, or the complete data
        if (cs_only) return(subset(data, cs_cleo))
        return(data)
        
    } else {
        ## data have no rows, FIXME add expected columns and attr()ibutes
        return(data)
    }
}

##' @export
cleo <- function(analysis,
                 chrom, pos_start, pos_end, pos, 
                 hgncid, ensemblid, rs, surround = 500000,
                 priorsd = 1., priorc = 1e-5, cs_size = 0.95, 
                 dbc = getOption("gtx.dbConnection", NULL)) {
    ## database connection is checked within fm_cleo.data() call

    ## If in future, priorsd and priorc are columns in analyses
    ## *here* is the place we would query and use those values

    ## call fm_cleo with cs_only = FALSE
    data <- fm_cleo(analysis = analysis,
                    chrom = chrom, pos_start = pos_start, pos_end = pos_end, pos = pos,
                    hgncid = hgncid, ensemblid = ensemblid, rs = rs, surround = surround,
                    priorsd = priorsd, priorc = priorc, cs_size = cs_size, cs_only = FALSE,
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
        pci1 <- aggregate(datai[, 'pp_cleo', drop = FALSE], 
                         by = list(signal = datai$signal),
                         FUN = sum, na.rm = TRUE)
        pci1 <- pci1[order(pci1$pp_cleo, decreasing = TRUE), ]
        pci <- pci1$pp_cleo
        names(pci) <- pci1$signal
        ## FIXME do setdiff with names(attr(data, 'nullpp_cleo')) to add true zeros
    } else {
        pci <- NULL # FIXME should be all zero with signals from names(attr(data, 'nullpp_cleo'))
    }

    ## aggregate over signals
    ## FIXME, this throws error if data has no rows
    gtxlog('Aggregating over signals')
    t0 <- as.double(Sys.time())
    dsumpp <- aggregate(data[ , 'pp_cleo', drop = FALSE], 
                         by = with(data, list(chrom = chrom, pos = pos, ref = ref, alt = alt)),
                         FUN = sum, na.rm = TRUE)
    dorcs <- aggregate(data[ , 'cs_cleo', drop = FALSE], 
                        by = with(data, list(chrom = chrom, pos = pos, ref = ref, alt = alt)),
                        FUN = any)
    dataa <- merge(merge(dsumpp, dorcs, by = c('chrom', 'pos', 'ref', 'alt'), all = TRUE),
                   vep, by = c('chrom', 'pos', 'ref', 'alt'), all.x = TRUE, all.y = FALSE)
    dataa <- dataa[with(dataa, order(chrom, pos, ref, alt)), ]
    t1 <- as.double(Sys.time())
    gtxlog('Aggregated in ', round(t1 - t0, 3), 's.')

    datas = subset(datas, cs_cleo) # drop to signalwise CS
    
    ##Preserve original attributes, note this list should be introspected not hard coded (FIXME)
    for (a in c('analysis', 'region', 'params_cleo', 'nullpp_cleo')) {
        attr(datai, a) <- attr(data, a)
        attr(dataa, a) <- attr(data, a)
    }

    nullpp <- attr(data, 'nullpp_cleo')
    nullpp <- nullpp[order(nullpp)]
    
    return(list(nsignals = length(attr(data, 'nullpp_cleo')),
                nullpp = nullpp,
                altpp = 1. - nullpp, # need to think of a better name for this 'alt model ppc'
                pci = pci, 
                signals = datas,
                aggregate = dataa))
}
