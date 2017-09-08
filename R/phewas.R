phewas.data <- function(chrom, pos, rs,
                        analysis, analysis_not, 
                        description_contains,
                        phenotype_contains,
                        ncase_ge,
                        ncohort_ge,
                        dbc = getOption("gtx.dbConnection", NULL)) {
    gtxdbcheck(dbc)

    a1 <- gtxanalyses(analysis = analysis, analysis_not = analysis_not, 
                      description_contains = description_contains,
                      phenotype_contains = phenotype_contains,
                      ncase_ge = ncase_ge, ncohort_ge = ncohort_ge,
                      has_access_only = TRUE, 
                      dbc = dbc)
    v1 <- sqlWrapper(dbc,
                     sprintf('SELECT chrom, pos, ref, alt FROM sites_ukb WHERE %s;',
                             gtxwhere(chrom = chrom, pos = pos, rs = rs))) # default uniq = TRUE
    res <- do.call(rbind, lapply(unique(a1$results_db), function(results_db) {
        sqlWrapper(getOption('gtx.dbConnection'), 
                   sprintf('SELECT analysis, feature, beta, se, pval, rsq, freq FROM %s.gwas_results WHERE %s AND chrom=\'%s\' AND pos=%s AND ref=\'%s\' AND alt=\'%s\';',
                           sanitize(results_db, type = 'alphanum'),
                           gtxwhat(analysis = a1$analysis),
                           sanitize1(v1$chrom, values = c(as.character(1:22), "X", "Y")),
                           sanitize1(v1$pos, type = "int"),
                           sanitize1(v1$ref, type = "ACGT+"),
                           sanitize1(v1$alt, type = "ACGT+")),
                   uniq = FALSE)
    }))
    res <- res[order(res$pval), ]
    return(res)
}
