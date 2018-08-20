#' Given an Ensembl ID, identify all significant QTLs
#'
#' 
find_qtl_bygene <- function(ensemblid, cutoff=2e-5,
        dbc = getOption("gtx.dbConnection", NULL)) {
    gtxdbcheck(dbc)

    qtl_analyses_sql = 'SELECT analysis FROM analyses WHERE entity_type NOT IN ("", "NA")'
    qtl_analyses = dbGetQuery(dbc, qtl_analyses_sql)

    qtls_sql = paste0("SELECT * FROM gwas_results WHERE entity='", 
                      ensemblid, 
                      "' AND analysis in ('", 
                      "", paste(qtl_analyses$analysis, collapse="', '"),
                      "') AND pval<", cutoff)
#    print(qtls_sql)

    qtls = dbGetQuery(dbc, qtls_sql)
    qtls
}
#options(gtx.dbConnection = dbConnect(odbc(), dsn = 'impaladsn'))
#dbGetQuery(getOption('gtx.dbConnection'), 'USE ukbiobank;')
#find_qtl_bygene("ENSG00000137033")
