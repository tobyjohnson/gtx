#' Given an Ensembl ID, identify all significant QTLs
#'
#' This is an utility query to find all the QTLs for a specific gene.
#' It assumes that all the analyses for which entity_type is not null are QTLs
#'
#' @param ensemblid	Ensembl Id to search
#' @param cutoff	Only return QTLs below this p-value threshold
#'
#' @export
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
