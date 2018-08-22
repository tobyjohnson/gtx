#' Given an Ensembl ID, identify all significant QTLs
#'
#' This is an utility query to find all the QTLs for a specific gene.
#' It assumes that all the analyses for which entity_type is not null are QTLs
#'
#' @param geneId	Gene to search. Recommended to use Ensembl ID. The function will try to infer gene symbols, but result may not be correct.
#' @param cutoff	Only return QTLs below this p-value threshold
#'
#' @export
find_qtl_bygene <- function(geneId, cutoff=2e-5,
        dbc = getOption("gtx.dbConnection", NULL)) {
    gtxdbcheck(dbc)

    if (grepl("^ENS", geneId) ) {
	# geneId is already an Ensembl Id. Nothing more to do
	ensemblid = geneId
    } else {
	# geneId is not an Ensembl Id. Since all eQTLs use Ensembl Ids, this need to be converted.
        get_ensembl_query = paste0("SELECT ensemblid FROM genes WHERE hgncid = '", geneId, "' LIMIT 1")
#        print(get_ensembl_query)
        ensemblid = dbGetQuery(dbc, get_ensembl_query)
    }
#    print(ensemblid)
    if (!grepl("^ENS", ensemblid)) {
	stop(paste("Could not identify an Ensembl ID for", geneId, "- please specify an Ensembl ID manually"))
    }

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
