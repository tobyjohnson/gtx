#' Given an Ensembl ID, identify all significant QTLs
#'
#' 
find_qtl_bygene <- function(ensemblid, cutoff=2e-5
			    dbc = getOption("gtx.dbConnection", NULL)) {
	gtxdbcheck(dbc)
	# TODO: Could call gtxentity to identify the entity type

	qtl_analyses_sql = 'SELECT analysis FROM analyses WHERE entity_type NOT IN ("", "NA")'
	qtl_analyses = sqlWrapper(dbc, qtl_analyses_sql)

	qtls_sql = sprintf("SELECT * FROM gwas_results WHERE entity=%s AND analysis in ('%s') AND pval<%s", ensemblid, paste(qtl_analyses, collapse="','"), cutoff)
	print(qtls_sql)

	qtls = sqlWrapper(dbc, qtls_sql)
	qtls
}

