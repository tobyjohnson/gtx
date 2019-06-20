##########
#phewas
getSQLArgsPhewas00 <- function(chrom, pos, ref, alt, rs, dbc, uniq, zrok){
  sqlText <- getSQLQueryPhewas00(chrom, pos, ref, alt, rs)
  allArgs <- list(dbc = dbc,
                  sql = sqlText,
                  uniqu = uniq,
                  zrok = zrok)
  
  return(allArgs)
}

getSQLQueryPhewas00 <- function(chrom, pos, ref, alt, rs){
  queryText <- sprintf('SELECT chrom, pos, ref, alt, rsid AS rs FROM sites WHERE %s;',
                       gtxwhere(chrom = chrom, pos = pos, ref = ref, alt = alt, rs = rs))
  
  return(queryText)
}


##
