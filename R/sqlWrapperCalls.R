##########
#phewas
phewasGetV1 <- function(chrom, pos, ref, alt, rs, dbc,
                        uniq, zrok, connectionType) {
  
  connectionArguments <- getSQLArgsPhewas00(chrom, pos, ref, alt, rs, dbc,
                                            uniq, zrok)
  
  v1Data <- getDataFromDB(connectionType = connectionType, connectionArguments = connectionArguments)
  
  
  
  return(v1Data)
  
}

getSQLArgsPhewas00 <-
  function(chrom, pos, ref, alt, rs, dbc, uniq, zrok) {
    sqlText <- getSQLQueryPhewas00(chrom, pos, ref, alt, rs)
    allArgs <- list(
      dbc = dbc,
      sql = sqlText,
      uniq = uniq,
      zrok = zrok
    )
    
    return(allArgs)
  }

getSQLQueryPhewas00 <- function(chrom, pos, ref, alt, rs) {
  queryText <-
    sprintf(
      'SELECT chrom, pos, ref, alt, rsid AS rs FROM sites WHERE %s;',
      gtxwhere(
        chrom = chrom,
        pos = pos,
        ref = ref,
        alt = alt,
        rs = rs
      )
    )
  
  return(queryText)
}


##############################
#res nearby
getResLoop <- function(results_db, v1, nearby, uniq, zrok, connectionType) {
  connectionArguments <- getSQLArgsPhewas01(results_db, v1, nearby, uniq, zrok)
  
  loopRes <- getDataFromDB(connectionType = connectionType, connectionArguments)
  
  return(loopRes)
}

getSQLArgsPhewas01 <- function(results_db, v1, nearby, uniq, zrok) {
  sqlText <- getSQLQueryPhewas01(results_db, v1, nearby)
  allArgs <- list(
    dbc = getOption('gtx.dbConnection'),
    sql = sqlText,
    uniq = uniq,
    zrok = zrok
  )
  
  return(allArgs)
}

getSQLQueryPhewas01 <- function(results_db, v1, nearby) {
  queryText <- sprintf(
    'SELECT analysis, entity, min(pval) AS pval_nearby, count(pval) AS num_pvals \
    FROM %sphewas_results \
    WHERE %s AND pval IS NOT NULL GROUP BY analysis, entity;',
    results_db,
    gtxwhere(
      chrom = v1$chrom,
      pos_ge = v1$pos - nearby,
      pos_le = v1$pos + nearby
    )
  )
  
  return(queryText)
}



####################
#res

getRes <- function(results_db, v1, uniq, zrok, connectionType) {
  connectionArguments <-
    getSQLArgsPhewas02(results_db, v1, uniq, zrok)
  
  loopRes <-
    getDataFromDB(connectionType = 'connectionType', connectionArguments)
  
  return(loopRes)
}


getSQLArgsPhewas02 <- function(results_db, v1, uniq, zrok) {
  sqlText <- getSQLQueryPhewas02(results_db, v1)
  allArgs <- list(
    dbc = getOption('gtx.dbConnection'),
    sql = sqlText,
    uniq = uniq,
    zrok = zrok
  )
  
  return(allArgs)
}

getSQLQueryPhewas02 <- function(results_db, v1) {
  queryText <- sprintf('SELECT analysis, entity, beta, se, pval, rsq, freq \
                                FROM %sphewas_results \
                                WHERE %s AND pval IS NOT NULL;',
                       results_db,
                       do.call(gtxwhere, v1[ , c('chrom', 'pos', 'ref', 'alt')]))
  
  return(queryText)
}


####################
#res in other case
getResCase2 <- function(results_db, v1, w1, uniq, zrok, connectionType){
  connectionArguments <- getSQLArgsPhewas03(results_db, v1, w1, uniq, zrok)
  
  loopRes <- getDataFromDB(connectionType = connectionType, connectionArguments)
  
  return(loopRes)
}

getSQLArgsPhewas03 <- function(results_db, v1, w1, uniq, zrok) {
  sqlText <- getSQLQueryPhewas03(results_db, v1, w1)
  allArgs <- list(
    dbc = getOption('gtx.dbConnection'),
    sql = sqlText,
    uniq = uniq,
    zrok = zrok
  )
  
  return(allArgs)
}

getSQLQueryPhewas03 <- function(results_db, v1) {
  queryText <- sprintf(
    'SELECT phewas_results.analysis, entity, beta, se, pval, rsq, freq \
              FROM %sphewas_results LEFT JOIN analyses_tags USING (analysis) \
              WHERE %s AND pval IS NOT NULL AND %s;',
    results_db,
    do.call(gtxwhere, v1[, c('chrom', 'pos', 'ref', 'alt')]),
    w1
  )
  
  return(queryText)
}


#######################
#res_nearby other case

getResCase3 <- function(results_db, v1, nearby, w1, uniq, zrok, connectionType){
  
  connectionArguments <- getSQLArgsPhewas04(results_db, v1, nearby, w1, uniq, zrok)
  
  loopRes <- getDataFromDB(connectionType = connectionType, connectionArguments)
  
  return(loopRes)
}

getSQLArgsPhewas04 <- function(results_db, v1, nearby, w1, uniq, zrok){
  sqlText <- getSQLQueryPhewas04(results_db, v1, nearby, w1)
  allArgs <- list(
    dbc = getOption('gtx.dbConnection'),
    sql = sqlText,
    uniq = uniq,
    zrok = zrok
  )
  
  return(allArgs)
}

getSQLQueryPhewas04 <- function(results_db, v1, nearby, w1) {
  queryText <- sprintf(
    'SELECT analysis, entity, min(pval) AS pval_nearby, count(pval) AS num_pvals \
                FROM %sphewas_results \
                WHERE %s AND pval IS NOT NULL AND %s \
                GROUP BY analysis, entity;',
    results_db,
    gtxwhere(
      chrom = v1$chrom,
      pos_ge = v1$pos - nearby,
      pos_le = v1$pos + nearby
    ),
    w1
  )
  
  return(queryText)
}
