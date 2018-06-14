options(gtx.dbConnection = odbcConnect('impaladsn'))
sqlQuery(getOption('gtx.dbConnection'), 'USE ukbiobank;')

