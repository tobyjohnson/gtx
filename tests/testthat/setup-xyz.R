# This file is executed before all tests start
# See the file teardown-xyz, which be executed after the last test, and will close these connections

options(gtx.dbConnection = dbConnect(odbc(), dsn = 'impaladsn'))
dbGetQuery(getOption('gtx.dbConnection'), 'USE ukbiobank;')