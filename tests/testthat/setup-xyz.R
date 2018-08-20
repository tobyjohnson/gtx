# This file is executed before all tests start
# See the file teardown-xyz, which be executed after the last test, and will close these connections
library(futile.logger)

flog.appender(appender.file("unit_tests.log"))

options(gtx.dbConnection = dbConnect(odbc(), dsn = 'impaladsn'))
dbGetQuery(getOption('gtx.dbConnection'), 'USE ukbiobank;')