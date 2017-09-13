library(gtx)
library(testthat)

suppressPackageStartupMessages(library(gtx))
options(gtx.dbConnection = odbcConnect('impaladsn'))
sqlQuery(getOption('gtx.dbConnection'), 'USE ukbiobank;')



test_that("Testing GTEx analysis querying", 

 	lung.analyses = gtxanalyses(description_contains = c('lung'))
	lung.analyses.noliu17 = gtxanalyses(description_contains = c('lung'), analysis_not = 'liu17_lc')

	expect_that(nrow(lung.analyses) > 9)
	expect_that(nrow(lung.analysis.noliu17) < nrow(lung.analyses))
	)


	  
