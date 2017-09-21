library(gtx)
library(testthat)

suppressPackageStartupMessages(library(gtx))
options(gtx.dbConnection = odbcConnect('impaladsn'))
sqlQuery(getOption('gtx.dbConnection'), 'USE ukbiobank;')



test_that("Testing GTEx analysis querying", {

 	lung.analyses = gtxanalyses(description_contains = c('lung'))

	lung.analyses.noliu17 = gtxanalyses(description_contains = c('lung'), analysis_not = 'liu17_lc')

	expect_equal(nrow(lung.analyses)>9, T)
	expect_equal(nrow(lung.analyses.noliu17) < nrow(lung.analyses), T)
})

test_that("Phewas Querying", {
	  phewas1 <- phewas.data(rs = 'rs11001819', 
		 description_contains = c('lung', 'copd', 'fev1'), analysis_not = 'liu17_lc')

	  expect_that(nrow(phewas1), is_more_than(2))
	  
}	  )
