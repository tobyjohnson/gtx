context("QTLs")
setup({
	options(gtx.dbConnection = odbcConnect('impaladsn'))
	sqlQuery(getOption('gtx.dbConnection'), 'USE ukbiobank;')
}
)
teardown(
	odbcCloseAll()
)

test_that("Zhernakova analysis entry exists", {
  expect_equal(nrow(gtxanalyses("eqtl_gene_blood_zhernakova")), 1)
})
	

test_that("IPSDSN analysis entry exists ", {
  expect_equal(nrow(gtxanalyses("eqtl_blood_ipsdsn")), 1)
})
		
