context("coloc")

setup({
	options(gtx.dbConnection = odbcConnect('impaladsn'))
	sqlQuery(getOption('gtx.dbConnection'), 'USE ukbiobank;')
}
)
teardown(
	odbcCloseAll()
)

test_that("test works", {
  expect_equal(2 * 2, 4)
})
