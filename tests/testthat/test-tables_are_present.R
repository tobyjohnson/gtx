context("tables")

test_that("impala ukbiobank.analyses works", {
  skip('Not working on gtx2_dev.  Test may need modification.')
  
  expect_equal(nrow(sqlQuery(getOption('gtx.dbConnection'), "SELECT * FROM ukbiobank.analyses LIMIT 10")), 10)
})

test_that("impala ukbiobank.gwas_results works", {
  skip('Not working on gtx2_dev.  Test may need modification.')
  
  expect_equal(nrow(sqlQuery(getOption('gtx.dbConnection'), "SELECT * FROM ukbiobank.gwas_results LIMIT 10")), 10)
})
