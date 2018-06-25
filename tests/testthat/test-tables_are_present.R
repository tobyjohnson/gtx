context("tables")

test_that("impala ukbiobank.analyses works", {
  expect_equal(nrow(sqlQuery(getOption('gtx.dbConnection'), "SELECT * FROM ukbiobank.analyses LIMIT 10")), 10)
})

test_that("impala ukbiobank.gwas_results works", {
  expect_equal(nrow(sqlQuery(getOption('gtx.dbConnection'), "SELECT * FROM ukbiobank.gwas_results LIMIT 10")), 10)
})
