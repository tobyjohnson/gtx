context("my-test")

test_that("impala analyses works", {
  expect_equal(nrow(sqlQuery(getOption('gtx.dbConnection'), "SELECT * FROM ukbiobank.analyses LIMIT 10")), 10)
})
