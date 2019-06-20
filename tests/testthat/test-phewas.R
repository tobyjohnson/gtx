context("phewas")

# NOTE: Cannot guarantee exact count, as more data will be loaded
#       Will improve test when GTX supports tests against static data
test_that("handles missing p-values", {
  expect_gte(nrow(phewas(rs='rs144979264')), 5878)
})
	
