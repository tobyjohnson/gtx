context("mranalysis")
library(DescTools)

test_that("Non-empty output if searching for one word", {
  res <- extract_study_info("asthma")
  expect_gt(nrow(res),0)
})

test_that("Non-empty output search by chrom:pos", {
  #Extract study meta-information
  study <- extract_study_info("asthma")
  
  #Extract instruments
  res <- extract_exposure(analyses = study)
  expect_gt(nrow(res),0)
  
  #Extract instruments with rsIDs
  res_rsid <- extract_exposure(analyses = study, rsid = TRUE)
  expect_equal(ncol(res_rsid)-ncol(res),1)
})

