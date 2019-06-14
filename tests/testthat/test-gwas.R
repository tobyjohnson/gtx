context("gwas")

test_that("Unknown error", {
  options('gtx.debug'=TRUE)
  results <- g1 <- gwas('UKB50k_TA_copdExacerbation_NI_HES_SPIROMETRY',
                        style = 'manhattan', pval_thresh = 5e-06)
})