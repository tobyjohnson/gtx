context("regionplot")

test_that("basic regionplot works", {
  skip('Not working on gtx2_dev.  Test may need modification.')
  
  expect_success(regionplot("gsk_ukb500k_imputed_interim_bolt_f_3571_0_0_f_BIN_Back_pain_for_3_months_vs_no_pain", chrom = '6', pos = 32566577, surround = 5e3, ensemblid = 'ENSG00000013573'))
})