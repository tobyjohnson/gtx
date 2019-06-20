context("regionplot")

test_that("basic regionplot works", {
  skip('Not working on gtx2_dev.  Test may need modification.')
  
  expect_success(regionplot("gsk_ukb500k_imputed_interim_bolt_f_3571_0_0_f_BIN_Back_pain_for_3_months_vs_no_pain", chrom = '6', pos = 32566577, surround = 5e3, ensemblid = 'ENSG00000013573'))
})

test_that("UB-131 Handle no LD data when selecting index variant by smallest pval", {
  flog.info('UB-131')
  regionplot('GSK500kV2b_TA_copdMapLoose_noAsthma_primAndSec_a1', hgncid = 'IL5', case_emac_ge = 25, style = 'ld')
})

test_that("UB-135 Handle no LD data when selecting index variant by RSID", {
  flog.info('UB-135')
  regionplot('GSK500kV3_Eosinophill_percentage', rs = 'rs575996679', case_emac_ge = 25, style = 'ld')
})

test_that("UB-138 Confirm expected posterior probabilities & proper credible set annotation", {
  input    <- read_csv(file = system.file("testdata", package = "gtx", "UB-138_input.csv.gz"))
  expected <- read_csv(file = system.file("testdata", package = "gtx", "UB-138_expected.csv.gz"))
  actual   <- fm_signal(input, priorsd = 1, priorc = 1e-5, cs_size = 0.95, cs_only = FALSE)
  
  expect_identical(dim(actual), dim(expected))
  expect_equal(actual, expected) # 
}
)

test_that("UB-66", {
  flog.info('UB-66')
  regionplot(analysis = 'GSK500kV2b_HES_Bronchiectasis_a1', chrom = '6', pos = 32565465, maf_ge = 0.05, rsq_ge = 0.1, style = 'ld')
})

test_that("UB-107", {
  flog.info('UB-107')
  regionplot('GSK500kV2b_TA_copdMapLoose_primAndSec_a1', 
             chrom = 10, pos_start = 8840700, pos_end = 9264716, style = 'classic') 
})
