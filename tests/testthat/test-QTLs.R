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
  expect_equal(nrow(gtxanalyses("eqtl_gene_ipsdsn")), 1)
})

test_that("IPSDSN coloc  HLA-DRB5", {

	res = coloc('gsk_ukb500k_imputed_interim_bolt_f_3571_0_0_f_BIN_Back_pain_for_3_months_vs_no_pain', 
		       'eqtl_gene_ipsdsn', chrom = '6', pos = 32566577, surround = 400e3, 
		       entity = 'ENSG00000198502', style="none")
	expect_gt(res$nvariants, 10) # should be 18
})

test_that("IPSDSN coloc SOX5", {
	res = coloc('gsk_ukb500k_imputed_interim_bolt_f_3571_0_0_f_BIN_Back_pain_for_3_months_vs_no_pain', 
		     'eqtl_gene_ipsdsn', chrom = '6', pos = 32566577, surround = 400e3, 
		     entity = 'SOX5', style="none")
	expect_gt(res$nvariants, 0) 
})
