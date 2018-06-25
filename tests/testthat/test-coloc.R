context("coloc")

test_that("coloc backpain", {

	res = coloc('gsk_ukb500k_imputed_interim_bolt_f_3571_0_0_f_BIN_Back_pain_for_3_months_vs_no_pain', 
		       'GSK500kV2b_Back_pain_3_mos_vs_no_pain_a2', chrom = '6', pos = 32566577, surround = 400e3, 
		       entity = 'ENSG00000198502', style="none")
	expect_gt(res$nvariants, 10) # should be 18
})

