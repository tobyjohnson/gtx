context("QTLs")

test_that("Zhernakova analysis entry exists", {
  expect_equal(nrow(gtxanalyses("eqtl_gene_blood_zhernakova")), 1)
})
	

test_that("IPSDSN analysis entry exists ", {
  expect_equal(nrow(gtxanalyses("eqtl_gene_ipsdsn")), 1)
})

test_that("IPSDSN coloc PTEN", {
  skip('Not working on gtx2_dev.  Test may need modification.')
  
	res = coloc('eqtl_gene_blood_zhernakova', 
		       'eqtl_gene_ipsdsn', chrom = '10', pos = 87863113, surround = 400e3, 
		       entity = 'ENSG00000171862', style="none")
	expect_gt(res$nvariants, 10) # should be 18
})

test_that("IPSDSN coloc MGAT3", {
  skip('Not working on gtx2_dev.  Test may need modification.')
  
	res = coloc('eqtl_gene_blood_zhernakova', 
		     'eqtl_gene_ipsdsn', chrom = '22', pos = 39883229, surround = 400e3, 
		     entity = 'MGAT3', style="none")
	expect_gt(res$nvariants, 0) 
})


test_that("find_qtl_bygene function with Ensembl id", {
	geneId <- "ENSG00000137033"
	gene_eQTLs <- find_qtl_bygene(geneId) 
	expect_gt(nrow(gene_eQTLs), 0)
	expect_gt(length(unique(gene_eQTLs$analysis)), 0)

})
test_that("find_qtl_bygene function with Gene Symbol", {

	geneId <- "IL33"
	gene_eQTLs <- find_qtl_bygene(geneId) 
	expect_gt(nrow(gene_eQTLs), 0)
	expect_gt(length(unique(gene_eQTLs$analysis)), 0)
})
test_that("find_qtl_bygene function FAILS with wrong gene symbol", {
	geneId <- "das asd as das das"
	expect_error(find_qtl_bygene(geneId) )
})


