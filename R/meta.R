## Function to implement two-way fixed effects inverse variance weighted
## GWAS meta-analysis, on the fly

meta <- function(analysis1, analysis2,
                 entity, entity1, entity2, 
                 dbc = getOption("gtx.dbConnection", NULL)) {
  gtxdbcheck(dbc)
  
  ## substitute generic entity for entity1 and entity2 if needed
  if (missing(entity1) && !missing(entity)) entity1 <- entity
  if (missing(entity2) && !missing(entity)) entity2 <- entity
  
  ## Determine entity, if required, for each analysis
  xentity1 <- gtxentity(analysis1, entity = entity1)
  xentity2 <- gtxentity(analysis2, entity = entity2)
  
  ## Get association statistics, inner JOIN (silently drop variants not present in both datasets...)
  ## consider options for FULL, LEFT and RIGHT JOINs
  res <- sqlWrapper(dbc, 
                    sprintf('SELECT 
                                 t1.chrom, t1.pos, t1.ref, t1.alt,
                                 t1.beta AS beta1, t1.se AS se1 
                                 t2.beta AS beta2, t2.se AS se2 
                             FROM 
                                 (SELECT
                                      chrom, pos, ref, alt, beta, se, 
                                  FROM %s.gwas_results 
                                  WHERE %s %s
                                 ) AS t1 
                                 JOIN 
                                 (SELECT 
                                      chrom, pos, ref, alt, beta, se, 
                                  FROM %s.gwas_results 
                                  WHERE %s %s
                                 ) AS t2
                                 USING (chrom, pos, ref, alt);',
                            gtxanalysisdb(analysis1), 
                            gtxwhat(analysis1 = analysis1), # analysis1= argument allows only one analysis
                            if (!is.null(xentity1)) sprintf(' AND feature=\'%s\'', xentity1$entity) else '', # FIXME will change to entity
                            gtxanalysisdb(analysis2), 
                            gtxwhat(analysis1 = analysis2), # analysis1= argument allows only one analysis
                            if (!is.null(xentity2)) sprintf(' AND feature=\'%s\'', xentity2$entity) else ''  # FIXME will change to entity
                            ),
                    uniq = FALSE) # expect >=1 rows
  res <- data.table(res)

  ## Compute precision of meta-analysis result (=1/Variance)
  res[ , prec := se1^-2 + se2^-2] # a sum over input GWAS, can be kept as a running total with vector add for each input GWAS in turn
  ## Compute meta-analysis beta, as inverse variance weighted average of input betas
  res[ , beta := (beta1*se1^-2 + beta2*se2^-2)/prec] # numerator is a sum over input GWAS, running total as above, then divide by prec
  ## Compute meta-analysis standard error
  res[ , se := prec^-0.5]
  ## Compute meta-analysis P-value
  res[ , pval := pnorm(-abs(beta)/se)*2]
  ## Drop prec column
  res[ , prec := NULL]

  ## FIXME add calculation of weighted average rsq and freq, weighted using ncohort(=ncase+ncontrol) from input metadata
  
  return(res)
  ## Note, should be options to return only the meta-analysis results (chrom,pos,ref,alt,beta,se,pval,rsq,freq)
  ##
}
