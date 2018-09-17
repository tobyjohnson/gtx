## Function to implement two-way fixed effects inverse variance weighted
## GWAS meta-analysis, on the fly

## FIXME would be more flexible to allow the option of different rsq_ge, maf_ge, emac_ge, case_emac_ge for each analysis

meta <- function(analysis1, analysis2,
                 entity, entity1, entity2,
                 rsq_ge, maf_ge, emac_ge, case_emac_ge, 
                 dbc = getOption("gtx.dbConnection", NULL)) {
  gtxdbcheck(dbc)

  ## Get metadata from table analyses for the two analyses
  amd <- sqlWrapper(getOption('gtx.dbConnection_cache_analyses', dbc),  
                    sprintf('SELECT analysis, ncase, ncontrol, ncohort, has_freq, has_rsq 
                             FROM analyses WHERE %s;',
                             gtx:::gtxwhat(analysis = c(analysis1, analysis2))),
                    uniq = FALSE)
  ## if ncohort is NA, replace with ncase+ncontrol, should not be needed with new schema definition
  amd <- within(amd, ncohort <- ifelse(!is.na(ncohort), ncohort, ncase + ncontrol))
  xn1 <- subset(amd, analysis == analysis1, select = 'ncohort', drop = TRUE) # FIXME throw error if NA
  xn2 <- subset(amd, analysis == analysis2, select = 'ncohort', drop = TRUE) # FIXME throw error if NA

  ## set up WHERE clauses for filtering
  xfilter1 <- gtx:::gtxfilter(rsq_ge = rsq_ge, maf_ge = maf_ge, emac_ge = emac_ge, case_emac_ge = case_emac_ge, analysis = analysis1)
  xfilter1 <- if (identical(xfilter1, '()')) '' else paste('AND', xfilter1) # FIXME should not be needed because gtxfilter should be fixed to return '(True)'
  xfilter2 <- gtx:::gtxfilter(rsq_ge = rsq_ge, maf_ge = maf_ge, emac_ge = emac_ge, case_emac_ge = case_emac_ge, analysis = analysis2)
  xfilter2 <- if (identical(xfilter2, '()')) '' else paste('AND', xfilter2)
  
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
                                 t1.beta AS beta1, t1.se AS se1, 
                                 t1.rsq AS rsq1, t1.freq AS freq1, 
                                 t2.beta AS beta2, t2.se AS se2, 
                                 t2.rsq AS rsq2, t2.freq AS freq2
                             FROM 
                                 (SELECT
                                      chrom, pos, ref, alt, beta, se, rsq, freq
                                  FROM %sgwas_results 
                                  WHERE %s %s AND %s AND pval IS NOT NULL
                                 ) AS t1 
                                 JOIN 
                                 (SELECT 
                                      chrom, pos, ref, alt, beta, se, rsq, freq
                                  FROM %sgwas_results 
                                  WHERE %s %s AND %s AND pval IS NOT NULL
                                 ) AS t2
                                 USING (chrom, pos, ref, alt);',
                            gtxanalysisdb(analysis1), 
                            gtxwhat(analysis1 = analysis1), # analysis1= argument allows only one analysis
                            xfilter1, 
                            xentity1$entity_where, 
                            gtxanalysisdb(analysis2), 
                            gtxwhat(analysis1 = analysis2), # analysis1= argument allows only one analysis
                            xfilter2, 
                            xentity2$entity_where 
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

  ## Calculate weighted average rsq and freq, weighted using ncohort(=ncase+ncontrol) from input metadata
  ## FIXME make robust if has_freq or has_rsq is 0 (=FALSE) for any input dataset
  res[ , rsq := (xn1*rsq1 + xn2*rsq2)/(xn1 + xn2)]
  res[ , freq := (xn1*freq1 + xn2*freq2)/(xn1 + xn2)]

  ## FIXME add option to order res by chrom/pos/ref/alt, make option since it will be expensive
  
  return(res)
  ## Note, should be options to return only the meta-analysis results (chrom,pos,ref,alt,beta,se,pval,rsq,freq)
  ##

  ## FIXME, since R cannot currently write back the full set of summary stats,
  ## we need a way to plug this into existing
  ## tools, e.g. a gwas() like summarization of top hits, Manhattan and QQ plots
  ## e.g. selection a region only and feed into regionplot(), coloc() etc.
  
}
