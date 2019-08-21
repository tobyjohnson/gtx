#' Extracts top instruments based on the chosen p-value (the instruments are distance-pruned)
#' 
#' This function allows the users to pre-process and extract the 
#' required data from the GWA studies database, which then
#' can be used to run Mendelian Randomization and other related
#' analyses.
#' 
#' @author Valeriia Haberland \email{valeriia.haberland@@bristol.ac.uk}
#' @param p (optional): p-value threshold for the top findings where 5e-08 (default). Note! the instruments with a more lenient p-value threshold cannot be extracted. 
#' @param analyses (optional): the study information: at least should have a column with analysis IDs; entity_type IDs if applicable.
#' @param str (optional): a keyword/keyphrase to choose the studies of interest; by default searches all studies
#' @param rsid (optional): logical; if FALSE (default) returns only chrom:position; if TRUE includes rsIDs (Warning! SNPs without rsIDs may not be returned); 
#' @param dbc (default): other parameters
#' @return the dataframe with exposure instruments (Warning! they have to be formatted and possibly clumped to be used in the MR analysis)
#'
#' @export
#'
extract_exposure <- function(p=5e-08,analyses=NULL,str='',rsid=FALSE,dbc=getOption("gtx.dbConnection", NULL)){

if(rsid==FALSE){
  
  if(is.null(analyses)){
  
  res_exp<-NULL  
  ## Extract the exposure instruments given the p-value threshold
  ## and search keywords (optional)
  ## Return all (default) or chosen GWASs based on the keywords
  res_exp <- odbc::dbGetQuery(dbc, paste0("SELECT g.pos_index,g.chrom,g.ref_index,g.alt_index,g.beta_index,g.se_index,
                                             g.pval_index,g.freq_index,g.analysis,g.entity,a.description,a.phenotype,
                                             a.unit,a.ncase,a.ncontrol,a.ncohort FROM gwas_results_top_hits_with_entity AS g
                                           INNER JOIN analyses AS a ON g.analysis=a.analysis
                                           WHERE g.pval_index < ", p," AND a.description LIKE '%", str, "%'
                                           ORDER BY g.analysis,g.entity;"))
  
  if( (is.null(res_exp)) | (dim(res_exp)[1] == 0) ) message("Warning! No result has been returned by the query.")
  
  }else{
    
  res_exp<-NULL
  ## Extract the exposure instruments given the p-value threshold
  ## and the list of analysis IDs (required) and entity IDs (if applicable)
  ## Return chosen GWASs based on the analysis IDs
  if(is.data.frame(analyses)){
 
  #Create a string of analyses
  s<-paste(unique(analyses$analysis), collapse="','")
  s<-paste0("'",s,"'")
  
  res_exp<-odbc::dbGetQuery(dbc, paste0("SELECT g.pos_index,g.chrom,g.ref_index,g.alt_index,g.beta_index,
                                          g.se_index,g.pval_index,g.freq_index,g.analysis,g.entity,a.description,a.phenotype,
                                          a.unit,a.ncase,a.ncontrol,a.ncohort FROM gwas_results_top_hits_with_entity AS g
                                         INNER JOIN analyses AS a ON g.analysis=a.analysis
                                         WHERE g.pval_index < ", p," AND g.analysis IN (",s,")
                                         ORDER BY g.analysis,g.entity;"))

  #Choose only the right <analysis,entity> pairs
  if( (!is.null(res_exp)) & (dim(res_exp)[1] > 0) ){
    
     ids<-unique(analyses[,which(names(analyses) %in% c("analysis","entity_type","entity"))])
     if(is.data.frame(ids)){
       colnames(ids)[2]<-"entity"
       res_exp<-dplyr::left_join(res_exp,ids)
     }
     else{
       message("Warning! Entity cannot be identified; 'analyses' dataframe is probably missing 'entity' column.")
     }
   }
   else{
     message("Warning! No result has been returned by the query.")
   }
  }
  else{
    stop("Error: 'analyses' [data.frame] should contain at least one column with analysis IDs if searching for specific studies.")
    }
  }
}
  else{
    
    message("Warning! SNPs without rsIDs may not be returned.")
    
    if(is.null(analyses)){
      
      res_exp<-NULL  
      ## Extract the exposure instruments given the p-value threshold
      ## and search keywords (optional)
      ## Return all (default) or chosen GWASs based on the keywords
      res_exp <- odbc::dbGetQuery(dbc, paste0("SELECT g.pos_index,g.chrom,s.rsid,g.ref_index,g.alt_index,g.beta_index,g.se_index,
                                                g.pval_index,g.freq_index,g.analysis,g.entity,a.description,a.phenotype,
                                                a.unit,a.ncase,a.ncontrol,a.ncohort FROM gwas_results_top_hits_with_entity AS g
                                               INNER JOIN analyses AS a ON g.analysis=a.analysis
                                               INNER JOIN sites AS s ON g.pos_index=s.pos AND g.chrom=s.chrom 
                                                        AND g.ref_index=s.ref AND g.alt_index=s.alt 
                                               WHERE g.pval_index < ", p," AND a.description LIKE '%", str, "%'
                                               ORDER BY g.analysis,g.entity;"))
      
     if( (is.null(res_exp)) | (dim(res_exp)[1] == 0) ) message("Warning! No result has been returned by the query.")
        
    }else{
      
      res_exp<-NULL
      ## Extract the exposure instruments given the p-value threshold
      ## and the list of analysis IDs (required) and entity IDs (if applicable)
      ## Return chosen GWASs based on the analysis IDs
      if(is.data.frame(analyses)){
        
        #Create a string of analyses
        s<-paste(unique(analyses$analysis), collapse="','")
        s<-paste0("'",s,"'")
        
        res_exp <- odbc::dbGetQuery(dbc, paste0("SELECT g.pos_index,g.chrom,s.rsid,g.ref_index,g.alt_index,g.beta_index,
                                                  g.se_index,g.pval_index,g.freq_index,g.analysis,g.entity,a.description,a.phenotype,
                                                  a.unit,a.ncase,a.ncontrol,a.ncohort FROM gwas_results_top_hits_with_entity AS g
                                                 INNER JOIN analyses AS a ON g.analysis=a.analysis
                                                 INNER JOIN sites AS s ON g.pos_index=s.pos AND g.chrom=s.chrom 
                                                       AND g.ref_index=s.ref AND g.alt_index=s.alt 
                                                 WHERE g.pval_index < ", p," AND g.analysis IN (",s,")
                                                 ORDER BY g.analysis,g.entity;"))
        
        #Choose only the right <analysis,entity> pairs
        if( (!is.null(res_exp)) & (dim(res_exp)[1] > 0) ){
          
          ids<-unique(analyses[,which(names(analyses) %in% c("analysis","entity_type","entity"))])
          if(is.data.frame(ids)) {
            colnames(ids)[2]<-"entity"
            res_exp<-dplyr::left_join(res_exp,ids)
          }
          else{
            message("Warning! Entity cannot be identified; 'analyses' dataframe is probably missing 'entity' column.")
          }
        }
        else{
          message("Warning! No result has been returned by the query.")
        }
      }
      else{
        stop("Error: 'analyses' [data.frame] should contain at least one column with analysis IDs where fuzzy=FALSE.")
      }
    }
  }
  
  return(res_exp)
}

#' Extracts the outcome instruments based on the exposure SNPs
#' 
#' This function allows the users to pre-process and extract the 
#' required data from the GWA studies database, which then
#' can be used to run Mendelian Randomization and other related
#' analyses.
#' 
#' @author Valeriia Haberland \email{valeriia.haberland@@bristol.ac.uk}
#' @param snp_list (required if the exposure instruments extracted from the custom file/existing dataframe): the dataframe which at least should have either two columns position (1) and chromosome (2)
#' or one column with rsIDs.
#' @param analyses_exp (required if the exposure instruments extracted from the database): the exposure study information; should have at least one column: analysis and if applicable: entity. 
#' @param analyses_out (required): the outcome study information; should have at least one column: analysis and if applicable: entity.
#' @param p (optional): the p-value threshold to extract the exposure instruments from the database (if analyses_exp is not NULL).
#' @param rsid (optional): logical; if FALSE (default) returns SNPs identified only by chrom:position; if TRUE includes rsIDs (Warning! SNPs without rsIDs may not be returned); 
#' @param search (optional): if "pos:chrom" (default) searches the instruments in the database by their chromosome position; alternatively, searches by their rsIDs (rsid should be "TRUE")
#' (only applicable to the searches for external files with the exposure instruments, i.e. snp_list is not NULL).  
#' @param dbc (default): other parameters 
#' @return the dataframe with outcome instruments
#' 
#' @export
#' 
#' @import DescTools
#' 
extract_outcome<-function(snp_list=NULL,analyses_exp=NULL,analyses_out,p=5e-08,rsid=FALSE,search="pos:chrom",dbc=getOption("gtx.dbConnection", NULL)){

message("Warning! Ensure that DescTools library is attached.")
  
if(rsid==FALSE){
  
  ## Extract the data based on the input 
  ## from the custom file
  if( is.null(analyses_exp) ) {
    
   if( !is.null(snp_list) & is.data.frame(snp_list) ){
    
    #Create a string of analyses
     
    s<-paste(unique(analyses_out$analysis), collapse="','")
    s<-paste0("'",s,"'")
    
    #Ensure the right data.type
    snp_list[1]<-as.integer(unlist(snp_list[1]))
    snp_list[2]<-as.character(unlist(snp_list[2]))
    
    #Choose only the unique instruments
    snp_list<-unique(snp_list[,c(1,2)])
    
    if(ncol(snp_list)>1){
      ## Extract the outcome instruments using SNPs,
      ## obtained from the custom file/dataframe
      ## Query needs an optimization (in development)
      res <- apply(snp_list,1,function(X) odbc::dbGetQuery(dbc, paste0("SELECT g.pos,g.chrom,g.ref,g.alt,g.beta,g.se,
                                                                         g.pval,g.freq,g.analysis,g.entity,a.description,a.phenotype,
                                                                         a.unit,a.ncase,a.ncontrol,a.ncohort FROM gwas_results as g
                                                                        INNER JOIN analyses as a ON g.analysis=a.analysis
                                                                        WHERE g.pos=", X[1]," AND g.chrom='", X[2],"' AND g.analysis IN (",s,")
                                                                        ORDER BY g.analysis,g.entity;")))
      
      #Combine all results
      res_out<-NULL
        for(i in 1:length(res))
          res_out<-rbind(res_out,res[[i]])
    
      #Throw a warning if the result is empty
      #Choose only the right <analysis,entity> pairs
      if( (!is.null(res_out)) & (dim(res_out)[1] > 0) ){
      
        ids<-unique(analyses_out[,which(names(analyses_out) %in% c("analysis","entity_type","entity"))])
      
        if(is.data.frame(ids)){
          colnames(ids)[2]<-"entity"
          res_out<-dplyr::left_join(res_out,ids)
        
        }else{ message("Warning! Entity cannot be identified; 'analyses_out' dataframe is probably missing 'entity' column.") }
      }else{ message("Warning! No result has been returned by the query.") }
    }else{ stop("Error: snp_list should have at least two columns in this order: position (1) and chromosome (2).") }
  }else{ stop("Error: snp_list [data.frame] cannot be empty and should have at least two columns in this order: position (1) and chromosome (2).") }
}else{
    
    if(!is.null(snp_list)) message("Warning! These instruments will not be considered; analyses_exp should be NULL.")
    
    res_out<-NULL
    
    if(is.data.frame(analyses_exp)){  
      
      #Create a string of the exposure analyses
      s_exp<-paste(unique(analyses_exp$analysis), collapse="','")
      s_exp<-paste0("'",s_exp,"'")
      
      #Create a string of the outcome analyses
      s_out<-paste(unique(analyses_out$analysis), collapse="','")
      s_out<-paste0("'",s_out,"'")
      
      ## Extract the outcome instruments using SNPs
      res_out <- odbc::dbGetQuery(dbc, paste0("SELECT g.pos,g.chrom,g.ref,g.alt,g.beta,g.se,
                                                g.pval,g.freq,g.analysis,g.entity,a.description,a.phenotype,
                                                a.unit,a.ncase,a.ncontrol,a.ncohort FROM gwas_results as g
                                              INNER JOIN (SELECT DISTINCT gg.pos_index,gg.chrom FROM gwas_results_top_hits_with_entity AS gg
                                                      WHERE gg.pval_index < ", p," AND gg.analysis IN (",s_exp,")) AS out 
                                              ON g.pos=out.pos_index AND g.chrom=out.chrom
                                              INNER JOIN analyses AS a ON g.analysis=a.analysis
                                              WHERE g.analysis IN (",s_out,")
                                              ORDER BY g.analysis,g.entity;"))
      
      #Choose only the right <analysis,entity> pairs
      if( (!is.null(res_out)) & (dim(res_out)[1] > 0) ){
          
          ids<-unique(analyses_out[,which(names(analyses_out) %in% c("analysis","entity_type","entity"))])
          if(is.data.frame(ids)){
            colnames(ids)[2]<-"entity"
            res_out<-dplyr::left_join(res_out,ids)
            
          }else{ message("Warning! Entity cannot be identified; 'analyses_out' dataframe is probably missing 'entity' column.") }
        }else{ message("Warning! No result has been returned by the query.") }
      }else{stop("Error: 'analyses_exp' [data.frame] should contain at least one column with 'analysis'.") }
    }
 }
  else{
    
    message("Warning! SNPs without rsIDs may not be returned.")
    
    ## Extract the data based on the input 
    ## from the custom file
    if(is.null(analyses_exp)) {
      
      if(!is.null(snp_list) & is.data.frame(snp_list)){
        
        #Create a string of analyses
        s<-paste(unique(analyses_out$analysis), collapse="','")
        s<-paste0("'",s,"'")
          
        ## Extract the outcome instruments using SNPs,
        ## obtained from the custom file/dataframe
        ## Query needs an optimization (in development)  
        if(ncol(snp_list) > 1 & search=="pos:chrom") {
          
            #Ensure the right data.type
            snp_list[1]<-as.integer(unlist(snp_list[1]))
            snp_list[2]<-as.character(unlist(snp_list[2]))
          
            #Choose only the unique instruments
            snp_list<-unique(snp_list[,c(1,2)])
          
            #Search by chromosome positions
            res <- apply(snp_list,1,function(X) odbc::dbGetQuery(dbc, paste0("SELECT g.pos,g.chrom,s.rsid,g.ref,g.alt,g.beta,g.se,
                                                                                g.pval,g.freq,g.analysis,g.entity,a.description,a.phenotype,
                                                                                a.unit,a.ncase,a.ncontrol,a.ncohort FROM gwas_results as g
                                                                              INNER JOIN sites AS s USING(pos,chrom,ref,alt)
                                                                              INNER JOIN analyses as a ON g.analysis=a.analysis
                                                                              WHERE g.pos=", X[1]," AND g.chrom='", X[2],"' AND g.analysis IN (",s,")
                                                                              ORDER BY g.analysis,g.entity;")))
          }else if(ncol(snp_list)<=1 & search=="pos:chrom") {
                stop("Error: snp_list should have at least two columns in this order: position (1) and chromosome (2).")
          }else if(ncol(snp_list)>0 & search!="pos:chrom") {
            
            #Ensure the right data.type and correct format of rsIDs
            n<-which(tolower(colnames(snp_list)) %like any% c("%rs%","%snp%"))
            snp_list<-as.data.frame(unique(snp_list[[n]]))
            colnames(snp_list)<-"SNP"
            
            #Format SNPs for DB
            snp_list$SNP<-gsub("rs","",unlist(snp_list$SNP))
            snp_list$SNP<-as.integer(unlist(snp_list$SNP))

            #Search by rsIDs
            res <- apply(snp_list,1,function(X) odbc::dbGetQuery(dbc, paste0("SELECT g.pos,g.chrom,ids.rsid,g.ref,g.alt,g.beta,g.se,
                                                                               g.pval,g.freq,g.analysis,g.entity,a.description,a.phenotype,
                                                                               a.unit,a.ncase,a.ncontrol,a.ncohort FROM gwas_results AS g
                                                                              INNER JOIN analyses AS a USING(analysis)
                                                                              INNER JOIN (SELECT * FROM sites AS s WHERE s.rsid IS NOT NULL) AS ids 
                                                                                  ON g.pos=ids.pos AND g.chrom=ids.chrom AND g.ref=ids.ref AND g.alt=ids.alt
                                                                              WHERE ids.rsid=",X[1]," AND g.analysis IN (",s,")
                                                                              ORDER BY g.analysis,g.entity;")))
          }else{
            stop("Error: snp_list should have either two columns in this order: position (1) and chromosome (2) or just one column with rsIDs (1).")
          }
          
          #Combine all results
          res_out<-NULL
          for(i in 1:length(res))
            res_out<-rbind(res_out,res[[i]])
          
          #Throw a warning if the result is empty
          #Choose only the right <analysis,entity> pairs
          if( (!is.null(res_out)) & (dim(res_out)[1] > 0) ){
            
            ids<-unique(analyses_out[,which(names(analyses_out) %in% c("analysis","entity_type","entity"))])
            
            if(is.data.frame(ids)){
              colnames(ids)[2]<-"entity"
              res_out<-dplyr::left_join(res_out,ids)
              
            }else{ message("Warning! Entity cannot be identified; 'analyses_out' dataframe is probably missing 'entity' column.") }
          }else{ message("Warning! No result has been returned by the query.") }
      }else{ stop("Error: snp_list [data.frame] cannot be empty and should have at least two columns in this order: position (1) and chromosome (2).") }
    }else{
      
      if(!is.null(snp_list)) message("Warning! These instruments will not be considered; analyses_exp should be NULL.")
      
      res_out<-NULL

      if(is.data.frame(analyses_exp)){  
        
        #Create a string of the exposure analyses
        s_exp<-paste(unique(analyses_exp$analysis), collapse="','")
        s_exp<-paste0("'",s_exp,"'")
        
        #Create a string of the outcome analyses
        s_out<-paste(unique(analyses_out$analysis), collapse="','")
        s_out<-paste0("'",s_out,"'")
        
        ## Extract the outcome instruments using SNPs
        res_out <- odbc::dbGetQuery(dbc, paste0("SELECT g.pos,g.chrom,s.rsid,g.ref,g.alt,g.beta,g.se,
                                                   g.pval,g.freq,g.analysis,g.entity,a.description,a.phenotype,
                                                   a.unit,a.ncase,a.ncontrol,a.ncohort FROM gwas_results as g
                                                 INNER JOIN (SELECT DISTINCT gg.pos_index,gg.chrom FROM gwas_results_top_hits_with_entity AS gg
                                                      WHERE gg.pval_index < ", p," AND gg.analysis IN (",s_exp,")) AS out 
                                                 ON g.pos=out.pos_index AND g.chrom=out.chrom
                                                 INNER JOIN sites AS s ON g.pos=s.pos AND g.chrom=s.chrom 
                                                            AND g.ref=s.ref AND g.alt=s.alt
                                                 INNER JOIN analyses AS a ON g.analysis=a.analysis
                                                 WHERE g.analysis IN (",s_out,")
                                                 ORDER BY g.analysis,g.entity;"))
        
        #Choose only the right <analysis,entity> pairs
        if( (!is.null(res_out)) & (dim(res_out)[1] > 0) ){
          
          ids<-unique(analyses_out[,which(names(analyses_out) %in% c("analysis","entity_type","entity"))])
          if(is.data.frame(ids)){
            colnames(ids)[2]<-"entity"
            res_out<-dplyr::left_join(res_out,ids)
            
          }else{ message("Warning! Entity cannot be identified; 'analyses_out' dataframe is probably missing 'entity' column.") }
        }else{ message("Warning! No result has been returned by the query.") }
      }else{stop("Error: 'analyses_exp' [data.frame] should contain at least one column with 'analysis'.") }
    }
  }
  
  return(res_out)
}

#' Formats the data for MR analysis.
#' 
#' Formats the database-extracted data for the Mendelian Randomization analysis (MR-Base)
#' 
#' @author Valeriia Haberland \email{valeriia.haberland@@bristol.ac.uk}
#' @param data (required): the dataframe with the exposure or outcome instruments (output of: extract_exposure or extract_outcome)
#' @param type (optional): specifies the data type, i.e. 'exposure' (default) or outcome.
#' @param rsid (optional): logical; if FALSE (default) uses pos:chrom_ref_alt as SNP identifier; if TRUE uses rsIDs as SNP identifier. 
#' @return the formatted dataframe suitable for the MR analysis
#' 
#' @export
#' @import TwoSampleMR
#' 
format_for_mr<-function(data,type="exposure",rsid=FALSE){
  
if(rsid==FALSE){
  
  if(type=="exposure"){
  
    data<-data %>% dplyr::mutate(SNP=paste0(pos_index,":",chrom,"_",ref_index,"_",alt_index),id=paste0(analysis,":",entity))  
    data<-data[toupper(data$ref_index) %in% c('A','G','C','T'),]
    data<-data[toupper(data$alt_index) %in% c('A','G','C','T'),]

    #Format the exposure for the MR analysis
    data<-TwoSampleMR::format_data(data,type="exposure",phenotype_col = "id", 
                          snp_col="SNP",beta_col="beta_index",se_col="se_index",
                          pval_col="pval_index",eaf_col="freq_index",effect_allele_col="alt_index",
                          other_allele_col="ref_index",samplesize_col ="ncohort")
  
  }
   else{
  
    data<-data %>% dplyr::mutate(SNP=paste0(pos,":",chrom,"_",ref,"_",alt),id=paste0(analysis,":",entity))  
    data<-data[toupper(data$ref) %in% c('A','G','C','T'),]
    data<-data[toupper(data$alt) %in% c('A','G','C','T'),]
    
    #Format the outcome for the MR analysis 
    data<-TwoSampleMR::format_data(data,type="outcome",phenotype_col = "id", 
                         snp_col="SNP",beta_col="beta",se_col="se",
                         pval_col="pval",eaf_col="freq",effect_allele_col="alt",
                         other_allele_col="ref",samplesize_col ="ncohort")
   }
  }
  else{
    
    data<-data %>% dplyr::mutate(id=paste0(analysis,":",entity))
    data<-data[!is.na(data$rsid),]
    data<-data[!is.null(data$rsid),]
    data$rsid<-paste0("rs",data$rsid)
    
    if(type=="exposure"){
      
      data<-data[toupper(data$ref_index) %in% c('A','G','C','T'),]
      data<-data[toupper(data$alt_index) %in% c('A','G','C','T'),]
      
      #Format the exposure for the MR analysis
      data<-TwoSampleMR::format_data(data,type="exposure",phenotype_col = "id", 
                           snp_col="rsid",beta_col="beta_index",se_col="se_index",
                           pval_col="pval_index",eaf_col="freq_index",effect_allele_col="alt_index",
                           other_allele_col="ref_index",samplesize_col ="ncohort")
      
    }
    else{
        
      data<-data[toupper(data$ref) %in% c('A','G','C','T'),]
      data<-data[toupper(data$alt) %in% c('A','G','C','T'),]
      
      #Format the outcome for the MR analysis 
      data<-TwoSampleMR::format_data(data,type="outcome",phenotype_col = "id", 
                           snp_col="rsid",beta_col="beta",se_col="se",
                           pval_col="pval",eaf_col="freq",effect_allele_col="alt",
                           other_allele_col="ref",samplesize_col ="ncohort")
    } 
  }

 return(data)
}

#' Runs the Mendelian Randomization (MR) and sensitivity analyses
#' 
#' This function runs MR and sensitivity analyses (TwoSampleMR) using the pre-formatted datasets (either format_data (TwoSampleMR) or format_for_mr (gtx)) where SNPs should be uniquely identified by either chromosome positions & alleles or rsIDs.
#' 
#' @author Valeriia Haberland \email{valeriia.haberland@@bristol.ac.uk}
#' @param exp (required): the formatted exposure instruments (output of: format_for_mr()) 
#' @param out (required): the formatted outcome instruments (output of: format_for_mr())
#' @param rsid (optional): logical; if FALSE (default) uses pos:chrom_ref_alt as a SNP identifier; 
#' if TRUE uses rsIDs as a SNP identifier (the automated clumping is enabled for this option).
#' @param clump (optional): logical; if TRUE (default) enables an automated pruning of the exposure instruments 
#' (only applicable if rsid=TRUE); if FALSE disables an automated pruning of the exposure instruments.  
#' @param hetero (optional): logical; enables/diables the heterogeneity testing (TRUE by default).
#' @param pleio (optional): logical; enables/disables the pleiotropy testing (TRUE by default).
#' @param single (optional): logical; enables/disables the single SNP MR analysis (TRUE by default).
#' @return the list of dataframes: 
#' (1) MR results
#' (2) Heterogeneity analysis
#' (3) Pleiotropy analysis
#' (4) Single SNP MR analysis
#' (5) Harmonisation data
#' 
#' @export
#' @import TwoSampleMR
#' 
#' @examples To access the dataframe, e.g. \dontrun{df[[1]]} - MR results
#' 
run_mr_gsk<-function(exp,out,rsid=FALSE,clump=TRUE,hetero=TRUE,pleio=TRUE,single=TRUE){
  
  #Pruning of the exposure instruments
  if(rsid==TRUE){
   
   if(clump==TRUE) exp<-TwoSampleMR::clump_data(exp)
   else message("Warning! An automated LD-pruning of the exposure instruments is disabled.")
    
    }else{
    message("Warning! An automated LD-pruning of the exposure instruments is disabled for the 'pos:chrom_ref_alt' identifiers.")
 }
  
  #Harmonise the data
  dat<-TwoSampleMR::harmonise_data(exp,out,action=2)
  
  #Run the MR analysis
  mr_res<-TwoSampleMR::mr(dat)
  
  if(hetero==TRUE)
    #Run the heterogeneity analysis
    mr_hetero<-TwoSampleMR::mr_heterogeneity(dat)
  else
    mr_hetero<-NULL
  
  if(pleio==TRUE)
    #Run the pleiotropy analysis
     mr_pleio<-TwoSampleMR::mr_pleiotropy_test(dat)
  else
     mr_pleio<-NULL
  
  if(single==TRUE)
    #Run the single SNP analysis
    mr_single<-TwoSampleMR::mr_singlesnp(dat)
  else
    mr_single<-NULL
  
  #Create a list of all results
  my.list<-list(mr_res,mr_hetero,mr_pleio,mr_single,dat)
  
 return(my.list) 
}

#' Visualisation of MR results
#' 
#' This methods creates three types of plots for a chosen exposure and outcome.
#' 
#' @author Valeriia Haberland \email{valeriia.haberland@@bristol.ac.uk}
#' @param mrres (required): the list of dataframes from run_for_mr method
#' @param scatter (optional): if TRUE creates a scatter plot for an exposure/outcome pair
#' @param forest (optional): if TRUE creates a forest plot for an exposure/outcome pair
#' @param funnel (optional): if TRUE creates a funnel plots for an exposure/outcome pair 
#' @param exp (optional): if "internal", the plot labels will be formatted assuming the exposure is 
#' from the internal source; it should be changed to "external" wherever is applicable.
#' @param out (optional): if "internal", the plot labels will be formatted accordingly for the outcome
#' as for the exposure; it should be changed to "external" wherever is applicable.
#' @return the list of three types of plots for one exposure vs one outcome:
#' 1) scatter plot
#' 2) forest plot
#' 3) funnel plot 
#' 
#' @import TwoSampleMR ggplot2 stringr
#' 
#' @export
#' 
#' @examples To access the plots, e.g. \dontrun{plots[[1]]} - a scatter plot
#' 
visual_for_mr <- function(mrres,scatter=TRUE,forest=TRUE,funnel=TRUE,exp="internal",out="internal",dbc=getOption("gtx.dbConnection", NULL)){
  
  ##Format plot labels for the internal data;
  ##if "external" is stated for "exp" (exposure) or "out" (outcome), 
  ##then the labels will not be formatted for exp or out respectively
  labels<-format_plot_labels(mrres, exp, out, dbc = dbc)
  exp_label<-labels[1]
  out_label<-labels[2]
  
  ##Build a scatter plot
  if(scatter==TRUE) {
    p1 <- TwoSampleMR::mr_scatter_plot(mrres[[1]], mrres[[5]])
    p1 <- p1[[1]] + labs(x = stringr::str_wrap(paste0("SNP effect on ",exp_label), width = 40), y=stringr::str_wrap(paste0("SNP effect on ",out_label), width = 35))
  } else {
    p1 <- NULL
    message("Warning! Scatter plots will not be created.")
  }
  
  ##Build a forest plot
  if(forest==TRUE & !is.null(mrres[[4]])){
    p2 <- TwoSampleMR::mr_forest_plot(mrres[[4]])
    p2 <- p2[[1]] + labs(x = stringr::str_wrap(paste0("MR effect size for ",exp_label," on ",out_label), width = 70))
  } else {
    p2 <- NULL
    message("Warning! Forest plots will not or cannot be created.")
  }
  
  ##Build a funnel plot
  if(funnel==TRUE & !is.null(mrres[[4]])){
    p3 <- TwoSampleMR::mr_funnel_plot(mrres[[4]])
  } else {
    p3 <- NULL
    message("Warning! Funnel plots will not or cannot be created.")
  }
  
  ##Return the list of plots
  list.plots <- list(p1,p2,p3)
  
  return(list.plots)
}

#'Format plot labels for internal data
#'
#'@author Valeriia Haberland \email{valeriia.haberland@@bristol.ac.uk}
#'@param mrres (required): MR results as derived from run_mr_gsk()
#'@param exp_flag (optional): a flag that indicates a source of the exposure, i.e. internal or external;
#'labels are fomatted for the internal data only
#'@param out_flag (optional): a flag that indicates a source of the outcome, i.e. internal or external;
#'labels are formatted for the internal data only
#'@return labels for the chosen exposure and outcome; if external, they will not be modified.
#'
#'@keywords internal
#'
format_plot_labels<-function(mrres,exp_flag="internal",out_flag="internal",dbc=getOption("gtx.dbConnection", NULL)){
  
  ##Reformat labels for internal usage
  if(exp_flag=="internal"){
    exp_id<-gsub(':.+','',unique(mrres[[1]]$exposure))
    exp_label<-gtxanalysis_label(analysis = exp_id, nlabel = FALSE, dbc = dbc)
  } else {
    exp_label<-as.character(unique(mrres[[1]]$exposure))
  }
  
  if(out_flag=="internal"){
    out_id<-gsub(':.+','',unique(mrres[[1]]$outcome))
    out_label<-gtxanalysis_label(analysis = out_id, nlabel = FALSE, dbc = dbc)
  } else {
    out_label<-as.character(unique(mrres[[1]]$outcome))
  }
  
  labels<-c(exp_label,out_label)
  
  return(labels)
}  

#' Search study by keyword.
#' 
#' Finds a unique identifier for the study using a keyword or keyphrase. NOTE! keyword or keyphrase search is quite naive at this moment and may not return desirable results.
#' 
#' @author Valeriia Haberland \email{valeriia.haberland@@bristol.ac.uk}
#' @param str (required): a keyword/keyphrase to search for the study of interest
#' @param dbc (default): other parameters 
#' @return the study meta-information which then can be used in extract_outcome()
#' or extract_exposure() methods
#' 
#' @keywords internal
#' 
extract_study_info <- function(str,dbc=getOption("gtx.dbConnection", NULL)) {
  
  #Extract all records with the given description
  studies<-odbc::dbGetQuery(dbc, paste0("SELECT * FROM analyses WHERE description LIKE '%", str, "%';"))
  
  if( nrow(studies) < 1 ) {
    message("No results were returned! Please, try a shorter phrase.")  
  }
  
  return(studies)
}
