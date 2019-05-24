#' PICCOLO: Colocalization of GWAS and eQTL/pQTL signals without summary stats
#'
#'   PICCOLO estimates the probability that two genetic signals between GWAS and e/pQTLs 
#' . are shared using colocalization approach. A proper colocalization analysis requires 
#' . full summary stats for both the GWAS and QTL studies. Unfortunately, studies frequently 
#'   do not share full summary stats which makes a traditional colocalization analysis impossible..
#' . In this case, PICCOLO enable to do a colocalization using PICS.
#'      
#' @param rs SNP rsid. Single string or vector.
#' @param chrom SNP chromosome. Single string or vector. 
#' @param pos SNP position. Single string or vector.
#' @param pval the association p-value for the SNP. Single string or vector.
#' @param ancestry Ancestry information: EUR/ASN/AFR. Currently, only EUR is available. OPTIONAL!!
#' @param indication Associated phenotype or trait info. OPTIONAL!!
#' @param dbc Database connection. Default: getOption("gtx.dbConnection", NULL) 
#'   
#' @return  
#'  a data frame containing the result of the
#'  PICS calculation.
#' 
#' @author Kijoung Song \email{kys21207@@gsk.com}
#'
#' @export

piccolo <- function(chrom,pos,rs,pval,ancestry,indication,dbc=getOption("gtx.dbConnection", NULL)){
  
  # check database connection
  gtxdbcheck(dbc)
  
  # Check input columns
  gtx_debug("piccolo | validating input.")
  if((missing(rs) & missing(chrom) & missing(pos)) | (missing(rs) & missing(chrom)) |(missing(rs) & missing(pos))) { stop("piccolo | must specify either: rs or chr/pos", call. = FALSE) }
  if(!missing(rs) & !missing(chrom) & !missing(pos))  { stop("piccolo | Duplicate input: input either rs or chr/pos", call. = FALSE) }
  if(missing(pval))     { stop("piccolo | input is missing: pval", call. = FALSE) }
  if(!missing(rs)){
    if(missing(ancestry)) { ancestry <- c(rep("EUR", length(rs))) }
    if(missing(indication)) { indication <- c(rep(NA_character_, length(rs)))}
    tmp.list <- list(rs,pval,ancestry,indication)
    if(all(sapply(tmp.list,length)==length(tmp.list[[1]]))){
      snpID <- rs
      input <- tibble(snpID, pval, ancestry, indication)
    }else{
      stop("piccolo1 | the length of each input is different.", call. = FALSE)
    }
    # exclude any duplicated input
    input <- input %>% distinct(snpID, indication, .keep_all = TRUE)
    
  }else if(!missing(chrom) & !missing(pos)){
    if(missing(ancestry)) { ancestry <- c(rep("EUR", length(pos))) }
    if(missing(indication)) { indication <- c(rep(NA_character_, length(pos)))}
    tmp.list <- list(chrom,pos,pval,ancestry,indication)
    if(all(sapply(tmp.list,length)==length(tmp.list[[1]]))){
      snpID <- paste(chrom,pos,sep=":")
      input <- tibble(snpID, pval, ancestry, indication)
    }else{
      stop("piccolo2 | the length of each input is different.", call. = FALSE)
    }
    # exclude any duplicated input
    input <- input %>% distinct(snpID, indication, .keep_all = TRUE)
    
  }
  #print(input)
  input$ancestry <- ifelse(!is.na(input$ancestry) , as.character(input$ancestry),"EUR")
  input <- input[!is.na(input$snpID) & !is.na(input$pval) ,]
  
  gwas.pics <- tryCatch(pics_calc(input),error=function(e) NULL)
  if(is.null(gwas.pics)) { stop("All SNP IDs are not available in the current LD reference", call. = FALSE)}
  
  idx.snps <- gwas.pics %>% group_by(chrom_idx,pos_idx,pval_idx) %>% summarise(start = min(pos),end = max(pos))
  idx.nearst.genes <- 
    idx.snps %>% 
    group_by(chrom_idx,pos_idx,pval_idx) %>% 
    do( sqlWrapper(dbc,paste0("SELECT ensemblid FROM genes WHERE genetype = 'protein_coding' AND chrom = '", .$chrom_idx, "' AND (",
                              "(pos_start >= ",     .$start - 5e5, " AND pos_start <= ", .$end + 5e5, ")",
                              " OR (pos_end >= ",   .$start - 5e5, " AND pos_end <= ",   .$end + 5e5, ")", 
                              " OR (pos_start >= ", .$start - 5e5, " AND pos_end <= ",   .$end + 5e5, ")",
                              ")"),uniq = F, zrok = TRUE))
  
  idx.nearst.genes$ID <- paste(idx.nearst.genes$chrom_idx,idx.nearst.genes$pos_idx,idx.nearst.genes$pval_idx,sep="_")
  #gwas.pics$ID <- paste(gwas.pics$chrom_idx,gwas.pics$pos_idx,gwas.pics$pval_idx,sep="")
  
  
  if(length(unique(idx.nearst.genes$ensemblid)) > 1000){
    ens2List <- split(unique(idx.nearst.genes$ensemblid),cut(1:length(unique(idx.nearst.genes$ensemblid)),ceiling(length(unique(idx.nearst.genes$ensemblid))/1000),F))
  }else{
    ens2List <- split(unique(idx.nearst.genes$ensemblid),1)
  }
  tmp <- list()
  for(j in 1:length(ens2List)){
    # TODO fix database
    tmp[[j]]<-sqlWrapper(dbc,paste0("SELECT * FROM pics_qtls WHERE entity IN ('",paste(ens2List[[j]],collapse="', '"),"') "),uniq = F, zrok = FALSE)
  }
  pics.qtls <- int_fastDoCall("rbind",tmp)
  pics.qtls$ID <- paste(pics.qtls$tissue,pics.qtls$entity,pics.qtls$rsid_idx,sep="_")
  
  tmp00 <- list()
  n1 <- 1
  for(j in unique(idx.nearst.genes$ID)){
    dta1 <- subset(gwas.pics, paste(chrom_idx,pos_idx,pval_idx,sep="_") == j, c("pos","pics","snpID_idx","rsid_idx","pval_idx","chrom_idx","ancestry","indication"))
    dta1 <- dta1[order(-dta1$pics),]
    names(dta1) <- c("pos1","pics1","snpID","rsid","pval","chrom","ancestry","indication")
    
    tmp <- subset(idx.nearst.genes, ID == j)
    x <- subset(pics.qtls, entity %in% tmp$ensemblid)
    tmp01 <- list()
    n2 <- 1
    for(k in unique(x$ID)){
      dta2 <- subset(x, ID %in% k)
      dta2 <- dta2[order(-dta2$pics),]
      dta2 <- dta2[,c("pos","pics","rsid_idx","pval_idx","tissue","hgncid","entity","pmid_idx")]
      names(dta2) <- c("pos2","pics2","eqtl_rsid","eqtl_pval","tissue","hgnc_symbol","ensembl_ID","pubmed_ID")
      tmp <- int_coloc_pics_lite(dta1,dta2,pics1="pics1",pics2="pics2",rsid1="pos1",rsid2="pos2")
      tmp01[[n2]] <- bind_cols(dta1[1,],dta2[1,],tmp)
      n2 <- n2+1
    }
    tmp00[[n1]] <- int_fastDoCall("rbind",tmp01)
    n1 <- n1+1
  }
  res <- do.call("rbind",tmp00)
  res <- res[,c("snpID","rsid","pval","chrom","pos1","ancestry","indication","eqtl_rsid","eqtl_pval","pos2","hgnc_symbol","ensembl_ID","tissue","pubmed_ID","H3","H4")]
  names(res) <- c("gwas_input","gwas_rsid","gwas_pval","gwas_chrom","gwas_pos","ancestry","indication","qtl_rsid","qtl_pval","qtl_pos","hgncid","ensemblid","tissue","pmid","H3","H4")
  check.missing.gwas.snp <- subset(input, !(input$snpID %in% unique(res$gwas_input)),c("snpID","pval","ancestry","indication") )
  names(check.missing.gwas.snp) <- c("gwas_input","gwas_pval","ancestry","indication")
  if(nrow(check.missing.gwas.snp) >= 1) res <- int_sbind(check.missing.gwas.snp,res)
  
  return(res)
}

#' PICS calculation using the linkage information of 10,000 UKB samples (EUR only)
#' 
#'   The PICS algorithm calculates the most likely causal SNPs given the observed association signal
#' . at a locus. For an associated locus, enter the most highly-associated SNP 
#'   (referred to as the index SNP) and the strength of association. Using the linkage information
#'   in 10,000 UKB samples, the algorithm identifies the SNPs that are most likely to be the causal variants
#'   responsible for the association (PICS_Probability).
#'
#'   See \strong{Genetic and Epigenetic Fine-Mapping of Causal Variants in Autoimmune Disease} 
#'   by Kyle Kai-How Farh,et al. Nature 518, 337–343 (19 February 2015) at
#'   \url{http://www.nature.com/nature/journal/v518/n7539/full/nature13835.html#close}
#'
#' @param data Data frame with columns rs|chr:pos, pval, ancestry and indication. Use the same column names for your data. However, ancestry and indication are optional.
#' @param dbc Database connection. Default: getOption("gtx.dbConnection", NULL) 
#'   
#' @return  
#'  a data frame containing the result of the
#'  PICS calculation.
#' 
#' @author Kijoung Song \email{kys21207@@gsk.com}
#'
#' @export

pics_calc <- function(index.data,dbc=getOption("gtx.dbConnection", NULL)){
  gtx_debug("pics_calc | starting pics calc")
  rs.snpid <- subset(index.data, grepl("rs", index.data$snpID))
  rs.snpid.tmp <- gsub("rs", "", rs.snpid$snpID)
  .snpid <- subset(index.data, !grepl("rs", index.data$snpID))
  dta.ext <- NULL
  if (nrow(rs.snpid) > 0) {
    if (length(unique(rs.snpid.tmp)) > 1000) {
      rsList <- split(unique(rs.snpid.tmp), cut(1:length(unique(rs.snpid.tmp)), 
                                                ceiling(length(unique(rs.snpid.tmp))/1000), F))
    } else {
      rsList <- split(unique(rs.snpid.tmp), 1)
    }
    tmp01 <- list()
    for(i in 1:length(rsList)){
      x <- tryCatch(sqlWrapper(dbc, paste0("SELECT chrom,pos,rsid FROM sites WHERE rsid IN (", 
                                           paste(rsList[[i]], collapse = ","), ")"), uniq = F, 
                               zrok = FALSE), error = function(e) NULL)
      if (!is.null(x)) {
        x = x[!duplicated(x$rsid), ]
        x$rsid1 <- paste("rs", dta.ext$rsid, sep = "")
        x$snpID <- dta.ext$rsid1
        x$rsid <- NULL
        tmp01[[i]] <- x
      } else {
        gtx_warn("pics_calc: all rsid's are missing, check your rsid's carefully!")
      }
    }
    dta.ext <- int_fastDoCall("rbind", tmp01)
  }
  if (nrow(.snpid) > 0) {
    .snpid$chrom <- unlist(lapply(strsplit(as.character(.snpid$snpID), 
                                           ":"), function(x) x[[1]]))
    .snpid$pos <- as.numeric(unlist(lapply(strsplit(as.character(.snpid$snpID), 
                                                    ":"), function(x) x[[2]])))
    tmp00 <- list()	
    n <- 1	
    for (i in unique(.snpid$chrom)) {
      .snpid.sub <- subset(.snpid, chrom == i)
      if (length(unique(.snpid.sub$pos)) > 1000) {
        rsList <- split(unique(.snpid.sub$pos), cut(1:length(unique(.snpid.sub$pos)), 
                                                    ceiling(length(unique(.snpid.sub$pos))/1000), F))
      } else {
        rsList <- split(unique(.snpid.sub$pos), 1)
      }
      tmp01 <- list()
      for(j in 1:length(rsList)){
        x = tryCatch(sqlWrapper(dbc, paste0("SELECT chrom,pos,rsid FROM sites WHERE chrom = '", 
                                            i, "' AND pos IN (", paste(rsList[[j]], collapse = ","), 
                                            ")"), uniq = F, zrok = FALSE), error = function(e) NULL)
        if (!is.null(x)) {
          x <- x %>% distinct(chrom, pos, .keep_all = TRUE) %>% 
            mutate(snpID = paste(chrom, pos, sep = ":"), 
                   rsid1 = paste("rs", rsid, sep = "")) %>% 
            select(-rsid)
          tmp01[[j]] <- x
        }
      }
      tmp00[[n]] <- int_fastDoCall("rbind", tmp01)
      n <- n + 1
    }
    res <- int_fastDoCall("rbind", tmp00)
    if (is.null(res)) 
      gtx_warn("pics_calc: all chr:pos' are missing, check your chr:pos' carefully!")
    if (is.null(dta.ext)) {
      dta.ext <- res
    } else {
      dta.ext <- rbind(dta.ext, res)
    }
  }
  if (!is.null(dta.ext)) {
    index.data = merge(index.data, dta.ext, by = "snpID")
  } else {
    stop("pics_calc: all rs's are missing, carefully check the input data, specially rs.", 
         call. = FALSE)
  }
  # exclude chrom = "X" and "XY" 
  index.data <- subset(index.data, chrom %in% as.character(seq(1,22))) 
  
  tmp.ld <- index.data
  tmp.ld$chrom2 <- tmp.ld$chrom
  tmp.ld$pos2 <- tmp.ld$pos
  tmp.ld$r <- 1
  tmp.ld$r2 <- 1
  tmp.ld <- tmp.ld[, c("snpID", "rsid1", "pval", "chrom", "pos", 
                       "chrom2", "pos2", "ancestry", "indication", "r", "r2")]
  names(tmp.ld) <- c("snpID", "rsid1", "pval", "chrom1", "pos1", 
                     "chrom2", "pos2", "ancestry", "indication", "r", "r2")
  tmp00 <- list()
  n <- 1
  for (i in unique(index.data$chrom)) {
    sub.dta <- subset(index.data, chrom %in% i)
    if (length(unique(sub.dta$pos)) > 1000) {
      pos1List <- split(unique(sub.dta$pos), cut(1:length(unique(sub.dta$pos)), 
                                                 ceiling(length(unique(sub.dta$pos))/1000), F))
    } else {
      pos1List <- split(unique(sub.dta$pos), 1)
    }
    tmp01 <- list()
    for (j in 1:length(pos1List)) {
      sub.ld = sqlWrapper(dbc, paste0("SELECT * FROM ld WHERE pos1 IN (", 
                                      paste(pos1List[[j]], collapse = ","), ") and chrom1='", 
                                      i, "'"), uniq = F, zrok = FALSE)
      sub.ld[!duplicated(sub.ld[, c("pos2", "pos1")]),]
      tmp01[[j]] <- sub.ld
    }
    tmp00[[n]] <- int_fastDoCall("rbind", tmp01)
    n <- n+1
  }
  all.ld <- int_fastDoCall("rbind", tmp00)
  all.ld$r2 <- all.ld$r^2
  all.ld <- all.ld[order(all.ld$pos1, -all.ld$r2), ]
  all.ld <- merge(all.ld, index.data, by.x = c("chrom1", "pos1"), 
                  by.y = c("chrom", "pos"), all.x = T)
  all.ld <- int_sbind(tmp.ld, all.ld)
  all.ld <- subset(all.ld, r2 > 0.5)
  all.ld$pval = ifelse(all.ld$pval == 0, 1e-250, all.ld$pval)
  all.ld$SD = sqrt(1 - abs(all.ld$r)^6.4) * sqrt(-log10(all.ld$pval))/2
  all.ld$Mean = all.ld$r2 * (-log10(all.ld$pval))
  all.ld$prob = ifelse(all.ld$SD == 0, 0.8, 1 - pnorm(-log10(all.ld$pval), 
                                                      all.ld$Mean, all.ld$SD))
  prob.sum <- all.ld %>% group_by(chrom1, pos1, pval) %>% summarise(prob_sum = sum(prob))
  all.ld = merge(all.ld, prob.sum, by = c("chrom1", "pos1", 
                                          "pval"))
  all.ld$pics = all.ld$prob/all.ld$prob_sum
  tmp00 <- list()
  n <- 1
  for (i in unique(all.ld$chrom2)) {
    sub.dta <- subset(all.ld, chrom2 %in% i)
    if (length(unique(sub.dta$pos2)) > 1000) {
      pos2List <- split(unique(sub.dta$pos2), cut(1:length(unique(sub.dta$pos2)), 
                                                  ceiling(length(unique(sub.dta$pos2))/1000), F))
    } else {
      pos2List <- split(unique(sub.dta$pos2), 1)
    }
    tmp01 <- list()        
    for (j in 1:length(pos2List)) {
      tmp <- sqlWrapper(dbc, paste0("SELECT chrom,pos,rsid FROM sites WHERE pos IN (", 
                                    paste(pos2List[[j]], collapse = ","), ") and chrom='", 
                                    i, "'"), uniq = F, zrok = FALSE)
      tmp1 <- subset(tmp, !is.na(rsid))
      tmp1 <- subset(tmp1, !duplicated(paste(chrom, pos, 
                                             sep = "_")) & !duplicated(rsid))
      tmp2 <- subset(tmp, is.na(rsid))
      tmp3 <- subset(tmp2, !(paste(chrom, pos, sep = "_") %in% 
                               paste(tmp1$chrom, tmp1$pos, sep = "_")))
      tmp01[[j]] <- rbind(tmp1, tmp3)
    }
    tmp00[[n]] <- int_fastDoCall("rbind", tmp01)
    n <- n+1
  }
  snp.anno <- int_fastDoCall("rbind", tmp00)
  snp.anno$snp2 <- ifelse(!is.na(snp.anno$rsid), paste("rs", 
                                                       snp.anno$rsid, sep = ""), snp.anno$rsid)
  names(snp.anno) <- c("chrom2", "pos2", "rsid", "snp2")
  pics.result <- merge(all.ld, snp.anno, by = c("chrom2", "pos2"))
  pics.result <- pics.result[order(pics.result$snpID, -pics.result$r2), 
                             c("chrom2", "pos2", "snp2", "pics", "ancestry", "indication", 
                               "snpID", "rsid1", "pval", "chrom1", "pos1")]
  names(pics.result) <- c("chrom", "pos", "rsid", "pics", "ancestry", 
                          "indication", "snpID_idx", "rsid_idx", "pval_idx", "chrom_idx", 
                          "pos_idx")
  return(pics.result)}

#' int_coloc_pics_lite 
#'   Test for colocalization of two PICS sets
#' 
#' @return only H3 & H4 posteriors
#' @param data1,data2  PICS sets from read.pics or download.pics
#' @param pics1,pics2  column name to pull PICS prob from. Default = "PICS_probability"
#' @param rsid1,rsid2  column name to pull rsid from.      Default = "Linked_SNP"
#' @param rounded   Decimal points to round posteriors to
#' @param priorc1   Prior probability for colocalization with siganl for data1  Default = 1e-4
#' @param priorc2   Prior probability for colocalization with siganl for data2 Default = 1e-4
#' @param priorc12  Prior probability for colocalization of both signals.   Default = 1e-5
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}  


int_coloc_pics_lite <- function(data1,
                                data2,
                                pics1    = "PICS_probability", # column header for poster probabilities in data1
                                pics2    = "PICS_probability", # column header for poster probabilities in data2
                                rsid1    = "Linked_SNP",       # column header for snps in LD in data1
                                rsid2    = "Linked_SNP",       # column header for snps in LD in data2
                                rounded  = 6,
                                priorc1  = 1e-4, 
                                priorc2  = 1e-4, 
                                priorc12 = 1e-5
) {
  stopifnot(exists("data1") & exists("data2"))
  if(is.logical(data1)){
    if(is.na(data1)){
      gtx_warn("int_coloc_pics_lite: data1 is NA, skipping coloc.")
      return(list(results = NA, nvariants = NA))
    }
  }
  if(is.logical(data2)){
    if(is.na(data2)){
      gtx_warn("int_coloc_pics_lite: data2 is NA, skipping coloc.")
      return(list(results = NA, nvariants = NA))
    }
  }
  pics <- int_harmonize_pics(data1, 
                             data2, 
                             opts <- data.frame(rsid1 = rsid1, rsid2 = rsid2, pics1 = pics1, pics2 = pics2, stringsAsFactors = FALSE))
  
  nv <- dim(pics)[1]
  res <- data.frame(prior = norm1(c(priorc1*priorc2*nv*(nv - 1), priorc12*nv)),
                    bf    = c((sum(pics[[1]])*sum(pics[[2]]) - sum(pics[[1]]*pics[[2]]))/(nv*(nv - 1)), 
                              sum(pics[[1]]*pics[[2]])/nv))
  res$bf <- res$bf/res$bf[1]
  res$posterior <- norm1(res$prior*res$bf)
  if (is.finite(rounded)) {
    res$posterior = round(res$posterior, rounded)
  }
  return(data.frame(H3=res$posterior[1], H4=res$posterior[2]))
}

#' int_harmonize_pics
#' 

int_harmonize_pics <- function(data1,
                               data2, 
                               opts = data.frame(pics1 = "PICS_probability",
                                                 pics2 = "PICS_probability",
                                                 rsid1 = "Linked_SNP",
                                                 rsid2 = "Linked_SNP",
                                                 stringsAsFactors = FALSE)
){
  ids <- unique(c(data1[[opts$rsid1]], data2[[opts$rsid2]]))
  tmp <- as.data.frame(matrix(data = NA, nrow = length(ids), ncol = 2))
  pp1 <- if (opts$pics1==opts$pics2) paste(opts$pics1, ".1", sep = "") else opts$pics1
  pp2 <- if (opts$pics1==opts$pics2) paste(opts$pics2, ".2", sep = "") else opts$pics2
  colnames(tmp) <- c(pp1, pp2)
  for(n in 1:length(ids)){
    tmp[[pp1]][n] <- if(!is.na(match(ids[n], data1[[opts$rsid1]]))) data1[which(data1[[opts$rsid1]]==ids[n]),][[opts$pics1]][1] else 0
    tmp[[pp2]][n] <- if(!is.na(match(ids[n], data2[[opts$rsid2]]))) data2[which(data2[[opts$rsid2]]==ids[n]),][[opts$pics2]][1] else 0 
  }
  res <- as.data.frame(cbind(
    norm1(tmp[[pp1]], log = FALSE),
    norm1(tmp[[pp2]], log = FALSE)
  ))
  colnames(res) <- c(pp1, pp2)
  rownames(res) <- ids
  return(res)
}



#' Combine Shape Objects
#'
#'   Combine shape objects into one shape object. It works analogous to rbind.
#'
#' @author  Martijn Tennekes
#'


int_sbind = function(x, y, fill=NA) {
  int_sbind.fill = function(d, cols){ 
    for(c in cols)
      d[[c]] = fill
    d
  }
  
  x = int_sbind.fill(x, setdiff(names(y),names(x)))
  y = int_sbind.fill(y, setdiff(names(x),names(y)))
  
  rbind(x, y)
}

#' An Alternative To The Internal Do.Call
#'  
#'   The do.call can be somewhat slow, especially when working with large objects. 
#'   This function is based upon the suggestions from Hadley Wickham on the R mailing list, 
#'   see here. Also thanks to Tommy at StackOverflow for suggesting how to handle double 
#' . and triple colon operators, ::, further enhancing the function.
#'  
#' @param what   either a function or a non-empty character string naming the function to be called.
#' @param args   a list of arguments to the function call. The names attribute of args gives the argument names.
#' @param quote  a logical value indicating whether to quote the arguments.
#' @param envir  an environment within which to evaluate the call. This will be most useful if what is a character string and the arguments are symbols or quoted expressions.
#'
#' @author  Max Gorden 
#' 


int_fastDoCall <- function(what, args, quote = FALSE, envir = parent.frame()){
  if (quote)
    args <- lapply(args, enquote)
  
  if (is.null(names(args))){
    argn <- args
    args <- list()
  }else{
    # Add all the named arguments
    argn <- lapply(names(args)[names(args) != ""], as.name)
    names(argn) <- names(args)[names(args) != ""]
    # Add the unnamed arguments
    argn <- c(argn, args[names(args) == ""])
    args <- args[names(args) != ""]
  }
  
  if (class(what) == "character"){
    if(is.character(what)){
      fn <- strsplit(what, "[:]{2,3}")[[1]]
      what <- if(length(fn)==1) {
        get(fn[[1]], envir=envir, mode="function")
      } else {
        get(fn[[2]], envir=asNamespace(fn[[1]]), mode="function")
      }
    }
    call <- as.call(c(list(what), argn))
  }else if (class(what) == "function"){
    f_name <- deparse(substitute(what))
    call <- as.call(c(list(as.name(f_name)), argn))
    args[[f_name]] <- what
  }else if (class(what) == "name"){
    call <- as.call(c(list(what, argn)))
  }
  
  eval(call,
       envir = args,
       enclos = envir)
}