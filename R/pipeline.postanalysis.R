postanalysis.pipeline<- function(configFile){
  #load config
  config <- read.table(configFile, sep ="=", as.is = T, strip.white = TRUE,stringsAsFactors = FALSE, quote = "", row.names = 1)
  options(gtxpipe.project = config["project", 1])
  options(gtxpipe.user = config["user", 1])
  options(gtxpipe.email = config["email", 1])
  options(gtxpipe.date = as.character(getOption("gtxpipe.date", format(Sys.Date(), "%Y-%b-%d")))[1])
  
  if (!is.null(config["chunkMb", 1]) & !is.na(config["chunkMb", 1])) {
    chunkMb <- as.numeric(as.character(config["chunkMb", 1]))
  } else {
    chunkMb <- 4
  }
  genodir <-  config["genotypes",1]
  threshold.MAF <- as.numeric(as.character(config["threshold.MAF", 1]))
  threshold.Rsq <- as.numeric(as.character(config["threshold.Rsq", 1]))
  plot <- as.logical(config["plot", 1])
  flanking <- as.numeric(as.character(config["flanking", 1]))
  
  ## Get location of gtx library and load
  if (!is.null(config["gtxloc", 1]) & !is.na(config["gtxloc", 1])) {
    gtxloc = config["gtxloc",1]
  } else {
    gtxloc = "/GWD/appbase/projects/statgen/GXapp/G-P_assoc_pipeline/R-packages/x86_64-unknown-linux-gnu-library/3.0"
  }
  .libPaths(gtxloc)
  library(gtx)
  library(survival)
  library(ordinal)
  
  if (!is.null(config["GC", 1]) & !is.na(config["GC", 1])) {
    GC<- as.logical(config["GC",1])
  } else {
    GC <- T # use GCed results
    
  }
  
  
  
  #set working directory, should be the folder with analyses, input, outputs subfolders
  if(!is.null(config["dir", 1]) & !is.na(config["dir", 1]))  wkdir <-config["dir", 1] 
  else     {wkdir <- getwd()}
  setwd(wkdir)
  
  if(plot) {
    dir.create("plot")
    load("/GWD/appbase/projects/statgen/GXapp/G-P_assoc_pipeline/GDCgtx/package_update/lili/recombHapMapII_GRCh37.RData")#change to gtx location
    load("/GWD/appbase/projects/statgen/GXapp/G-P_assoc_pipeline/GDCgtx/package_update/lili/hg19.GENCODE19.RData")#change to gtx location
    source("/GWD/appbase/projects/statgen/GXapp/G-P_assoc_pipeline/GDCgtx/package_update/lili/regionplot.R")
  }
  
  #obtain list of subjects for LD calculation
  ld.subj <- NULL
  if(!is.null(config["ldsubjects", 1]) & !is.na(config["ldsubjects", 1]) & file.exists(config["ldsubjects", 1])){
    ld.subj<- read.table(config["ldsubjects", 1], as.is = T)[[1]]
    message("Using subjects from ", config["ldsubjects", 1], " for LD calculation.")
  }  else{ 
    message("Using all subjects for LD calculation.")
  }
  
  #load post analysis file
  models <- NULL
  varlist.all <- NULL
  if(!is.null(config["postanalysis", 1]) & !is.na(config["postanalysis", 1])&file.exists(config["postanalysis", 1])) {
    models<- read.table(config["postanalysis",1], na.strings=c('NA',''),colClasses="character", sep = "\t", 
                        as.is = T, header = T, strip.white = TRUE, quote = "", fill = TRUE)
    for (i in 1:nrow(models)) {
      if (!is.na(models[i,"varlist"])) {
        curr <- read.table( models[i,"varlist"], stringsAsFactors = FALSE, fill = TRUE, sep = "\t")[ , 1]
        models[i,"varlist"] = paste(curr, collapse = " ")
        varlist.all<- unique(c(varlist.all, curr))
      }
      #Need to convert NA values to empty strings else parse attempts in gtxpipe function will fail
      if (is.na(models[i,"contrasts"]) || is.null(models[i,"contrasts"])) {
        models[i,"contrasts"] = ""
      }
      if (is.na(models[i,"groups"]) || is.null(models[i,"groups"])) {
        models[i,"groups"] = ""
      }
    }
    ## Reorder and re-label columns to match what expected by gtxpipe
    models <- models[c("analysis",  "groups", "contrasts", "varlist")]
    names(models)[1] <- "model"
  }
  if(is.null(models)) {
    message("Missing postanalysis file!")
    exit()
  }
  
  #obtain dose data and LD info for all variants
  dir.create("data")
  chrBeginEndDir<-"/GWD/appbase/projects/statgen/RD-MDD-GX_PUBLIC/1KG/share.sph.umich.edu/1000genomes/fullProject/2013.05.02" #/GWD/appbase/projects
  #chunkMb<- 4 #7 #4
  dose.all <-Reduce(function(x, y) merge(x, y, all = T), 
                    lapply(1:length(varlist.all), function (idx){
    prefix.curr <- paste("data", gsub("[:*]",".", varlist.all[idx]), sep = "/")
    file.dose.curr <- paste(prefix.curr, "dose.txt", sep = ".")
    if(file.exists(file.dose.curr) & file.exists(paste(prefix.curr, "info.flanking.txt", sep = "."))) {
      message(Sys.time(), ": Loading dosage data from ",file.dose.curr)
      curr$dose <- read.delim(file.dose.curr, header = T, as.is = T, check.names = F)
    }else {
      message(Sys.time(), ": Obtaining dosage data for ",varlist.all[idx])
      curr <- getDoseandLD(varlist.all[idx], flanking, ld.subj, genodir, chunkMb, chrBeginEndDir)
      message(Sys.time(), ": writing info and LD data to ",paste(prefix.curr, "info.flanking.txt", sep = "."))
      write.table(curr$LD, file = paste(prefix.curr, "info.flanking.txt", sep = "."), sep = "\t", row.names = F, quote = F)
      message(Sys.time(), ": writing dosage data to ",file.dose.curr)
      write.table(curr$dose, file = file.dose.curr, sep = "\t", row.names = F, quote = F)
    }
    return (curr$dose)
  })) 
  
  #for each model
  #save phenogeno file per model per group
  #save assoc result per model for all groups and contrasts
  #plot per variant: 
  for( idx in 1:nrow(models)){
    model <- models[idx, "model"]
    varlist <- unlist(strsplit(models[idx, "varlist"], " "))
    groups <- tokenise.whitespace(models[idx, "groups"])
    acontrasts <- tokenise.whitespace(models[idx, "contrasts"])
    acontrasts.bad <- sapply(acontrasts, function(contrast1) return(length(unlist(strsplit(contrast1, "/"))) != 2))
    if (any(acontrasts.bad)) {
      stop('Model "', models[idx, "model"], '" has contrasts(s) ',
           paste(acontrasts[acontrasts.bad], collapse = ", "),' not in format group1/group2')
    }
    
    #obtain meta data for each group:modelname, groupname, modelcall, phenotype, association results
    agroups <- unique(c(groups, unlist(strsplit(acontrasts, "/"))))
    gres <- lapply(agroups, function(group){
      curr<- padata.permodelgroup(paste("analyses",model, group,ifelse(GC, "ALL.out.gc.txt.gz", "ALL.out.txt.gz"), sep ="/" ),  
                                  varlist, flanking,  threshold.MAF, threshold.Rsq, GC = GC)
      cols.curr<- c("pvalue", "beta", "SE",  "analysed.Freq1", "analysed.Rsq", "MAF")
      if(GC) cols.curr<- c(cols.curr, c("pvalue.GC",  "SE.GC"))
      for(c in cols.curr) 
        setnames(curr$assoc,c, paste(c, curr$group, sep = ".") )
      idx.all <- names(dose.all)[-1] %in% varlist
      if(any(idx.all))
        curr$phenotype<- merge(curr$phenotype, dose.all[c(T, names(dose.all)[-1] %in% varlist)], 
                               by.x = names(curr$phenotype)[1], by.y = names(dose.all)[1], all.x = T)
      else message(Sys.time(), ": No dosage data were found for ",model, " in ", group, " for ", paste(varlist, collapse = ", "), ".")
      message(Sys.time(), ": writing pheno/geno data to ",paste("data/PhenoGeno.",model, ".", group, ".txt" , sep = ""))
      write.table(curr$phenotype, file = paste("data/PhenoGeno.",model, ".", group, ".txt" , sep = ""), sep = "\t", row.names = F, quote = F)
      return (curr)
    })
    names(gres)<- agroups
    
    #merge association results from all groups
    #c("CHROM", "POS", "SNP","Al1", "Al2") in common and c("pvalue", "beta", "SE",  "analysed.Freq1", "analysed.Rsq", "MAF") for each group
    assoc <- Reduce(function(x, y) merge(x, y, all = T, by = c( "SNP","Al1", "Al2","CHROM", "POS"), allow.cartesian = T), 
                    lapply(gres, function(padata){return (padata$assoc)}))   
    #get contrast result
    if(length(acontrasts)>0 & !GC) 
      lapply(acontrasts, function(contrast){
        groups <- unlist(strsplit(contrast, "/"))
        group1 <- groups[1]
        group2 <- groups[2]
        chi2<- (assoc[,paste("beta", group1, sep ="."), with = FALSE ] - 
                  assoc[, paste("beta", group2, sep ="."), with = FALSE])^2/
          (assoc[,paste("SE", group1, sep ="."), with = FALSE]^2 + 
             assoc[,paste("SE", group2, sep ="."), with = FALSE]^2)
        assoc[, paste("contrastpvalue", group1,group2, sep ="."):= pchisq(chi2, df = 1, lower.tail = FALSE)]})
    if(length(acontrasts)>0 & GC) {
      lapply(acontrasts, function(contrast){
        groups <- unlist(strsplit(contrast, "/"))
        group1 <- groups[1]
        group2 <- groups[2]
        file.curr<- paste("analyses/",model,"/ALL.", group1, "_vs_", group2, ".out.gc.txt.gz", sep ="" )
        message( Sys.time(), " Obtaining results from ", file.curr)
        curr <- read.table(gzfile(file.curr), quote = "", comment.char = "", header = TRUE, stringsAsFactors = FALSE)
        #SNP	pvalue.GC_LOSIR3005_vs_HISIR3005
        curr<-curr[!is.na(curr[[2]]) & curr$SNP %in% assoc$SNP, ]
        #curr<- data.table(curr[!is.na(curr[[2]]) & curr$SNP %in% assoc$SNP, ])
        #Need to solve the issue with data.table handling duplicated snps
        #setkey(curr, SNP)
        #curr<- curr[assoc$SNP, ] 
        #assoc[, paste("pvalue.GC.", group1,"_vs_", group2, sep =""):= curr[,2, with= FALSE ]] 
        assoc[, paste("pvalue.GC.", group1,"_vs_", group2, sep =""):= curr[["pvalue.GC"]][match(assoc$SNP, curr$SNP)]]})
    }
    
    setkey(assoc, SNP)
    message(Sys.time(), ": writing association result to ", paste("data/assoc.",model,  ".txt" , sep = ""))
    write.table(assoc[assoc$SNP %in% varlist, ], file = paste("data/assoc.",model,  ".txt" , sep = ""), sep = "\t", row.names = F, quote = F)
    message(Sys.time(), ": writing association result with flanking region to ", paste("data/assoc.",model,  ".flanking.txt" , sep = ""))
    write.table(assoc, file = paste("data/assoc.",model,  ".flanking.txt" , sep = ""), sep = "\t", row.names = F, quote = F)
    
    if(plot) 
      lapply(varlist, function(snp){
        prefix.curr <-  gsub("[:*]",".", snp)
        #get LD
        ld<- read.table( paste("data/", prefix.curr, ".info.flanking.txt", sep = ""), check.names = F, as.is = T, sep ="\t", header = T)
        
        #assoc result
        chr <- unlist(strsplit(snp, ":"))[1]
        pos <- as.numeric(as.character(unlist(strsplit(snp, ":"))[2]))
        assoc.curr<- assoc[assoc$CHROM == chr & assoc$POS >= pos-flanking & assoc$POS <= pos + flanking,   ]
        assoc.curr[, R2_LD:=ld$R2_LD[match(assoc.curr$SNP, ld$SNP)]]
        assoc.curr[assoc.curr$CHROM == chr & assoc.curr$POS == pos & assoc.curr$SNP != snp,  POS:= pos + 1]
        if(nrow(assoc.curr[snp,nomatch=0])>0 &nrow(assoc.curr) > 0) {
          prefix.curr <-  paste( prefix.curr,model, sep =".")
          pdf.file <- paste("plot/",prefix.curr, ".pdf", sep = "")
          message(paste(Sys.time(), ": Plotting to", pdf.file))
          pdf(file=pdf.file, height=8.5, width=11)
          
          alleles<- unlist(assoc.curr[snp, c( "Al1", "Al2"), with = FALSE])
          freq1<- assoc.curr[snp, grep("analysed.Freq1",names(assoc.curr)), with= FALSE]
          rsq<- assoc.curr[snp, grep("analysed.Rsq",names(assoc.curr)), with = FALSE]
          col.pvals<- names(assoc.curr)[grep(ifelse(GC, "pvalue.GC", "pvalue"), names(assoc.curr))]
          #region plot
          main <- paste(model,  snp, sep =": ")
          main.sub<- paste("Alleles=", paste(alleles, collapse="/", sep =""), 
                           ", Freq1=", format(min(freq1, na.rm = T), digit = 3),"~",format(max(freq1, na.rm = T), digit = 3), 
                           ", Rsq_imp=", format(min(rsq, na.rm = T), digit = 3),"~",format(max(rsq, na.rm = T), digit = 3), 
                           sep = "")
          if(length(freq1) ==1)
            main.sub<- paste("Alleles=", paste(alleles, collapse="/", sep =""), 
                             ", Freq1=", format(freq1, digit = 3),  ", Rsq_imp=", format(rsq,  digit = 3),sep = "")
          
          par(mfcol =c(1,1),mar = c(4, 4, 4, 4) + 0.1, oma = c(6,0,6,0)) 
          regionplot.multiP (assoc.curr, chr, pos, flanking = flanking, col.r2="R2_LD",
                             col.pvals = col.pvals, 
                             plabel = gsub(".","/", gsub(ifelse(GC,"pvalue.GC.", "pvalue."),"", col.pvals, fixed = T), fixed = T), 
                             gencode = genetable, recomb=recomb,
                             main = main, main.sub = main.sub, miny = 8, crity = -log10(5e-08))
          mtext(paste("Protocol:", getOption("gtxpipe.protocol", "NA")), side = 3, line = 4, adj = 0.025, outer = T)
          mtext("Page 1 of 1", side = 3, line =4, adj = 0.975, outer = T)
          mtext(paste("Population:", "ITT"), side = 3, line = 3, adj = 0.025, outer = T) # hard coded to ITT
          mtext(paste("Data as of:", getOption("gtxpipe.date", "NA")), side = 3, line = 3, adj = 0.975, outer = T)
          mtext(getOption("gtxpipe.project", NA), side = 3, line = 0, adj = 0.5, outer = T, cex = 1.2, font = 2)
          message(paste(Sys.time(), ": Region plot saved to ", pdf.file))
          
          #pheno geno plot
          par(mfrow = c(1,length(agroups)),mar = c(4, 4, 4, 4) + 0.1, oma = c(3,1,6,0))
          if(length(agroups)>2) par(mfrow=c(ceiling(length(agroups)/2), 2)) 
          
          lapply(agroups, function(group){
            detail <- getCallDetail(gres[[group]]$model)
            data.plot<- gres[[group]]$phenotype[c(detail$pheno, snp)]
            data.plot<-data.plot[complete.cases(data.plot),]
            data.plot[[snp]] <- dose2geno(data.plot[[snp]], alleles=alleles)
            beta<- unlist(assoc.curr[snp, paste("beta", group, sep ="."),with = FALSE])
            se<-unlist(assoc.curr[snp, paste(ifelse(GC, "SE.GC","SE"), group, sep ="."),with = FALSE])
            ci<-beta+c(-1, 1) *qnorm(0.975)*se
            main <- paste(model," in ",group, ": ", snp, sep ="")
            main.sub<- paste(c(paste("Freq_",alleles[1], sep =""),  "Rsq_imp",ifelse(GC, "P.GC", "P")),  
                             format(assoc.curr[snp, paste(c("analysed.Freq1", "analysed.Rsq", 
                              ifelse(GC, "pvalue.GC","pvalue")), group, sep ="."), with = FALSE], digit = 3), 
                             sep = "=", collapse = ", ")
            if(detail$type == "qt") {
              qtplot(detail$pheno, snp, data.plot, xlab = NULL)
              main.sub<- paste(main.sub,  ", eff_",alleles[1],  "(CI)=", format(beta, digit = 2), "(",
                               paste(format(ci, digit = 2), collapse =","), ")", sep ="")
            }
            if(detail$type %in% c("binary", "ordinal","negbin")) {
              otplot(detail$pheno, snp, data.plot, xlab = NULL)
              main.sub<- paste(main.sub,  ", OR_",alleles[1],  "(CI)=", format(exp(beta), digit = 2), "(",
                               paste(format(exp(ci), digit = 2), collapse =","), ")", sep ="")
            }
            if(detail$type == "survival") {
              data.plot$srv <- Surv(as.numeric(data.plot[[detail$pheno[1]]]), as.numeric(data.plot[[detail$pheno[2]]]))
              kmplot("srv", snp, data.plot, xlab =detail$pheno[1], ulab = "Survival Probability" )
              main.sub<- paste(main.sub,  ", HR_",alleles[1],  "(CI)=", format(exp(beta), digit = 2), "(",
                               paste(format(exp(ci), digit = 2), collapse =","), ")", sep ="")
            }
            title(main, line = 2)  
            title(main.sub, line = 0.5, cex.main = 0.8, col.main = "blue")
          })
          
          mtext(paste("Protocol:", getOption("gtxpipe.protocol", "NA")), side = 3, line = 4, adj = 0.025, outer = T)
          mtext("Page 1 of 1", side = 3, line =4, adj = 0.975, outer = T)
          mtext(paste("Population:", "ITT"), side = 3, line = 3, adj = 0.025, outer = T) # hard coded to ITT
          mtext(paste("Data as of:", getOption("gtxpipe.date", "NA")), side = 3, line = 3, adj = 0.975, outer = T)
          mtext(getOption("gtxpipe.project", NA), side = 3, line = 1, adj = 0.5, outer = T, cex = 1.2, font = 2)
          if(length(acontrasts)>0){
            if(!GC) contrastP<-unlist(assoc.curr[snp, paste("contrastpvalue.", gsub("/",".",acontrasts, fixed = T), sep =""), with = FALSE])
            if(GC) contrastP<-unlist(assoc.curr[snp, paste("pvalue.GC.", gsub("/","_vs_",acontrasts, fixed = T), sep =""), with = FALSE])
            
            mtext(paste(ifelse(GC, "P.GC.", "P."),acontrasts, "=", format(contrastP, digit = 3), sep ="", collapse =", "), 
                  side = 3, line = 0, adj = 0.5, outer = T, col="blue")
          }
          message(paste(Sys.time(), ": pheno/geno plot saved to ", pdf.file))
          graphics.off()
        }else message(Sys.time(), ": ", snp, " was ignored as no results found in ", model )
      })
  }
}

padata.permodelgroup <- function(retFile, varlist, flanking,  threshold.MAF, threshold.Rsq, GC = T){
  retFile <- tryCatch(as.character(retFile)[1], error = function(e) "") # sanitize
  if (!file.exists(retFile)) stop(retFile, '" does not exist')
  ##get assoc result 
  message( Sys.time(), " Obtaining results from ", retFile)
  res<- getResults(retFile, varlist, flanking, GC = GC )
  #c("CHROM", "POS", "SNP", "pvalue", "beta", "SE", "Al1", "Al2", "analysed.Freq1", "analysed.Rsq", "MAF"
  res<- res[MAF >= threshold.MAF & analysed.Rsq >= threshold.Rsq, ]
  
  adir <- dirname(retFile) # analysis directory
  
  #read in analysis-dataset.csv for phenotypes, model (type of endpoints)
  adataf <- file.path(adir, "analysis-dataset.csv") # analysis data filename
  if (!file.exists(adataf)) stop('Analysis dataset "', adataf, '" does not exist')
  message(Sys.time()," Obtaining analysis dataset from ", adataf)
  adata0 <- scan(file.path(adir, "analysis-dataset.csv"), sep = "\n", character(0), quiet = TRUE)
  
  ## obtain model name and group name
  mm <- which(substr(adata0, 1, 28)=="# Analysis dataset for model")
  tmp<- unlist(strsplit(adata0[mm], " ", fixed = TRUE))
  mname <- gsub("\"","", tmp[grep("model", tmp)+1], fixed = T)
  gname <- gsub("\"","", tmp[grep("group", tmp)+1], fixed = T)
  rm(tmp)
  
  ## match call, error if more than 1
  mc <- which(substr(adata0, 1, 8) == '# call: ')
  stopifnot(identical(length(mc), 1L))
  qcall <- parse(text = substring(adata0[mc], 9)) # quoted call
  ## match transformations (can be any number) and build into dataframe
  mt <- which(substr(adata0, 1, 9) == '# where: ')
  if (length(mt) > 0) {
    gtxpipe.transformations <- do.call(rbind, lapply(strsplit(substring(adata0[mt], 10), ' <- ', fixed = TRUE), function(t1) {
      if (!identical(length(t1), 2L)) {
        stop('fatal error: "', t1, '" is not a transformation like "target <- fun"')
      }
      data.frame(targets = t1[1], fun = t1[2], stringsAsFactors = FALSE)
    }))
  } else {
    gtxpipe.transformations <- data.frame(targets = NULL, fun = NULL)
  }
  
  ## in read.csv should force some settings (stringsAsFactors = TRUE) just in case user options
  ## try to override
  adata <- read.csv(textConnection(adata0[substr(adata0, 1, 1) != "#"]), stringsAsFactors = TRUE)
  adata[[1]] <- as.character(adata[[1]])
  ## apply transformations
  if (nrow(gtxpipe.transformations) > 0) {
    for (idx in 1:nrow(gtxpipe.transformations)) {
      target <- gtxpipe.transformations$targets[idx]
      if (target %in% names(adata)) stop("Transformation overwrites existing variable")
      message(target, " <- ", gtxpipe.transformations$fun[idx]) # debugging message
      adata[[target]] <- eval(parse(text = gtxpipe.transformations$fun[idx]), envir = adata)
    }
  }
  
  padata<- list( name = mname,
                 group = gname, 
                 model = qcall, 
                 phenotype =adata,
                 assoc = res)
  class(padata) <- "postanalysisdata"
  return (padata)
}

##get phenotype name(s) and type from a formula call
getCallDetail<- function(qcall) {
  #qcall: coxph(formula = Surv(Pheno.SRVMO, Pheno.SRVCFLCD) ~ demo.AGE + demo.SEX + PC1)
  #qcall: glm(formula = Pheno.altcc ~ pop.TRTGRP + PC1, family = "binomial")
  tmp<- unlist(strsplit(as.character(qcall), split="[(+, =)\"]"))
  tmp<- tmp[nchar(tmp)>0]
  mtype <- tmp[1]
  if(length(grep("family", tmp))>0)
    mtype <- paste(mtype, tmp[grep("family", tmp)+1], sep = ".")
  
  pheno <- tmp[grep("~", tmp)-1]
  if(mtype == "coxph") pheno <- tmp[(grep("Surv", tmp)+1): (grep("~", tmp)-1)] 
  
  ptype.list <- c(lm = "qt",  
                  clm="ordinal", 
                  coxph="survival", 
                  glm.binomial="binary",
                  glm.nb="negbin")
  ptype <- ptype.list[mtype]
  return (list(pheno = pheno, type = ptype ))
}

##subset association results by varlist and flanking size in bp
getResults<- function(retFile, varlist, flanking, GC = T ){
  ## Reading results which were compiled across chunks during make call
  res1 <- read.table(gzfile(retFile), quote = "", comment.char = "", header = TRUE, stringsAsFactors = FALSE)
  ## Confirm expected columns present
  #CHROM	POS	SNP	beta.PBO3002	SE.PBO3002	pvalue.PBO3002	Al1	Al2	analysed.Freq1.PBO3002	analysed.Rsq.PBO3002	pvalue.GC.PBO3002	SE.GC.PBO3002
  if(GC){
    #CHROM	POS	SNP	beta.PBO3002	SE.PBO3002	pvalue.PBO3002	Al1	Al2	analysed.Freq1.PBO3002	analysed.Rsq.PBO3002	pvalue.GC.PBO3002	SE.GC.PBO3002
    i<- grep("beta",  names(res1))
    stopifnot(length(i)==1)
    g<- gsub("beta", "", names(res1)[i], fixed = T)
    if(nchar(g)> 0) names(res1)<- gsub(g, "", names(res1), fixed = T)
  }
  names(res1) <- gsub("X.", "", names(res1), fixed = T) 
  stopifnot(all(c("CHROM", "POS", "SNP", "pvalue", "beta", "SE", "Al1", "Al2", "analysed.Freq1", "analysed.Rsq") %in% names(res1)))
  
  ## Convert to data table to improve efficiency retaining only needed columns and rows with non-missing pvalue
  res1 <- data.table(res1[!is.na(res1$pvalue), ])
  res1[, chrpos:=paste(CHROM, POS, sep =":")]
  ## Sort by 
  setkey(res1, chrpos)
  
  varlist.expanded <- unique(unlist(lapply(varlist, function(snp){
    chrpos<- unlist(strsplit(as.character(snp), ":", fixed = TRUE))[1:2]
    return (paste(chrpos[1], as.numeric(as.character(chrpos[2]))+ c((-1*flanking):flanking), sep =":"))})))
  varlist.expanded<- varlist.expanded[varlist.expanded %in% res1$chrpos]
  res1<- res1[varlist.expanded,]
  res1<- res1[!is.na(res1$pvalue), ]
  res1[, MAF:= ifelse(analysed.Freq1> 0.5, 1-analysed.Freq1, analysed.Freq1) ]
  setkey(res1, SNP)
  return (res1[, chrpos:=NULL ])
}


##Obtain list of chunks (prefix) for a region 
getchunks<- function(chr, pos.start, pos.end, chunkMb, chrBeginEndDir){
  if(toupper(chr) == "X") chr = 23
  chrbegin<- read.table(file.path(chrBeginEndDir, "chrbegin.data"), as.is = T)[1:2]
  chrend<- read.table(file.path(chrBeginEndDir, "chrend.data"), as.is = T)[1:2]
  chr.start <- chrbegin[[2]][chrbegin[[1]]==chr]
  chr.end <- chrend[[2]][chrend[[1]]==chr]
  N.chunk<- ceiling((chr.end-chr.start+1)/(chunkMb*1e6))
  chunks <- data.frame(chunk = paste("chr", chr, "chunk", 1:N.chunk,  sep =""), 
                       start = (0:(N.chunk-1)) * chunkMb *1e6 + chr.start , 
                       end = (1:N.chunk) * chunkMb *1e6 + chr.start - 1, 
                       stringsAsFactors =FALSE)
  chunks.idx<- intersect(which(pos.start <= chunks$end), which(pos.end >= chunks$start))
  ret<- c(chunks$chunk[chunks.idx], paste("chr", chr, "chunkG", sep =""))
  if(chr == 23) ret<- c(gsub("chr23","chr23-F",  ret), gsub("chr23","chr23-M",  ret))
  #HLA region
  if(chr==6 ) ret<- c(ret, "Imputed_HLAalleles_AllSubjects_Additive" )
  return ( ret )
}

#dose data for snp
#LD in flanking region for snp
getDoseandLD<- function(snp, flanking, ld.subj, genodir, chunkMb, chrBeginEndDir) {
  chr <- unlist(strsplit(snp, ":"))[1]
  pos <- as.numeric(as.character(unlist(strsplit(snp, ":"))[2]))
  pos.start<- pos-flanking
  pos.end<- pos + flanking
  chunks.curr<- getchunks(chr, pos.start, pos.end, chunkMb, chrBeginEndDir)
  
  dose<- NULL
  info<- NULL
  for ( prefix in chunks.curr) {        
    file.info.curr<- paste(genodir, paste(prefix, "info", sep = "."), sep="/")
    if(!file.exists(file.info.curr))
      file.info.curr<- paste(genodir, paste(prefix, "info.gz", sep = "."), sep="/")
    if(!file.exists(file.info.curr)) {
      message(prefix, " is ignored as info file is missing")
      next
    }
    message(Sys.time(), ": Loading info data from ",file.info.curr)
    if(length(grep(".gz", file.info.curr))> 0)
      info.curr <- read.table(gzfile(file.info.curr), comment = "", quote = "", header = TRUE, as.is = TRUE)
    else info.curr <- read.table(file.info.curr, comment = "", quote = "", header = TRUE, as.is = TRUE)
    col.names <- c("CHROM", "POS", "Al1")
    if(length(grep("HLAalleles",prefix))>0) col.names <-c("CHROM", "POS", "HLAGENE","Al1" )
    chrpos <- read.table(text=as.character(info.curr$SNP), sep = ":", header = F, fill = T, col.names= col.names)[1:2]
    info.curr<- cbind(info.curr, chrpos)
    markers.flag<- info.curr$POS >= pos.start & info.curr$POS <= pos.end
    info.curr<- info.curr[markers.flag, ]
    
    if(nrow(info.curr)== 0) next
    
    file.dose.curr <- paste(genodir, paste(prefix, "dose.gz", sep = "."),sep="/")
    message(Sys.time(), ": Loading dosage data from ",file.dose.curr)
    dose.curr <- read.table(gzfile(file.dose.curr), comment = "", quote = "", header = FALSE,
                            colClasses = c(rep("character", 2), rep("numeric", nrow(info.curr))))
    dose.curr <- dose.curr[c(T, F, markers.flag)]
    message(Sys.time(), ": ", nrow(info.curr), " variants retained")
    names(dose.curr) <- c("USUBJID", as.character(info.curr$SNP))
    dose.curr$USUBJID <-unlist(lapply(dose.curr$USUBJID, function(data){
      return (unlist(strsplit(data, "->", fixed = T))[1])}))     
    if(is.null(dose))   dose<- dose.curr
    else dose <- cbind(dose,dose.curr[-1])
    info<- rbind(info, info.curr)
  }
  
  dose.snp<- dose[c(1, match(snp, names(dose)))]
  
  if(!is.null(ld.subj)) 
    dose<- dose[dose$USUBJID %in% ld.subj,]
  info$R2_LD<- unlist(lapply(info$SNP, function(asnp){
    return (cor(as.numeric(as.character(dose[[asnp]])), 
                as.numeric(as.character(dose[[snp]])),
                use = "pairwise.complete.obs")^2)}))
  info$IndexSNP_LD<- snp
  return (list(dose=dose.snp, LD = info))
}



#convert dosage to genotype, 0 for allele2/2, 2 for allele 1/1
dose2geno <- function(dose, alleles=c("A", "T")) {
  genos <- c(paste(alleles[2], alleles[2], sep = ""),
             ifelse(alleles[1] > alleles[2], paste(alleles[2], alleles[1], sep = ""),paste(alleles[1], alleles[2], sep = "")),
             paste(alleles[1], alleles[1], sep = ""))
  ret <- factor(genos[ifelse(dose >1.5, 2, ifelse(dose >0.5, 1, 0))+1], levels = genos)
  return (ret)
}


