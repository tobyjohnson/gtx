pipeline <- function(config = "config.R", slave) {

  ## Most options are supplied in a configuration file because:
  ## * This allows the same options to be easily passed to slave processes initiated through a parallel make
  ## * Required packages (e.g. survival, ordinal) can be loaded (instead of passing an option of package names)
  ## * "Complex" options can be set using small blocks of code instead of
  ##   unwieldy function arguments like x = {data(xx); c(xx, 1:4)}
  if (!file.exists(config)) stop ('Configuration file "', config, '" does not exist')
  source(config)

  usubjid <- getOption("clinical.usubjid", "USUBJID")
  
  genotype.path <- tryCatch(genotype.path, error = function(e) 'genotypes')
  
  if (!missing(slave)) {
    slave <- tryCatch(as.character(slave)[1], error = function(e) "") # sanitize
    adir <- dirname(slave) # analysis directory
    
    tmp <- unlist(strsplit(basename(slave), ".", fixed = TRUE))
    if (length(tmp) < 2 || !identical(tmp[length(tmp)], "done")) stop ('Bad slave "', slave, '"')
    job <- paste(tmp[-length(tmp)], collapse = ".") # job is basename with trailing .done stripped off
    rm(tmp)
    
    adataf <- file.path(adir, "analysis-dataset.csv") # analysis data filename
    if (!file.exists(adataf)) stop('Analysis dataset "', adataf, '" does not exist')
    adata0 <- scan(file.path(adir, "analysis-dataset.csv"), sep = "\n", character(0), quiet = TRUE)
    mc <- which(substr(adata0, 1, 8) == '# call: ') # match call, error if more than 1
    stopifnot(identical(length(mc), 1L))
    qcall <- parse(text = substr(adata0[mc], 9, nchar(adata0[mc]))) # quoted call
    ## in read.csv should force some settings (stringsAsFactors = TRUE) just in case user options
    ## try to override
    adata <- read.csv(textConnection(adata0[substr(adata0, 1, 1) != "#"]),
                      stringsAsFactors = TRUE)
    adata[[usubjid]] <- as.character(adata[[usubjid]])
    ## ? could sink to .done file, blockassoc should use message() not cat()
    res <- blockassoc(qcall = qcall, data = adata,
                      minimac = file.path(genotype.path, job),
                      threshold.MAF = .01, threshold.Rsq = 0)
    ## Require pvalue always, beta and SE if being used for interaction contrasts
    ## Consider writing some comments to this file (blockassoc should return character vector instead of printing messages)
    write.table(res[ , c("SNP", "beta", "SE", "pvalue")], # should be an option to set outputs, esp analysed.*
                row.names = FALSE, quote = FALSE, 
                file = gzfile(file.path(adir, paste(job, "out.gz", sep = "."))))
    return(invisible(NULL))
  }
        
  ## Not in slave mode

  pgx.models <- tryCatch(pgx.models, error = function(e) data.frame(NULL))
  if (nrow(pgx.models) < 1) stop("No PGx models specified.  Set pgx.models")
  ## check pgx.models is a data frame with required columns
  ## check pgx.models$model are valid directory names
  
  pgx.groups <- tryCatch(pgx.groups, error = function(e) data.frame(NULL))
  if (nrow(pgx.models) < 1) {
    warning("No PGx groups specified.  Using a single (ITT) group.  Set pgx.groups")
    pgx.groups <- data.frame(group = 'ITT', deps = 'pop.PNITT', fun = 'pop.PNITT == "Y"', stringsAsFactors = FALSE)
  }

  pgx.transformations <- tryCatch(pgx.transformations, error = function(e) data.frame(NULL))
  
  ## Parse and check groups and constrasts
  ## Compute multiple testing adjustment (total number of groups and contrasts)
  ## Compute total set of groups that need within-group analyses
  pgx.models <- cbind(pgx.models, 
                      do.call(rbind, lapply(1:nrow(pgx.models), function(modelid) {
                        contrasts1 <- tokenise.whitespace(pgx.models[modelid, "contrasts"])
                        contrasts1.bad <- sapply(contrasts1, function(contrast1) return(length(unlist(strsplit(contrast1, "/"))) != 2))
                        if (any(contrasts1.bad)) {
                          stop('PGx model "', pgx.models[modelid, "model"], '" has contrasts(s) ',
                               paste('"', contrasts1[contrasts1.bad], '"', sep = "", collapse = ", "),
                               ' not in format group1/group2')
                        }
                        groups1 <- tokenise.whitespace(pgx.models[modelid, "groups"])
                        agroups <- unique(c(groups1, unlist(strsplit(contrasts1, "/"))))
                        agroups.missing <- is.na(match(agroups, pgx.groups$group))
                        if (any(agroups.missing)) {
                          stop('PGx model "', pgx.models[modelid, "model"], '" has undefined group(s) ',
                               paste('"', agroups[agroups.missing], '"', sep = "", collapse = ", "))
                        }
                        return(data.frame(alpha1 = length(contrasts1) + length(groups1), 
                                          agroups = paste(agroups, collapse = " "), 
                                          stringsAsFactors = FALSE))
                      })))

  ## The set of all variables we need for actual analysis
  ## Transformations -> new deps

  ## Compute dependencies in terms of untransformed variables
  pgx.models$depsu <- sapply(pgx.models$deps, function(deps1) {
    paste(sapply(tokenise.whitespace(deps1), function(dep1) {
      mm <- match(dep1, pgx.transformations$targets)
      if (is.na(mm)) return(dep1)
      return(pgx.transformations$deps[mm])
    }), collapse = " ")})

  deps <- unique(tokenise.whitespace(c(pgx.models$depsu, pgx.groups$deps,
                                       "pop.TRTGRP", "demo.SEX", "demo.AGE", "demo.RACE", "demo.ETHNIC")))
  ## force in pop.PNITT even though this is not a mandated variable per dsm ?
  ## allow a force in list.  sort by unique(forcelist, deps)
  
  ## check pgx.groups is a data frame with required columns (group, deps, fun) of character class
  ## check pgx.groups$group are valid directory names
  clinical.derivations <- tryCatch(clinical.derivations, error = function(e) data.frame(NULL))

  ## which clindata do we actually need?
  ddeps <- unique(tokenise.whitespace(clinical.derivations[sapply(1:nrow(clinical.derivations), function(idx) {
    any(tokenise.whitespace(clinical.derivations[idx, "targets"]) %in% deps)
  }), "deps"]))

  cat("Require clinical datasets: ", paste(ddeps, collapse = ", "), "\n")
  
  clinical.path <- tryCatch(clinical.path, error = function(e) 'clinical')
  clindata <- clinical.import(clinical.path, only = ddeps)
  ### Check clinical datasets required were all read / import should have only and require options??

  message("Read clinical datasets OK")
  anal1 <- clinical.derive(clindata, clinical.derivations, only = deps)
  message("Computed derived variables OK")

#  inc1 <- with(anal1, eval(parse(text = getOption("clinical.subset", 'pop.PNITT == "Y"'))))

  ancestrypcs <- read.table(pgx.eigenvec, header = FALSE, sep = " ", as.is = TRUE)[ , -1] # drop first (FID) column
  names(ancestrypcs) <- c(usubjid, paste("PC", 1:(ncol(ancestrypcs) - 1), sep = ""))
  
  anal1$is.PGx <- anal1[[usubjid]] %in% ancestrypcs[[usubjid]]


  
  groupby <- do.call(c, lapply(1:nrow(pgx.groups), function(idx) {
    return(with(anal1, list(eval(parse(text = pgx.groups$fun[idx])),
                            eval(parse(text = pgx.groups$fun[idx])) & is.PGx)))
  }))
  
  names(groupby) <- paste(rep(pgx.groups$group, each = 2),
                          rep(c("ITT", "PGx"), times = nrow(pgx.groups)),
                          sep = ", ")

  pgx.groups$N.ITT <- sapply(groupby[paste(pgx.groups$group, "ITT", sep = ", ")], sum)
  pgx.groups$N.PGx <- sapply(groupby[paste(pgx.groups$group, "PGx", sep = ", ")], sum)

## Automagic report what arms (pop.TRTGRP groups) are actually included 
  pgx.groups <- cbind(pgx.groups, 
    do.call(rbind, lapply(groupby[paste(pgx.groups$group, "ITT", sep = ", ")], function(gl) {
      tmp <- table(anal1$pop.TRTGRP[gl])
      tmp <- tmp[tmp > 0]
      if (length(tmp) == 0) return(data.frame(arms = "None", adjust.arm = NA, stringsAsFactors = FALSE))
      if (length(tmp) == 1) return(data.frame(arms = tmp, adjust.arm = FALSE, stringsAsFactors = FALSE))
      return(data.frame(arms = paste(names(tmp), " (N=", tmp, ")", sep = "", collapse = ", "),
                        adjust.arm = TRUE, stringsAsFactors = FALSE))
    })))

  pgx.groups$N.notPGx <-   pgx.groups$N.ITT - pgx.groups$N.PGx
  
  ## write into source subdir
  pgx.groups[ , c("group", "arms", "N.ITT", "N.PGx")] # source item 1, can we call this disposition?
  # Percent PGx
  # apply(pgx.groups[ , c("N.notPGx", "N.PGx")], 1, prettypc)[2, , drop = TRUE]

  any(rowSums(as.data.frame(groupby)) > 1) # TRUE => overlapping groups
  ## rowSums all 0 or 1, no overlapping groups.  Any >1 implies overlapping groups
  ## coercing logical to integer

  ## allow a demographics sort list
  ## write into source subdir
  write.csv(demographics(anal1, by = groupby), file = "demographics.csv") # source item 2

  adir0 <- "analyses" # this should be an option

  sink("Makefile") # this should be an option

  cat("## Makefile generated by gtx pipeline()\n")
  cat(".PHONY:\tall\nall:\t\n\n") # check having all twice is okay
  alltarget <- c(NULL)
  
  cat("GENOTYPES := $(wildcard ", genotype.path, "/*.dose.gz)\n\n", sep = "")

## Write analysis datasets and Makefile
  for (modelid in 1:nrow(pgx.models)) {
    for (agroup1 in tokenise.whitespace(pgx.models[modelid, "agroups"])) {
      groupid <- match(agroup1, pgx.groups$group)
      trtgrp.cov <- if (pgx.groups$adjust.arm[groupid]) "pop.TRTGRP" else NULL
      adir <- file.path(adir0, pgx.models[modelid, "model"], pgx.groups[groupid, "group"])
      dir.create(adir, recursive = TRUE, showWarnings = FALSE)
      adata <- merge(anal1[groupby[[paste(agroup1, ", PGx", sep = "")]],
                           c(usubjid,  
                             tokenise.whitespace(pgx.models[modelid, "depsu"]),
                             trtgrp.cov)],
                     ancestrypcs,
                     all = FALSE)
      ## Remove subjects with any missing data on any dependency variable
      adata <- adata[which(apply(!is.na(adata), 1, all)), ]
      ## Apply transformations to *the analysis dataset*
      ## ONLY APPLY RELEVANT TRANSFORMATIONS depsu -> deps
      ## This should be done by slave...
      for (idx in 1:nrow(pgx.transformations)) {
        target <- pgx.transformations$targets[idx]
        if (target %in% names(adata)) stop("Transformation overwrites existing variable")
        adata[[target]] <- eval(parse(text = pgx.transformations$fun[idx]), envir = adata)
      }
      adata <- adata[ , unique(c(usubjid, 
                                 tokenise.whitespace(pgx.models[modelid, "deps"]),
                                 trtgrp.cov,
                                 names(ancestrypcs)))]

      sel <- apply(!is.na(adata), 1, all)
      adata <- adata[sel, , drop = FALSE]

      ## Analysis datasets will not preserve R classes such as Surv
      ## Everything except usubjid will be a factor

      ## Everything as paths relative to getwd() at runtime

      mtmp <- eval(parse(text = pgx.models[modelid, "fun"]), envir = adata)
      m0 <- update(mtmp, formula = as.formula(paste("~ . +", paste(c(trtgrp.cov, names(ancestrypcs)[-1]), collapse = "+"))))
      drop1(m0, test="Chisq")

      ## Should test and not overwrite, or clear out analyses uf updating
      
      sink(file.path(adir, "analysis-dataset.csv")) # sink inside sink
      cat('# Analysis dataset for model "', pgx.models[modelid, "model"], '" in group "', pgx.groups[groupid, "group"], '"\n', sep = '')
      ## first 8 characters MUST be '# call : '
      cat('# call: ', as.character(as.expression(m0$call)), '\n', sep = '')
      ## consider adding transformations like '# where: osMonths <- Surv(SRVMO, SRVCFLCD)'

      sink()
      suppressWarnings(write.table(adata, sep = ",", row.names = FALSE, 
                                   file = file.path(adir, "analysis-dataset.csv"),
                                   append = TRUE))
      
      cat('# Analysis for model "', pgx.models[modelid, "model"], '" in group "', pgx.groups[groupid, "group"], '"\n', sep = '')
      cat('MODEL', modelid, 'GROUP', groupid, ' := $(patsubst ', genotype.path, '/%.dose.gz,', adir, '/%.done,$(GENOTYPES))\n', sep = '')
      cat('$(MODEL', modelid, 'GROUP', groupid, '):\t', adir, '/%.done:\t', genotype.path, '/%.info.gz ', genotype.path, '/%.dose.gz\n', sep = '')
      cat('\t@(echo "library(gtx); pipeline(config = \\"', config, '\\", slave = \\"$@\\")" | R --vanilla && touch $@)\n\n', sep = '')
      ## delete output, delete done file, run R, touch done file
      ## 		mkdir -p a1; rm -f $@; sleep 60; uname -a >$@; date >>$@
      alltarget <- c(alltarget, paste('$(MODEL', modelid, 'GROUP', groupid, ')', sep = ''))
    }
  }

  cat("all:\t", paste(alltarget, collapse = " "), "\n", sep = "") 
  sink() # Makefile

  ## Call "make all" using option for make command
  ## SGE_ARCH=lx24-amd64 qmake -v PATH -cwd -l qname=dl580 -- --jobs=4
  ## SGE_ARCH=lx24-amd64 nohup qmake -v PATH -cwd -l qname=dl580 -- --jobs=256 &
  
  ## consider adding -l h_data=4G

  message("Run make now")


  ## index
  ## p-values for noncontrasts, GC, 
  ## contrasts, assume independent (user responsibility), b1-b2, sqrt(s12+s22), pval, gc
  ## for all, number non-NA, number <= 5e-8/alpha/2

  ## candidate variants, number <.05/K/alpha/2
  ## DONE!

  ## (no need to read info files if no hits...)

  ## Default filtering on Rsq (at e.g. 0.05) would eliminate invariant dosages 0.3,...,0.3,...
  return(invisible(NULL))
}
