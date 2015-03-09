## Using contrasts as a variable name shadows a base function !!

## option namespace should all be
## gtx.* or gtxpipe.*

## pgx.eigenvec is actually list of genotyped subjects plus any covariates desired

gtxpipe <- function(gtxpipe.models = getOption("gtxpipe.models"),
                    gtxpipe.groups = getOption("gtxpipe.groups"),
                    gtxpipe.derivations = getOption("gtxpipe.derivations"),
                    gtxpipe.transformations = getOption("gtxpipe.transformations", data.frame(NULL)),
                    gtxpipe.eigenvec) {
  ## arguments for project-specific


  message("gtxpipe() from package gtx version ", packageVersion("gtx"), " on ", R.version.string)
  
  usubjid <- as.character(getOption("gtx.usubjid", "USUBJID"))[1] # variable name for unique subject identifier
  ## Replace this with a function getusubjid that prints warning messages, applies make.names() etc
  
  ## options where read repeatedly in code, constant (in C const sense)
  ## and likely constant over different projects on same IT
  ## overwrite with failsafe defaults now
  options(gtxpipe.genotypes = as.character(getOption("gtxpipe.genotypes", "genotypes"))[1]) # directory containing genotype data
  options(gtxpipe.clinical = as.character(getOption("gtxpipe.clinical", "clinical"))[1]) # directory containing clinical data
  options(gtxpipe.analyses = as.character(getOption("gtxpipe.analyses", "analyses"))[1]) # top level directory for analyses
  options(gtxpipe.outputs = as.character(getOption("gtxpipe.outputs", "outputs"))[1]) # target directory for outputs
  dir.create(getOption("gtxpipe.outputs"), recursive = TRUE, showWarnings = FALSE) # throw error?

  ##
  ## Check gtxpipe.models
  ##
  if (missing(gtxpipe.models) || is.null(gtxpipe.models)) stop("No models specified.  You need to set gtxpipe.models")
  gtxpipe.models <- as.data.frame(gtxpipe.models)
  if (nrow(gtxpipe.models) < 1) stop("No models specified.  You need to have at least one row in gtxpipe.models")
  with(list(mn = setdiff(c("model", "deps", "fun"), names(gtxpipe.models))), 
       if (length(mn) > 0) stop("Models not correctly specified.  You need to have columns ",
                                paste(mn, collapse = ", "), " in gtxpipe.models"))
  gtxpipe.models$model <- gsub("\\s+", "", gtxpipe.models$model)
  ## Enforce alphanumeric names because they will be used as directory names
  with(list(bm = !grepl("^[A-Za-z0-9]+$", gtxpipe.models$model)),
       if (any(bm)) stop("Models ", paste(models[bm, "model"], collapse = ", "),
                         "have non-alphanumeric names.  You need to fix this in gtxpipe.models"))
  if (is.na(match("groups", names(gtxpipe.models)))) {
    warning("Models do not have analysis group(s) specified.  Defaulting to ITT")
    gtxpipe.models$groups <- "ITT"
    ## Note this group is always defined.  Should not always be called "ITT" though.  FIXME: Generalize all subjects concept.
  }
  if (is.na(match("contrasts", names(gtxpipe.models)))) {
    warning("Models do not have analysis contrast(s) specified.  Defaulting to none")
    gtxpipe.models$contrasts <- ""
  }
  if (is.na(match("cvlist", names(gtxpipe.models)))) {
    warning("Models do not have candidate variant(s) specified.  Defaulting to none")
    gtxpipe.models$cvlist <- ""
  }
  # FIXME: check column classes are all character?
  
  ##
  ## Check gtxpipe.groups
  ##
  if (missing(gtxpipe.groups) || is.null(gtxpipe.groups)) {
    warning("No groups specified.  Defaulting to a single (ITT) group.")
    gtxpipe.groups <- data.frame(group = 'ITT', deps = 'pop.PNITT', fun = 'pop.PNITT == "Y"', stringsAsFactors = FALSE)
  }
  gtxpipe.groups <- as.data.frame(gtxpipe.groups)
  if (nrow(gtxpipe.groups) < 1) stop("No groups specified.  You need to have at least one row in gtxpipe.groups")
  with(list(mn = setdiff(c("group", "deps", "fun"), names(gtxpipe.groups))), 
       if (length(mn) > 0) stop("Groups not correctly specified.  You need to have columns ",
                                paste(mn, collapse = ", "), " in gtxpipe.groups"))
  gtxpipe.groups$group <- gsub("\\s+", "", gtxpipe.groups$group)
  ## Enforce alphanumeric names because they will be used as directory names
  with(list(bg = !grepl("^[A-Za-z0-9]+$", gtxpipe.groups$group)),
       if (any(bg)) stop("Groups ", paste(groups[bm, "group"], collapse = ", "),
                         "have non-alphanumeric names.  You need to fix this in gtxpipe.groups"))
  ## Could we robustly assume that enrolment corresponds to *some* flag pop.PNXXXX == "Y" for some value of XXXX
  ## If no ITT group (assumed to be all subjects enrolled) is defined, add one
  if (is.na(match("ITT", gtxpipe.groups$group))) {
    gtxpipe.groups <- rbind(gtxpipe.groups,
                            data.frame(group = 'ITT', deps = 'pop.PNITT', fun = 'pop.PNITT == "Y"', stringsAsFactors = FALSE))
  }
  
  ##
  ## Check gtxpipe.derivations
  ##
  if (missing(gtxpipe.derivations) || is.null(gtxpipe.derivations)) stop("No derivations specified.  You need to set gtxpipe.derivations")
  gtxpipe.derivations <- as.data.frame(gtxpipe.derivations)  
  with(list(mn = setdiff(c("targets", "types", "deps", "data", "fun"), names(gtxpipe.derivations))), 
       if (length(mn) > 0) stop("Derivations not correctly specified.  You need to have columns ",
                                paste(mn, collapse = ", "), " in gtxpipe.derivations"))

  ##
  ## Check gtxpipe.transformations
  ##
  if (missing(gtxpipe.transformations)) gtxpipe.transformations <- NULL
  gtxpipe.transformations <- as.data.frame(gtxpipe.transformations)
  ## FIXME check code runs with no transformations
  
  ##
  ## Check groups and contrasts specified in gtxpipe.models are defined in gtxpipe.groups
  ## Compute modelwise multiple testing adjustment, alpha1=1/(total number of groups and contrasts)
  ## Compute analysis groups (agroups), the set of groups that require within-group analyses
  ##  e.g. a group not specified for a primary analysis but required for a contrast
  gtxpipe.models <- cbind(gtxpipe.models, 
                      do.call(rbind, lapply(1:nrow(gtxpipe.models), function(modelid) {
                        contrasts1 <- tokenise.whitespace(gtxpipe.models[modelid, "contrasts"])
                        ## make these local variables?
                        contrasts1.bad <- sapply(contrasts1, function(contrast1) return(length(unlist(strsplit(contrast1, "/"))) != 2))
                        if (any(contrasts1.bad)) {
                          stop('Model "', gtxpipe.models[modelid, "model"], '" has contrasts(s) ',
                               paste(contrasts1[contrasts1.bad], collapse = ", "),
                               ' not in format group1/group2')
                        }
                        groups1 <- tokenise.whitespace(gtxpipe.models[modelid, "groups"])
                        agroups <- unique(c(groups1, unlist(strsplit(contrasts1, "/"))))
                        agroups.missing <- is.na(match(agroups, gtxpipe.groups$group))
                        if (any(agroups.missing)) {
                          stop('Model "', gtxpipe.models[modelid, "model"], '" has undefined group(s) ',
                               paste('"', agroups[agroups.missing], '"', sep = "", collapse = ", "))
                        }
                        return(data.frame(alpha1 = 1/(length(contrasts1) + length(groups1)), 
                                          agroups = paste(agroups, collapse = " "), 
                                          stringsAsFactors = FALSE))
                      })))

  ## The set of all variables we need for actual analysis
  ## Transformations -> new deps

  ## Compute model dependencies in terms of untransformed variables
  gtxpipe.models$depsu <- sapply(gtxpipe.models$deps, function(deps1) {
    paste(sapply(tokenise.whitespace(deps1), function(dep1) {
      mm <- match(dep1, gtxpipe.transformations$targets)
      if (is.na(mm)) return(dep1)
      return(gtxpipe.transformations$deps[mm])
    }), collapse = " ")})

  deps <- unique(tokenise.whitespace(c(gtxpipe.models$depsu, gtxpipe.groups$deps,
                                       "pop.TRTGRP", "demo.SEX", "demo.AGE", "demo.RACE", "demo.ETHNIC")))
  ## force in pop.PNITT even though this is not a mandated variable per dsm ?
  ## allow a force in list.  sort by unique(forcelist, deps)

  with(list(bd = setdiff(deps, tokenise.whitespace(gtxpipe.derivations$targets))),
       if (length(bd) > 0) stop("Required variable(s) ", paste(bd, collapse = ", "),
                                " have no derivation rule(s) defined.  You need to add rules to gtxpipe.derivations"))
  message("gtxpipe: Required variables: ", paste(deps, collapse = ", "))
  
  ## which clinical datasets do we actually need?  clinical.derive *always* requires a pop dataset
  ddeps <- unique(c("pop", tokenise.whitespace(gtxpipe.derivations[sapply(1:nrow(gtxpipe.derivations), function(idx) {
    any(tokenise.whitespace(gtxpipe.derivations[idx, "targets"]) %in% deps)
  }), "deps"])))

  message("gtxpipe: Required datasets: ", paste(ddeps, collapse = ", "))

  message('gtxpipe: Looking for datasets in gtxpipe.clinical="', getOption("gtxpipe.clinical"), '"')
  clindata <- clinical.import(getOption("gtxpipe.clinical"), only = ddeps)
  
  with(list(bd = setdiff(ddeps, names(clindata))),
       if (length(bd) > 0) stop("Required dataset(s) not found"))
  ## import should have only and require options??
  
  #message("gtxpipe: Read clinical datasets OK")
  anal1 <- clinical.derive(clindata, gtxpipe.derivations, only = deps)
  message("gtxpipe: Computed derived variables OK")

#  inc1 <- with(anal1, eval(parse(text = getOption("clinical.subset", 'pop.PNITT == "Y"'))))

  if (!missing(gtxpipe.eigenvec) && file.exists(gtxpipe.eigenvec)) {
    ancestrypcs <- read.table(gtxpipe.eigenvec, header = FALSE, sep = " ", as.is = TRUE)[ , -1] # drop first (FID) column
    names(ancestrypcs) <- c(usubjid, paste("PC", 1:(ncol(ancestrypcs) - 1), sep = ""))
  } else {
    warning("PGx eigenvectors not available.  Continuing imperfectly")
    ancestrypcs <- anal1[ , usubjid, drop = FALSE]
  }
  ## is.PGx indicates PGx only if gtxpipe.eigenvec is one-to-one list of PGx subjects
  anal1$is.PGx <- anal1[[usubjid]] %in% ancestrypcs[[usubjid]]
  
  groupby <- do.call(c, lapply(1:nrow(gtxpipe.groups), function(idx) {
    return(with(anal1, list(eval(parse(text = gtxpipe.groups$fun[idx])),
                            eval(parse(text = gtxpipe.groups$fun[idx])) & is.PGx)))
  }))
  
  names(groupby) <- paste(rep(gtxpipe.groups$group, each = 2),
                          rep(c("ITT", "PGx"), times = nrow(gtxpipe.groups)),
                          sep = ", ")

  gtxpipe.groups$N.ITT <- sapply(groupby[paste(gtxpipe.groups$group, "ITT", sep = ", ")], sum)
  gtxpipe.groups$N.PGx <- sapply(groupby[paste(gtxpipe.groups$group, "PGx", sep = ", ")], sum)

  ## Automagic report what arms (pop.TRTGRP groups) are actually included 
  gtxpipe.groups <- cbind(gtxpipe.groups, 
    do.call(rbind, lapply(groupby[paste(gtxpipe.groups$group, "ITT", sep = ", ")], function(gl) {
      tmp <- table(anal1$pop.TRTGRP[gl])
      tmp <- tmp[tmp > 0]
      if (length(tmp) == 0) return(data.frame(arms = "None", adjust.arm = NA, stringsAsFactors = FALSE))
      if (length(tmp) == 1) return(data.frame(arms = tmp, adjust.arm = FALSE, stringsAsFactors = FALSE))
      return(data.frame(arms = paste(names(tmp), " (N=", tmp, ")", sep = "", collapse = ", "),
                        adjust.arm = TRUE, stringsAsFactors = FALSE))
    })))

  gtxpipe.groups$N.notPGx <-   gtxpipe.groups$N.ITT - gtxpipe.groups$N.PGx

  ## FIXME: Compute ACTUAL number enrolled,
  ## user can do this by making an ITT group
  ## but gtxpipe needs to know for the purpose of automagic report.
  ## Assume there will always be a group called "ITT" == enrolled???

  
  ## source item 1, can we call this disposition?
  ## write into source subdir
  print(gtxpipe.groups[ , c("group", "arms", "N.ITT", "N.PGx")])
  message('gtxpipe: Writing subject disposition to "', file.path(getOption("gtxpipe.outputs"), "01_subject_disposition.csv"), '"')
  write.csv(gtxpipe.groups[ , c("group", "arms", "N.ITT", "N.PGx")],
              file = file.path(getOption("gtxpipe.outputs"), "01_subject_disposition.csv"),
              row.names = FALSE)
  
              
                                        # Percent PGx
  # apply(gtxpipe.groups[ , c("N.notPGx", "N.PGx")], 1, prettypc)[2, , drop = TRUE]

  any(rowSums(as.data.frame(groupby)) > 1) # TRUE => overlapping groups
  ## rowSums all 0 or 1, no overlapping groups.  Any >1 implies overlapping groups
  ## coercing logical to integer

  ## allow a demographics sort list
  ## write into source subdir
  message('gtxpipe: Writing subject demographics to "', file.path(getOption("gtxpipe.outputs"), "02_subject_demographics.csv"), '"')
  write.csv(demographics(anal1, by = groupby),
            file = file.path(getOption("gtxpipe.outputs"), "02_subject_demographics.csv"),
            row.names = FALSE)

  ## There's a dataframe (here) and an option with the same name - FIXME
  gtxpipe.analyses <- do.call(rbind, lapply(1:nrow(gtxpipe.models), function(modelid) {
    agroups <- tokenise.whitespace(gtxpipe.models[modelid, "agroups"])
    acontrasts <- tokenise.whitespace(gtxpipe.models[modelid, "contrasts"])
    aprimary <- c(agroups %in% tokenise.whitespace(gtxpipe.models[modelid, "groups"]),
                  rep(TRUE, length(acontrasts)))
    return(data.frame(model = rep(gtxpipe.models[modelid, "model"], length(agroups) + length(acontrasts)),
                      group = c(agroups, acontrasts),
                      primary = aprimary,
                      stringsAsFactors = FALSE))
  }))

  adir0 <- "analyses" # this should be an option

  sink("Makefile") # this should be an option

  cat("## Makefile generated by gtxpipe()\n")
  cat(".PHONY:\tall\nall:\t\n\n") # check having all twice is okay
  alltarget <- c(NULL)
  
  cat("GENOTYPES := $(wildcard ", getOption("gtxpipe.genotypes"), "/*.dose.gz)\n\n", sep = "")

  ## Loop over analyses not nested loops, set N in gtxpipe.analyses
  ## Write analysis datasets and Makefile as side effects
  for (modelid in 1:nrow(gtxpipe.models)) {
    for (agroup1 in tokenise.whitespace(gtxpipe.models[modelid, "agroups"])) {
      groupid <- match(agroup1, gtxpipe.groups$group)
      trtgrp.cov <- if (gtxpipe.groups$adjust.arm[groupid]) "pop.TRTGRP" else NULL
      adir <- file.path(adir0, gtxpipe.models[modelid, "model"], gtxpipe.groups[groupid, "group"])
      dir.create(adir, recursive = TRUE, showWarnings = FALSE)
      adata <- merge(anal1[groupby[[paste(agroup1, ", PGx", sep = "")]],
                           c(usubjid,  
                             tokenise.whitespace(gtxpipe.models[modelid, "depsu"]),
                             trtgrp.cov)],
                     ancestrypcs,
                     all = FALSE)
      ## Remove subjects with any missing data on any dependency variable
      adata <- adata[which(apply(!is.na(adata), 1, all)), ]
      ## Apply transformations to *the analysis dataset*
      ## ONLY APPLY RELEVANT TRANSFORMATIONS depsu -> deps
      ## This should be done by slave...
      for (idx in 1:nrow(gtxpipe.transformations)) {
        target <- gtxpipe.transformations$targets[idx]
        if (target %in% names(adata)) stop("Transformation overwrites existing variable")
        adata[[target]] <- eval(parse(text = gtxpipe.transformations$fun[idx]), envir = adata)
      }
      adata <- adata[ , unique(c(usubjid, 
                                 tokenise.whitespace(gtxpipe.models[modelid, "deps"]),
                                 trtgrp.cov,
                                 names(ancestrypcs)))]

      sel <- apply(!is.na(adata), 1, all)
      adata <- adata[sel, , drop = FALSE]

      ## Analysis datasets will not preserve R classes such as Surv
      ## Everything except usubjid will be a factor

      ## Everything as paths relative to getwd() at runtime

      message("Fitting model ", gtxpipe.models[modelid, "model"], " in group ", agroup1)
      mtmp <- eval(parse(text = gtxpipe.models[modelid, "fun"]), envir = adata)
      if (length(c(trtgrp.cov, names(ancestrypcs)[-1])) > 0) {
        m0 <- update(mtmp, formula = as.formula(paste("~ . +", paste(c(trtgrp.cov, names(ancestrypcs)[-1]), collapse = "+"))))
      } else {
        m0 <- mtmp
      }

      ## Print model summary in Gx report.
      ## NOTE THAT WITH EMPTY MODELS drop1() throws an unhelpful error message
      ## drop1(m0, test="Chisq")

      ## Ensure that relevant options set in this gtxpipe() call are available
      ## to slave calls
      sink(file.path(adir, "options.R")) # sink inside sink
      cat('options(gtx.usubjid = "', usubjid, '")\n', sep = '')
      cat('options(gtxpipe.genotypes = "', getOption("gtxpipe.genotypes"), '")\n', sep = '') # will break if not scalar!
      cat('options(gtxpipe.threshold.MAF = ', getOption("gtxpipe.threshold.MAF", 0.01), ')\n', sep = '')
      cat('options(gtxpipe.threshold.Rsq = ', getOption("gtxpipe.threshold.Rsq", 0.01), ')\n', sep = '')

      sink() # options.R
      
      ## Should test and not overwrite, or clear out analyses if updating

      sink(file.path(adir, "analysis-dataset.csv")) # sink inside sink
      cat('# Analysis dataset for model "', gtxpipe.models[modelid, "model"], '" in group "', gtxpipe.groups[groupid, "group"], '"\n', sep = '')
      ## first 8 characters MUST be '# call : '
      cat('# call: ', as.character(as.expression(m0$call)), '\n', sep = '')
      ## consider adding transformations like '# where: osMonths <- Surv(SRVMO, SRVCFLCD)'

      sink()
      suppressWarnings(write.table(adata, sep = ",", row.names = FALSE, 
                                   file = file.path(adir, "analysis-dataset.csv"),
                                   append = TRUE))
      
      cat('# Analysis for model "', gtxpipe.models[modelid, "model"], '" in group "', gtxpipe.groups[groupid, "group"], '"\n', sep = '')
      cat('MODEL', modelid, 'GROUP', groupid, ' := $(patsubst ', getOption("gtxpipe.genotypes"), '/%.dose.gz,', adir, '/%.done,$(GENOTYPES))\n', sep = '')
      cat('$(MODEL', modelid, 'GROUP', groupid, '):\t', adir, '/%.done:\t', getOption("gtxpipe.genotypes"), '/%.info.gz ', getOption("gtxpipe.genotypes"), '/%.dose.gz\n', sep = '')
      cat('\t@(echo "library(gtx); gtx:::pipeslave(target = \\"$@\\")" | R --vanilla && touch $@)\n\n', sep = '')
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
  ## consider a series of make comments run in series
  message("Running make now")
  makesuccess <- system(getOption("gtxpipe.make", "make"))
  if (makesuccess != 0) {
    stop("make failed")
  }

  ## make can silently fail, need to check for .done files
  ## should know how many matched dose/info file pairs from earlier
  
  gtxpipe.results <- do.call(rbind, lapply(1:nrow(gtxpipe.models), function(modelid) {
    alpha <- .05 # this should be an option
    cvlist <- if (!is.na(gtxpipe.models[modelid, "cvlist"])) tokenise.whitespace(gtxpipe.models[modelid, "cvlist"]) else NULL
    if (length(cvlist) > 0) {
                                        ## alpha spend 0.5 for GWAS, 0.5 for CV
      thresh1 <- alpha*0.5*gtxpipe.models$alpha1[modelid]*1e-6
      thresh2 <- alpha*0.5*gtxpipe.models$alpha1[modelid]/length(cvlist)
    } else {
      thresh1 <- alpha*gtxpipe.models$alpha1[modelid]*1e-6
      thresh2 <- 0
    }
    
    agroups <- tokenise.whitespace(gtxpipe.models[modelid, "agroups"])
    res <- lapply(agroups, function(agroup1) {
      message("Collating results for model ", gtxpipe.models[modelid, "model"], " in group ", agroup1)
      adir <- file.path(adir0, gtxpipe.models[modelid, "model"], agroup1)
      res1 <- rbindlist(lapply(dir(adir, pattern = ".*\\.out\\.gz"), function(rfile) {
        tmp <- read.table(gzfile(file.path(adir, rfile)),
                          quote = "", comment.char = "", header = TRUE, stringsAsFactors = FALSE)
        stopifnot(all(c("SNP", "pvalue", "beta", "SE") %in% names(tmp)))
        return(data.table(tmp[!is.na(tmp$pvalue), c("SNP", "pvalue", "beta", "SE"), drop = FALSE]))
      }))
      setkey(res1, SNP)
      ## Makes sense to apply GC and not store redundant (non-GCed) results here
      ## pvalue from 'best' method, LRT or F test
      lambda <- res1[ , gclambda(pvalue)]
      setattr(res1, "lambda", lambda)
      invisible(if (lambda > 1.) {
        res1[ , pvalue := pchisq(qchisq(pvalue, df = 1, lower.tail = FALSE)/lambda, df = 1, lower.tail = FALSE)]
        ## overwrite, save memory
      })
      setnames(res1, "pvalue", "pvalue.GC") # not named by group to facilitate later calcs
      ## pvalue from Wald test
      lambdaWald <- res1[ , gclambda(pchisq((beta/SE)^2, df = 1, lower.tail = FALSE))]
      setattr(res1, "lambdaWald", lambdaWald)
      invisible(if (lambdaWald > 1.) {
        res1[ , SE := SE*sqrt(lambdaWald)] # overwrite, save memory
      })
      setnames(res1, "SE", paste("SE.GC", agroup1, sep = "."))
      setnames(res1, "beta", paste("beta", agroup1, sep = ".")) # named by groups, used in clever join in constrasts
      return(res1)
    })
    names(res) <- agroups
    
    contrasts1 <- tokenise.whitespace(gtxpipe.models[modelid, "contrasts"])
    resc <- lapply(contrasts1, function(contrast1) {
      message("Collating results for model ", gtxpipe.models[modelid, "model"], " for contrast ", contrast1)

      ## Are we better to do each contrast as a join, or lam the whole lot up and delete the unwanted columns later?
      
      groups <- unlist(strsplit(contrast1, "/"))
      stopifnot(identical(length(groups), 2L)) # internal error because checked this had length 2 previously
      group1 <- groups[1]
      group2 <- groups[2]
      
      ## More efficient to work with chi2 statistics, monotonic in P and using bespoke GC calculation
      ## ce = contrast expression
      ce <- parse(text = paste("list(chi2 = (beta.", group1, " - beta.", group2, ")^2/",
                    "(SE.GC.", group1, "^2 + SE.GC.", group2, "^2))", sep = ""))
      res1 <- res[[group1]][res[[group2]], eval(ce)]
      lambda <- res1[ , median(chi2, na.rm = TRUE)/qchisq(0.5, df = 1, lower.tail = FALSE)]
      setattr(res1, "lambda", lambda) # This is a Wald-like test but it's the primary one for the contrast
      invisible(if (lambda > 1.) {
        res1[ , pvalue.GC := pchisq(chi2/lambda, df = 1, lower.tail = FALSE)]
      } else {
        res1[ , pvalue.GC := pchisq(chi2, df = 1, lower.tail = FALSE)]
      })
      res1[ , chi2 := NULL]
      return(res1)
    })
    names(resc) <- contrasts1
    
    ## Note that
    resa <- c(res, resc)
                                        # is not fast
    
    ## Using the contrast names as part of column labels DOESN'T WORK
    return(do.call(rbind, lapply(names(resa), function(nn) {
      pcv <- resa[[nn]][cvlist , pvalue.GC]$pvalue.GC
      data.frame(model = gtxpipe.models[modelid, "model"],
                 group = nn,
                 lambda = attr(resa[[nn]], "lambda"),
                 thresh1 = thresh1,
                 hits1 = resa[[nn]][ , paste(sum(pvalue.GC <= thresh1, na.rm = TRUE), sum(!is.na(pvalue.GC)), sep = "/")],
                 ## hits1 !is.na() for denomin as constrasts have NA pvalues
                 hits2 = paste(sum(pcv <= thresh1, na.rm = TRUE), sum(!is.na(pcv)), sep = "/"),
                 thresh2 = thresh2,
                 stringsAsFactors = FALSE)
    })))
    ## Using primary threshold for non-primary analyses (specifically, groups used for contrasts but not of interest themselves
  }))
  
  ## Best place to compute per-analysis N???

  gtxpipe.analyses$index <- 1:nrow(gtxpipe.analyses)
  gtxpipe.results <- merge(gtxpipe.analyses, gtxpipe.results)
  gtxpipe.results <- gtxpipe.results[order(gtxpipe.results$index),]
  gtxpipe.results$index <- NULL


  message('gtxpipe: Writing top level summary results to "', file.path(getOption("gtxpipe.outputs"), "03_summary_results.csv"), '"')
  write.csv(gtxpipe.results,
            file = file.path(getOption("gtxpipe.outputs"), "03_summary_results.csv"),
            row.names = FALSE)
  
  ## contrasts, assume independent (user responsibility) but GC will roughly control for overlapping groups
  ## Note the GC coefficient is computed after groupwise GC
  

  ## (no need to read info files if no hits...)

  ## Default filtering on Rsq (at e.g. 0.05) would eliminate invariant dosages 0.3,...,0.3,...
  return(invisible(NULL))
}




## pipeslave, not exported but called from makefile as gtx:::pipeslave
pipeslave <- function(target) {
  target <- tryCatch(as.character(target)[1], error = function(e) "") # sanitize
  adir <- dirname(target) # analysis directory
  ofile <- file.path(adir, "options.R") # options file to read
  if (file.exists(ofile)) {
    source(ofile)
    ## should guarantee to set options gtx.usubjid, gtxpipe.genotypes,
    ## gtxpipe.threshold.MAF, and gtxpipe.threshold.Rsq
    ## but both have failsafe(?) fallbacks
  } else {
    warning('Options file "', ofile, '" does not exist')
  }
  usubjid <- getOption("gtx.usubjid", "USUBJID")
  tmp <- unlist(strsplit(basename(target), ".", fixed = TRUE))
  if (length(tmp) < 2 || !identical(tmp[length(tmp)], "done")) stop ('Bad target "', target, '"')
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
                    minimac = file.path(getOption("gtxpipe.genotypes", "genotypes"), job),
                    threshold.MAF = getOption("gtxpipe.threshold.MAF", 0.01),
                    threshold.Rsq = getOption("gtxpipe.threshold.Rsq", 0.01))
  ## Require pvalue always, beta and SE if being used for interaction contrasts
  ## Consider writing some comments to this file (blockassoc should return character vector instead of printing messages)
  write.table(res[ , c("SNP", "beta", "SE", "pvalue")], # should be an option to set outputs, esp analysed.*
              row.names = FALSE, quote = FALSE, 
              file = gzfile(file.path(adir, paste(job, "out.gz", sep = "."))))
  return(invisible(NULL))
}
        
