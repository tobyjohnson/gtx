## FIXME Using contrasts as a variable name shadows a base function !!

## option namespace should all be
## gtx.* or gtxpipe.*

## pgx.eigenvec is actually used as a list of genotyped subjects plus any covariates desired

## FIXME add direct hook for "user-derived" endpoints

gtxpipe <- function(gtxpipe.models = getOption("gtxpipe.models"),
                    gtxpipe.groups = getOption("gtxpipe.groups", data.frame(group = 'ITT', deps = 'pop.PNITT', fun = 'pop.PNITT', stringsAsFactors = FALSE)),
                    ## ugly to have this in the prototype (and hence verbatim in the man page)
                    gtxpipe.derivations = getOption("gtxpipe.derivations", {data(derivations.standard.IDSL); derivations.standard.IDSL}),
                    gtxpipe.transformations = getOption("gtxpipe.transformations", data.frame(NULL)),
                    gtxpipe.eigenvec,
                    stop.before.make = FALSE) {
  ## arguments for project-specific


  message("gtxpipe() from package gtx version ", packageVersion("gtx"), " on ", R.version.string)

  # R.version$os %in% c("linux-gnu", "cygwin")
  
  usubjid <- as.character(getOption("gtx.usubjid", "USUBJID"))[1] # variable name for unique subject identifier
  ## Replace this with a function getusubjid that prints warning messages, applies make.names() etc
  
  ## options where read repeatedly in code, constant (in C const sense)
  ## and likely constant over different projects on same IT
  ## overwrite with failsafe defaults now
  options(gtxpipe.genotypes = as.character(getOption("gtxpipe.genotypes", "genotypes"))[1]) # directory containing genotype data
  options(gtxpipe.clinical = as.character(getOption("gtxpipe.clinical", "clinical"))[1]) # directory containing clinical data
  options(gtxpipe.analyses = as.character(getOption("gtxpipe.analyses", "analyses"))[1]) # top level directory for analyses
  dir.create(getOption("gtxpipe.analyses"), recursive = TRUE, showWarnings = FALSE) # FIXME throw error if fails?
  options(gtxpipe.outputs = as.character(getOption("gtxpipe.outputs", "outputs"))[1]) # target directory for outputs
  dir.create(getOption("gtxpipe.outputs"), recursive = TRUE, showWarnings = FALSE) # FIXME throw error if fails?

  ## If user identity not specified, determine from USER environment variable
  options(gtxpipe.user = as.character(getOption("gtxpipe.user", paste("USER", Sys.getenv("USER", unset = "unknown"), sep = ":")))[1])

  
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
       if (any(bm)) stop("Models ", paste(gtxpipe.models[bm, "model"], collapse = ", "),
                         " have non-alphanumeric names.  You need to fix this in gtxpipe.models"))
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
    ## This would not be needed if default in the function arg default
    warning("No groups specified.  Defaulting to a single (ITT) group.")
    gtxpipe.groups <- data.frame(group = 'ITT', deps = 'pop.PNITT', fun = 'pop.PNITT', stringsAsFactors = FALSE)
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

  ## FIXME would be nice to use groupby to detect whether any contrasts use overlapping groups
  ## In theory, GC should fix this (albeit not optimally)
  
  ## The set of all variables we need for actual analysis
  ## Transformations -> new deps

  ## Compute model dependencies in terms of untransformed variables
  gtxpipe.models$depsu <- sapply(gtxpipe.models$deps, function(deps1) {
    paste(sapply(tokenise.whitespace(deps1), function(dep1) {
      mm <- match(dep1, gtxpipe.transformations$targets)
      if (is.na(mm)) return(dep1)
      return(gtxpipe.transformations$deps[mm])
    }), collapse = " ")})

  ## Should we do the same for group dependencies?
  
  deps <- unique(tokenise.whitespace(c(gtxpipe.models$depsu, gtxpipe.groups$deps,
                                       "pop.TRTGRP", "demo.SEX", "demo.AGE", "demo.RACE", "demo.ETHNIC")))
  ## force in pop.PNITT even though this is not a mandated variable per dsm ?
  ## allow a force in list.
  ## sort columns by unique(forcelist, deps) before computing demographics tables

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
  ## FIXME would be nicer if clinical.import had its own error checking; should it have "only" and "require" options??
  
  #message("gtxpipe: Read clinical datasets OK")
  anal1 <- clinical.derive(clindata, gtxpipe.derivations, only = deps)
  ## FIXME why is this so slow?
  message("gtxpipe: Computed derived variables OK")

  write.csv(anal1, 
            file = file.path(getOption("gtxpipe.outputs"), "subject_analysis_dataset.csv"))
            
  ## FIXME hook for user derived variables needed here (in case used in group defs etc)

  ## Even though we will re-apply transformations on subsets, they all should work on the complete analysis dataset
  if (nrow(gtxpipe.transformations) > 0) {
    for (idx in 1:nrow(gtxpipe.transformations)) {
      target <- gtxpipe.transformations$targets[idx]
      if (target %in% names(anal1)) stop("Transformation overwrites existing variable ", target)
      tryCatch(anal1[[target]] <- eval(parse(text = gtxpipe.transformations$fun[idx]), envir = anal1),
               error = function(e) 
               stop("Error in transformation '", target, " <- ", gtxpipe.transformations$fun[idx], "':\n", e$message))
    }
  }

  ## Order columns as follows:
  ##   USUBJID is first
  ##   if variable has a descriptor, include next, in the order appearing in clinical.descriptors
  ##   the rest
  anal1 <- anal1[ ,
                 order(ifelse(names(anal1) == usubjid, 0, 1),
                       match(names(anal1), names(getOption("clinical.descriptors", NULL))))]



  
#  inc1 <- with(anal1, eval(parse(text = getOption("clinical.subset", 'pop.PNITT'))))

  if (!missing(gtxpipe.eigenvec) && file.exists(gtxpipe.eigenvec)) {
    ancestrypcs <- read.table(gtxpipe.eigenvec, header = FALSE, sep = " ", as.is = TRUE)[ , -1] # drop first (FID) column
    names(ancestrypcs) <- c(usubjid, paste("PC", 1:(ncol(ancestrypcs) - 1), sep = ""))
  } else {
    warning("PGx eigenvectors not available.  Continuing imperfectly")
    ancestrypcs <- anal1[ , usubjid, drop = FALSE]
  }
  ## is.PGx indicates PGx only if gtxpipe.eigenvec is one-to-one list of PGx subjects
  anal1$is.PGx <- anal1[[usubjid]] %in% ancestrypcs[[usubjid]]

  ## gtxpipe needs to know number enrolled for the automagic report.
  ## Making an additional "ITT" group would have meant that number enrolled/PGx
  ## would automatically appear as a row in gtxpipe.groups, but then would
  ## have to program around unwanted columns appearing in the demographics table,
  ## and messy programming in the calculation of whether any groups overlap.
  ## Hence, storing this information separately.  Percent PGx can be calculated downstream
  ##
  ## Conceptually what we want is all subjects enrolled, which may not always
  ## be the same as the ITT population.  Could perhaps more robustly assume that
  ## enrolment corresponds to *some* flag pop.XXXX and use an option to specify?
  groupall <- with(anal1, list('All enrolled, ITT' = pop.PNITT,
                               'All enrolled, PGx' = pop.PNITT & is.PGx))

  ## FIXME tryCatch inside here
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
  ## NOTE we are assuming the user plays nice and if they define an overwriting transform this still works
  ## FIXME maybe we should not allow overwriting transforms?
  ## FIXME instead a clinical.ordering list of preferred orders // sort by count otherwise?
  gtxpipe.groups <- cbind(gtxpipe.groups, 
    do.call(rbind, lapply(groupby[paste(gtxpipe.groups$group, "ITT", sep = ", ")], function(gl) {
      tmp <- table(anal1$pop.TRTGRP[gl])
      tmp <- tmp[tmp > 0]
      if (length(tmp) == 0) return(data.frame(arms = "None", adjust.arm = NA, stringsAsFactors = FALSE))
      if (length(tmp) == 1) return(data.frame(arms = names(tmp), adjust.arm = FALSE, stringsAsFactors = FALSE))
      return(data.frame(arms = paste(names(tmp), " (N=", tmp, ")", sep = "", collapse = ", "),
                        adjust.arm = TRUE, stringsAsFactors = FALSE))
    })))

  snippets <- rbind(data.frame(value = format(sapply(groupall, sum)), # automatic row names
                               stringsAsFactors = FALSE),
                    data.frame(value = c(round(100*sum(groupall[[2]])/sum(groupall[[1]])),
                                 if (any(rowSums(as.data.frame(groupby)) > 1)) "overlapping" else "non-overlapping",
                                 if (any(gtxpipe.groups$adjust.arm, na.rm = TRUE)) "Yes" else "No"),
                               ## rowSums all 0 or 1, no overlapping groups.  Any >1 implies overlapping groups
                               ## Using Yes/No because output all text so downstream would have to parse even if logical
                               row.names = c("Overall PGx percent", "PGx group overlap", "PGx combines arms")))
  ## Should add the study name, groups, models, whether "efficacy" or "efficacy and safety" etc
 
  gtxpipe.groups$N.notPGx <-   gtxpipe.groups$N.ITT - gtxpipe.groups$N.PGx

  ## source item 2, can we call this disposition?
  ## write into source subdir
  rownames(gtxpipe.groups) <- gtxpipe.groups$group
#  print(gtxpipe.groups[ , c("group", "arms", "N.ITT", "N.PGx")])

  ## first call, no metadata
  metadata <- pipetable(gtxpipe.groups[ , c("group", "arms", "N.ITT", "N.PGx")],
                        "02", "subject_disposition",
                        "Disposition of subjects by PGx treatment group and availability of PGx data")
#  message('gtxpipe: Writing subject disposition to "', file.path(getOption("gtxpipe.outputs"), "02_subject_disposition.csv"), '"')
#  write.csv(gtxpipe.groups[ , c("group", "arms", "N.ITT", "N.PGx")],
#            file = file.path(getOption("gtxpipe.outputs"), "02_subject_disposition.csv"),
#            row.names = TRUE)
  
             
  # apply(gtxpipe.groups[ , c("N.notPGx", "N.PGx")], 1, prettypc)[2, , drop = TRUE]




  
  ## allow a demographics sort list
  ## write into source subdir
  metadata <- pipetable(demographics(anal1, by = groupby), 
                        "03", "subject_demographics",
                        "Demographics of subjects included in PGx analyses",
                        metadata)
#  message('gtxpipe: Writing subject demographics to "', file.path(getOption("gtxpipe.outputs"), "03_subject_demographics.csv"), '"')
#  write.csv(demographics(anal1, by = groupby),
#            file = file.path(getOption("gtxpipe.outputs"), "03_subject_demographics.csv"),
#            row.names = TRUE)

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

  ## Currently we compute analysis N's below.  We need to do that earlier
  ## because trying to analyse N=0 groups causes lots of problems.
  ## Need to solve circular argument that adjust.arm depends on whether multiple arms
  ## in the analysis dataset and size of analysis dataset depends on number of subjects
  ## with non-missing covariates (which may include arm).  (And note we may have studies with
  ## missing TRTGRP or ATRTGRP)
  
  adir0 <- getOption("gtxpipe.analyses")

  sink("Makefile") # name of this file should be an option

  cat("## Makefile generated by gtxpipe()\n")
  cat(".PHONY:\tall\nall:\t\n\n") # check having all twice is okay
  
  cat("GENOTYPES := $(wildcard ", getOption("gtxpipe.genotypes"), "/*.dose.gz)\n\n", sep = "")

  ## Loop over analyses not nested loops, set N in gtxpipe.analyses
  ## Write analysis datasets and Makefile as side effects
  analN <- do.call(rbind, lapply(1:nrow(gtxpipe.models), function(modelid) {
    return(do.call(rbind, lapply(tokenise.whitespace(gtxpipe.models[modelid, "agroups"]), function(agroup1) {
      groupid <- match(agroup1, gtxpipe.groups$group)
      ## FIXME next line throws an error if adjust.arm is NA (implies N=0 in ITT)
      trtgrp.cov <- if (gtxpipe.groups$adjust.arm[groupid]) "pop.TRTGRP" else NULL
      adir <- file.path(adir0, gtxpipe.models[modelid, "model"], gtxpipe.groups[groupid, "group"])
      dir.create(adir, recursive = TRUE, showWarnings = FALSE)
      odir = file.path(getOption("gtxpipe.outputs"), gtxpipe.models[modelid,"model"], gtxpipe.groups[groupid, "group"])
      dir.create(odir, recursive = TRUE, showWarnings = FALSE)
      ## Ensure that relevant options set in this gtxpipe() call are available
      ## to slave calls
      sink(file.path(adir, "options.R")) # sink inside sink
      cat('options(gtx.usubjid = "', usubjid, '")\n', sep = '')
      cat('options(gtxpipe.genotypes = "', getOption("gtxpipe.genotypes"), '")\n', sep = '') # will break if not scalar!
      cat('options(gtxpipe.threshold.MAF = ', getOption("gtxpipe.threshold.MAF", 0.01), ')\n', sep = '')
      cat('options(gtxpipe.threshold.Rsq = ', getOption("gtxpipe.threshold.Rsq", 0.01), ')\n', sep = '')
      ## FIXME need to pass through the names of any packages needed to run blockassoc inside pipeslave
      sink() # options.R
      


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
      for (idx in which(gtxpipe.transformations$targets %in% tokenise.whitespace(gtxpipe.models[modelid, "deps"]))) {
        target <- gtxpipe.transformations$targets[idx]
        if (target %in% names(adata)) stop("Transformation overwrites existing variable")
        message(target, " <- ", gtxpipe.transformations$fun[idx])
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
      cat('\t@(echo "library(gtx); gtx:::pipeslave(target = \\"$@\\")" | R --quiet --vanilla && touch $@)\n\n', sep = '')
      ## delete output, delete done file, run R, touch done file
      ## 		mkdir -p a1; rm -f $@; sleep 60; uname -a >$@; date >>$@
      ## Once all chunks "done", combine into a single tabix'd file
      ## Per qmake man page, need to combine multiple commands into single rule by separating with ;
      cat(adir, '/ALL.out.txt.gz: $(MODEL', modelid, 'GROUP', groupid, ')\n', sep='')
      cat("\tzgrep -h '^SNP' ", adir, "/*out.gz | head -n 1 | awk 'BEGIN{FS=", '"\\t";OFS="\\t";} {print "#CHROM","POS",$$0;}', "' >", 
          adir, '/ALL.out.txt ; \\\n', sep='')
      cat("\tzgrep -h -v '^SNP' ", adir, "/*out.gz | awk 'BEGIN{FS=", '"\\t";OFS="\\t";} {split($$1,coord,"[:_]"); print coord[1],coord[2],$$0;}', 
          "' | sort -T . -k 1,1 -k 2,2n >>", adir, '/ALL.out.txt ; \\\n', sep='')
      cat('\tbgzip -f ', adir, '/ALL.out.txt\n\n', sep='')
      cat(adir, '/ALL.out.txt.gz.tbi: ', adir, '/ALL.out.txt.gz\n', sep='')
      cat('\ttabix -f -b 2 -e 2 ', adir, '/ALL.out.txt.gz ; \\\n', sep='')
      cat('\trm ', adir, '/*out.gz\n\n', sep='')
      ## Run python script to generate genome-wide plots
      cat(odir, '/qq.png: ', adir, '/ALL.out.txt.gz ', adir, '/ALL.out.txt.gz.tbi\n', sep='')
      cat('\t', find.package('gtx'), '/python/plots.py -r ', adir, '/ALL.out.txt.gz -o ', odir, '\n\n', sep='')
      
      return(data.frame(model = gtxpipe.models[modelid, "model"], group = agroup1, N = nrow(adata),
                        makevar = paste(odir, '/qq.png', sep = ''),
                        stringsAsFactors = FALSE))
    })))
  }))

  cat("all:\t", paste(analN$makevar, collapse = " "), "\n", sep = "") 
  sink() # Makefile

  gtxpipe.analyses$index <- 1:nrow(gtxpipe.analyses)
  gtxpipe.analyses <- merge(gtxpipe.analyses, analN, all.x = TRUE, all.y = FALSE)

  ## FIXME Would be nice to automagically compute N or N1/N2 for contrasts.  Needs to be done within levels of model.

  metadata <- pipetable(snippets,
                        "01", "study_summary", "PGx study summary",
                        metadata)
#  write.csv(snippets, 
#            file = file.path(getOption("gtxpipe.outputs"), "01_study_summary.csv"),
#            row.names = TRUE)
  
  if (stop.before.make) return(invisible(NULL))
  
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

      ## Reading results which were compiled across chunks during make call
      res1 <- read.table(gzfile(file.path(adir, "ALL.out.txt.gz")),
                          quote = "", comment.char = "", header = TRUE, stringsAsFactors = FALSE)
      ## Confirm expected columns present
      stopifnot(all(c("SNP", "pvalue", "beta", "SE") %in% names(res1)))
      ## Convert to data table to improve efficiency retaining only needed columns and rows with non-missing pvalue
      res1 <- data.table(res1[!is.na(res1$pvalue), c("SNP", "pvalue", "beta", "SE"), drop = FALSE])
      ## Sort by SNP
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
    
    ## Note that ...
    resa <- c(res, resc) # ... is not fast
    
    ## Using the contrast names as part of column labels DOESN'T WORK
    return(do.call(rbind, lapply(names(resa), function(nn) {
      pcv <- resa[[nn]][cvlist , pvalue.GC]$pvalue.GC
      data.frame(model = gtxpipe.models[modelid, "model"],
                 group = nn,
                 lambda = attr(resa[[nn]], "lambda"),
                 thresh1 = thresh1,
                 hits1 = resa[[nn]][ , paste(sum(pvalue.GC <= thresh1, na.rm = TRUE), sum(!is.na(pvalue.GC)), sep = "/")],
                 ## hits1 !is.na() for denomin as constrasts have NA pvalues
                 hits2 = paste(sum(pcv <= thresh2, na.rm = TRUE), sum(!is.na(pcv)), sep = "/"),
                 thresh2 = thresh2,
                 stringsAsFactors = FALSE)
    })))
    ## Using primary threshold for non-primary analyses (specifically, groups used for contrasts but not of interest themselves
  }))
  
  ## Best place to compute per-analysis N???

  gtxpipe.results <- merge(gtxpipe.analyses, gtxpipe.results)
  gtxpipe.results <- gtxpipe.results[order(gtxpipe.results$index),]
  rownames(gtxpipe.results) <- 1:nrow(gtxpipe.results)
  gtxpipe.results$index <- NULL
  gtxpipe.results$makevar <- NULL
  
  ## set row names to something meaningful?
  metadata <- pipetable(gtxpipe.results,
                        "04", "summary_results", "PGx analysis summary",
                        metadata)

  message('gtxpipe: Writing source display metadata')
  write.csv(metadata,
            file.path(getOption("gtxpipe.outputs"), "lela_metadata"))
  
  ## contrasts, assume independent (user responsibility) but GC will roughly control for overlapping groups
  ## Note the GC coefficient is computed after groupwise GC

  message('gtxpipe: Generating report')
  ### 
  ### Make "short" report, should switch on presence/absence of positive results
  ###
  if (file.exists(file.path(getOption("gtxpipe.outputs"), "report-short.Rmd"))) {
    # do nothing
  } else {
    file.copy(system.file("templates/gtxpipe-report-negative.Rmd",
                          package = "gtx", mustWork = TRUE),
              file.path(getOption("gtxpipe.outputs"), "report-short.Rmd"))
  }
  tryCatch({
    suppressMessages(requireNamespace("knitr")) || stop("knitr package not available")
    oldwd <- getwd()
    setwd(getOption("gtxpipe.outputs"))
    knitr::knit2html("report-short.Rmd")
    setwd(oldwd)
  },
           error = function(e) {
             message("knitr::knit2html failed")
           })
  
  ## Note, no need to read info files if no hits...

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
              row.names = FALSE, quote = FALSE, sep = "\t",
              file = gzfile(file.path(adir, paste(job, "out.gz", sep = "."))))
  return(invisible(NULL))
}
        
pipetable <- function(data, number, filename, title, mdata) {
  if (missing(mdata)) mdata <- data.frame(NULL)

  path <- file.path(getOption("gtxpipe.outputs"), paste(number, filename, sep = "_"))
    message('gtxpipe: Writing ', title, ' to "', path, '.[csv|pdf]"')

  write.csv(data,
            file = paste(path, "csv", sep = "."),
            row.names = TRUE) # ???

  pdf(width = 8.3, height = 11.7,
      file = paste(path, "pdf", sep = "."))
  par(family="mono", cex.main = 1, cex.sub = 1)
  plot.new()
  plot.window(c(0, 1), c(0, 1))
  textgrid(data)
  title(main = title, sub = paste("Source table", number))
  dev.off()

  return(rbind(mdata,
               data.frame(Source_File_Name = paste(paste(number, filename, sep = "_"), "pdf", sep = "."), 
                          Display_Category = "PHARMACOGENETIC", 
                          Display_Type = "TABLE", 
                          Display_Number = number,
                          Title = title,
                          stringsAsFactors = FALSE)))
}
