## FIXME Using contrasts as a variable name shadows a base function !!

## option namespace should all be
## gtx.* or gtxpipe.*

## pgx.eigenvec is actually used as a list of genotyped subjects plus any covariates desired

## FIXME add direct hook for "user-derived" endpoints

gtxpipe <- function(gtxpipe.models = getOption("gtxpipe.models"),
                    gtxpipe.groups = getOption("gtxpipe.groups", data.frame(group = 'ITT', deps = 'pop.PNITT', fun = 'pop.PNITT', stringsAsFactors = FALSE)),
                    ## ugly to have this in the prototype (and hence verbatim in the man page)
                    ## Also R CMD CHECK thinks derivations.standard.IDSL is an unbound global variable
                    gtxpipe.derivations = getOption("gtxpipe.derivations", {data(derivations.standard.IDSL); derivations.standard.IDSL}),
                    gtxpipe.transformations = getOption("gtxpipe.transformations", data.frame(NULL)),
                    gtxpipe.eigenvec,
                    stop.before.make = FALSE) {
  ## arguments for project-specific

  data(gtx.version)
  message("gtxpipe() from package gtx version ", as.character(packageVersion("gtx")), " on ", R.version.string)
  message("gtx package build: ", gtx.version[1])
  # R.version$os %in% c("linux-gnu", "cygwin")

  ## Load any packages specified by gtxpipe.packages option, clean up ready to pass to pipeslave
  packages <- intersect(getOption('gtxpipe.packages', NULL), rownames(installed.packages()))
  for (package in packages) {
    message("gtxpipe: loading package '", package, "'")
    if (!require(package, character.only = TRUE, quietly = TRUE)) {
      stop("gtxpipe: fatal error, unable to load package '", package, "'")
    }
  }
  ## Note, any packages not available are silently ignored
  options(gtxpipe.packages = packages)
  
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

  ## If date not specified, determine from system date
  options(gtxpipe.date = as.character(getOption("gtxpipe.date", format(Sys.Date(), "%Y-%b-%d")))[1])
  
  ## Check genotypes directory contains at least one dose/info pair and enumerate for check against done files later
  doses = gsub('\\.dose.gz$','',list.files(path = getOption("gtxpipe.genotypes"), pattern = '\\.dose.gz$'))
  infos = gsub('\\.info.gz$','',list.files(path = getOption("gtxpipe.genotypes"), pattern = '\\.info.gz$'))
  chunks = intersect(doses,infos)
  if (length(chunks) < 1) stop("Needs to be at least one *.dose.gz/*.info.gz file pair in the genotypes directory [", 
                                getOption("gtxpipe.genotypes"), "] you can set this directory as the gtxpipe.genotypes option")
  
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
       if (any(bg)) stop("Groups ", paste(groups[bg, "group"], collapse = ", "),
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
                                       "pop.PNITT", "pop.TRTGRP", "demo.SEX", "demo.AGE", "demo.RACE", "demo.ETHNIC")))
  ## force in pop.PNITT even though this is not a mandated variable per dsm ?
  ## Code below for groupall assumes pop.PNITT exists - there will be an error if not FIXME
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

  ## FIXME ancestry PCs are not written here, note if merge, should merge all=TRUE since
  ## we will want clinical data for non-PGx subjects
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
  ## The definition used here should be linked to the "Population: ITT" header in displays
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

  snippets <- rbind(data.frame(value = c(
                                 getOption("gtxpipe.protocol", "NA"),
                                 getOption("gtxpipe.project", "NA"),
                                 getOption("gtxpipe.user", "NA"),
                                 getOption("gtxpipe.email", "NA"),
                                 getOption("gtxpipe.date", "NA"),
                                 R.version.string,
                                 as.character(packageVersion("gtx")),
                                 gtx.version[1]),
                               row.names = c("Protocol", "Project", "User", "Email", "DataAsOf", "R.version", "gtx.package.version", "gtx.package.build"),
                               stringsAsFactors = FALSE),
                    data.frame(value = format(sapply(groupall, sum)), # automatic row names
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
                        "subject_disposition",
                        "Subject disposition by PGx analysis group",
                        number = 2) # specifying number because first 4 tables are generated out-of-order
#  message('gtxpipe: Writing subject disposition to "', file.path(getOption("gtxpipe.outputs"), "02_subject_disposition.csv"), '"')
#  write.csv(gtxpipe.groups[ , c("group", "arms", "N.ITT", "N.PGx")],
#            file = file.path(getOption("gtxpipe.outputs"), "02_subject_disposition.csv"),
#            row.names = TRUE)
  
             
  # apply(gtxpipe.groups[ , c("N.notPGx", "N.PGx")], 1, prettypc)[2, , drop = TRUE]




  
  ## allow a demographics sort list
  ## write into source subdir
  metadata <- pipetable(demographics(anal1, by = groupby), 
                        "subject_demographics",
                        "Demographics of subjects included in PGx analyses",
                        metadata, number = 3) # specifying number because first 4 tables are generated out-of-order
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

  thisR <- getOption("gtxpipe.slaveR", file.path(R.home(component = "bin"), "R"))
  
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

      ## Ensure that relevant options set in this gtxpipe() call are available
      ## to slave calls
      sink(file.path(adir, "options.R")) # sink inside sink
      cat('options(gtx.usubjid = "', usubjid, '")\n', sep = '')
      cat('options(gtxpipe.genotypes = "', getOption("gtxpipe.genotypes"), '")\n', sep = '') # will break if not scalar!
      cat('options(gtxpipe.threshold.MAF = ', getOption("gtxpipe.threshold.MAF", 0.01), ')\n', sep = '')
      cat('options(gtxpipe.threshold.Rsq = ', getOption("gtxpipe.threshold.Rsq", 0.01), ')\n', sep = '')
      cat('options(gtxpipe.extended.output = ', getOption("gtxpipe.extended.output", TRUE), ')\n', sep = '')
      if (length(getOption("gtxpipe.packages")) > 0) {
        cat('options(gtxpipe.packages = c(',
            paste('\'', getOption("gtxpipe.packages"), '\'', sep = '', collapse = ', '),
            '))\n', sep = '')
      }
      sink() # options.R

      adata <- merge(anal1[groupby[[paste(agroup1, ", PGx", sep = "")]],
                           c(usubjid,  
                             tokenise.whitespace(gtxpipe.models[modelid, "depsu"]),
                             trtgrp.cov)],
                     ancestrypcs,
                     all = FALSE)
      ## Remove subjects with any missing data on any dependency variable
      adata <- adata[which(apply(!is.na(adata), 1, all)), ]
      ## Transformations must be applied to *the analysis dataset*
      ## ONLY APPLY RELEVANT TRANSFORMATIONS depsu -> deps
      ## Although transformations are applied by gtxpipeslave
      ##   (so that transformations can create R classes such as factor, Surv etc)
      ##   transformations also applied here to remove subjects missing-after-transformation
      ##   and to fit null model
      for (idx in which(gtxpipe.transformations$targets %in% tokenise.whitespace(gtxpipe.models[modelid, "deps"]))) {
        target <- gtxpipe.transformations$targets[idx]
        if (target %in% names(adata)) stop("Transformation overwrites existing variable")
        message(target, " <- ", gtxpipe.transformations$fun[idx])
        adata[[target]] <- eval(parse(text = gtxpipe.transformations$fun[idx]), envir = adata)
      }
      ## Remove any subjects with missing model deps
      ##   (note, missing depsu might be okay e.g. derivation returns missing which is converted to
      ##    numeric by transformation e.g. Surv2(t1, t2))
      sel <- apply(!is.na(adata[ , unique(c(usubjid, 
                                            tokenise.whitespace(gtxpipe.models[modelid, "deps"]), 
                                            trtgrp.cov,
                                            names(ancestrypcs)))]), 1, all)
      adata <- adata[sel, , drop = FALSE]

      ## Everything except usubjid will be a factor

      ## Everything as paths relative to getwd() at runtime

      message("gtxpipe: Fitting model ", gtxpipe.models[modelid, "model"], " in group ", agroup1)
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
      ## first 8 characters MUST be '# call : ' followed by call to fit null model
      cat('# call: ', as.character(as.expression(m0$call)), '\n', sep = '')

      ## optional lines for transformations, e.g. '# where: osMonths <- Surv(SRVMO, SRVCFLCD)'
      for (idx in which(gtxpipe.transformations$targets %in% tokenise.whitespace(gtxpipe.models[modelid, "deps"]))) {
        cat('# where: ', gtxpipe.transformations$targets[idx], ' <- ', gtxpipe.transformations$fun[idx], '\n', sep = '')
      }
      
      ## including candidate variant list to force analysis regardless of MAF and RSQR filters
      ## Note that ideally this should be passed through the options.R file, not in the model/data spec
      cat('# cvlist: ', gtxpipe.models[modelid, "cvlist"], '\n', sep = '')

      sink()
      ## Only write untransformed data since pipeslave will reapply transformations
      suppressWarnings(write.table(adata[ , unique(c(usubjid, 
                                                     tokenise.whitespace(gtxpipe.models[modelid, "depsu"]), 
                                                     trtgrp.cov,
                                                     names(ancestrypcs)))], 
                                   sep = ",", row.names = FALSE, 
                                   file = file.path(adir, "analysis-dataset.csv"),
                                   append = TRUE))

      #If cvlist specified, convert to bed file of regions +/- 500 kb for tabix extract of full genome results
      if (!is.na(gtxpipe.models[modelid, "cvlist"]) && gtxpipe.models[modelid, "cvlist"] != '') {
        pipebed(snplist = tokenise.whitespace(gtxpipe.models[modelid, "cvlist"]),
                flank = 500000, outfile = file.path(adir, "CV.bed"))
      }

      cat('# Analysis for model "', gtxpipe.models[modelid, "model"], '" in group "', gtxpipe.groups[groupid, "group"], '"\n', sep = '')
      cat('MODEL', modelid, 'GROUP', groupid, ' := $(patsubst ', getOption("gtxpipe.genotypes"), '/%.dose.gz,', adir, '/%.done,$(GENOTYPES))\n', sep = '')
      cat('$(MODEL', modelid, 'GROUP', groupid, '):\t', adir, '/%.done:\t', getOption("gtxpipe.genotypes"), '/%.info.gz ', getOption("gtxpipe.genotypes"), '/%.dose.gz\n', sep = '')
      cat('\t@(echo "library(gtx); gtx:::pipeslave(target = \\"$@\\")" | ', thisR, ' --quiet --vanilla && uname -a >$@)\n\n', sep = '')
      ## delete output, delete done file, run R, touch done file
      ## 		mkdir -p a1; rm -f $@; sleep 60; uname -a >$@; date >>$@
      ## Once all chunks "done", combine into a single tabix'd file
      ## Per qmake man page, need to combine multiple commands into single rule by separating with ;
      cat(adir, '/ALL.out.txt.gz: $(MODEL', modelid, 'GROUP', groupid, ')\n', sep='')
      cat("\tzgrep -h '^SNP' ", adir, "/*out.gz | head -n 1 | awk 'BEGIN{FS=", '"\\t";OFS="\\t";} {print "#CHROM","POS",$$0;}', "' | gzip >", 
          adir, '/ALL.out.txt.gz0 ; \\\n', sep='')
      cat('\tif [ -e "', adir, '/ALL.out.txt.gz" ]; then zgrep -v ', "'^#CHROM' ", adir, '/ALL.out.txt.gz | gzip >>', adir, '/ALL.out.txt.gz0 ; fi ; \\\n', sep='')
      cat("\tzgrep -h -v '^SNP' ", adir, "/*out.gz | awk 'BEGIN{FS=", '"\\t";OFS="\\t";} {split($$1,coord,"[:_]"); print coord[1],coord[2],$$0;}', 
          "' | gzip >>", adir, '/ALL.out.txt.gz0 ; \\\n', sep='')
      cat('\tzcat ', adir, '/ALL.out.txt.gz0 | sort -T . -k 1,1 -k 2,2n | uniq | bgzip -f >', adir, '/ALL.out.txt.gz ; \\\n', sep='')
      cat('\trm ', adir, '/ALL.out.txt.gz0\n\n', sep='')
      cat(adir, '/ALL.out.txt.gz.tbi: ', adir, '/ALL.out.txt.gz\n', sep='')
      cat('\ttabix -f -b 2 -e 2 ', adir, '/ALL.out.txt.gz ; \\\n', sep='')
      cat('\trm ', adir, '/*out.gz\n\n', sep='')
      #If cvlist defined, extract results
      if (!is.na(gtxpipe.models[modelid, "cvlist"]) && gtxpipe.models[modelid, "cvlist"] != '') {
        cat(adir, '/CV.out.txt.gz: ', adir, '/ALL.out.txt.gz.tbi\n',sep='')
        cat('\ttabix -hB ', adir, '/ALL.out.txt.gz ', adir, '/CV.bed | gzip > ', adir, '/CV.out.txt.gz\n\n', sep='')
        return(data.frame(model = gtxpipe.models[modelid, "model"], group = agroup1, N = nrow(adata),
                        makevar = paste(adir, '/CV.out.txt.gz', sep = ''),
                        stringsAsFactors = FALSE))
      }
      else {
        return(data.frame(model = gtxpipe.models[modelid, "model"], group = agroup1, N = nrow(adata),
                        makevar = paste(adir, '/ALL.out.txt.gz.tbi', sep = ''),
                        stringsAsFactors = FALSE))
      }
    })))
  }))

  cat("all:\t", paste(analN$makevar, collapse = " "), "\n", sep = "") 
  sink() # Makefile

  gtxpipe.analyses$index <- 1:nrow(gtxpipe.analyses)
  gtxpipe.analyses <- merge(gtxpipe.analyses, analN, all.x = TRUE, all.y = FALSE)

  ## FIXME Would be nice to automagically compute N or N1/N2 for contrasts.  Needs to be done within levels of model.

  metadata <- pipetable(snippets,
                        "study_summary", "PGx study summary",
                        metadata,
                        number = 1) # specifying number because first 4 tables are generated out-of-order
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
  makesuccess <- system(paste(getOption("gtxpipe.make", "make"), " 1>Makefile.out 2>Makefile.err", sep = ""))
  if (makesuccess != 0) {
    stop("make failed, check Makefile.out and Makefile.err")
  }


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
      message("gtxpipe: Collating results for model ", gtxpipe.models[modelid, "model"], " in group ", agroup1)
      adir <- file.path(adir0, gtxpipe.models[modelid, "model"], agroup1)

      ## make can silently fail, need to check for .done files
      dones = gsub('\\.done$','',list.files(path = adir,pattern = '\\.done$'))
      if (length(dones) < length(chunks)) stop("Not all chunks analysed for model [", gtxpipe.models[modelid, "model"], 
                                               "] in group [", agroup1, "]. Try re-running to capture the missing chunks [", 
                                               paste(setdiff(chunks,dones), collapse=", "), "]")

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

      plotdata <- rbind(snippets["Project", , drop = FALSE],
                        data.frame(value = c(lambda > 1., round(lambda, 4), 
                                     gtxpipe.models[modelid, "model"],
                                     agroup1,
                                     res1[ , sum(!is.na(res1$pvalue.GC))]), 
                                   row.names = c("GenomicControl", "Lambda",
                                     "Model",
                                     "Subgroup",
                                     "PValues"),
                                   stringsAsFactors = FALSE))
      
      ## Note, QQ and Manhattan plots are drawn *after* genomic control
      assign("metadata", pipeplot('res1[ , qq10(pvalue.GC, pch = 20)]',
                                  filename = paste("QQ", gtxpipe.models[modelid,"model"], agroup1, sep = "_"),
                                  title = paste("QQ plot for", gtxpipe.models[modelid,"model"], "in group", agroup1),
                                  metadata,
                                  number = 5, # *start* at 5 to leave space for 04_summary_results
                                  plotdata = plotdata,
                                  plotpar = list(mar = c(4, 4, 0, 0) + 0.1)),
             pos = parent.frame(n = 4))
      ## Have to use assign(..., pos = ) to update metadata from inside two levels of nested anonymous function
      assign("metadata", pipeplot('res1[ , manhattan(pvalue.GC, SNP, pch = 20, cex = 0.5)]', 
                                  filename = paste("Manhattan", gtxpipe.models[modelid,"model"], agroup1, sep = "_"),
                                  title = paste("Manhattan plot for", gtxpipe.models[modelid,"model"], "in group", agroup1),
                                  metadata,
                                  number = 5, # *start* at 5 to leave space for 04_summary_results
                                  plotdata = plotdata,
                                  plotpar = list(mar = c(4, 4, 0, 0) + 0.1)),
             pos = parent.frame(n = 4))
      ## Have to use assign(..., pos = ) to update metadata from inside two levels of nested anonymous function

      #Create bed file to tabix out genome-wide significant results
      if (sum(res1$pvalue.GC <= thresh1, na.rm = TRUE) > 0) {
        pipebed(snplist = res1[res1$pvalue.GC <= thresh1,SNP], flank = 500000, 
                outfile = file.path(adir, "GenomeWideSignif.bed"))
        #Run tabix extraction
        message("gtxpipe: Extracting genome-wide significant results now")
        tabixsuccess <- system(paste("tabix -hB ", adir , "/ALL.out.txt.gz ",
                                  adir , "/GenomeWideSignif.bed | gzip >", adir, 
                                  "/GenomeWideSignif.out.txt.gz", sep = ""))
        if (tabixsuccess != 0) {
          warning("gtxpipe: tabix extraction of genome-wide significant results failed")
        }
      }

      ## Apply GC to SEs using pvalues from Wald test
      lambdaWald <- res1[ , gclambda(pchisq((beta/SE)^2, df = 1, lower.tail = FALSE))]
      setattr(res1, "lambdaWald", lambdaWald)
      invisible(if (lambdaWald > 1.) {
        res1[ , SE := SE*sqrt(lambdaWald)] # overwrite, save memory
      })
      setnames(res1, "SE", "SE.GC")
      ## Set SE.GC to NA when mis-calibrated,
      ## working definition being Wald chisquare statistic > 2. times the LRT chisquare statistic
      res1[(beta/SE.GC)^2 / qchisq(pvalue.GC, df = 1, lower.tail = FALSE) > 2., SE.GC := NA]

      setnames(res1, "SE.GC", paste("SE.GC", agroup1, sep = "."))
      setnames(res1, "beta", paste("beta", agroup1, sep = ".")) # named by groups, used in clever join in constrasts
      return(res1)
    })
    names(res) <- agroups
    
    contrasts1 <- tokenise.whitespace(gtxpipe.models[modelid, "contrasts"])
    resc <- lapply(contrasts1, function(contrast1) {
      message("gtxpipe: Collating results for model ", gtxpipe.models[modelid, "model"], " for contrast ", contrast1)

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

      ## Note we can't use '/' separated contrasts in filenames
      
      plotdata <- rbind(snippets["Project", , drop = FALSE],
                        data.frame(value = c(lambda > 1., round(lambda, 4), 
                                     gtxpipe.models[modelid, "model"],
                                     paste(group1, group2, sep = " vs "), 
                                     res1[ , sum(!is.na(res1$pvalue.GC))]), 
                                   row.names = c("GenomicControl", "Lambda",
                                     "Model",
                                     "Contrast",
                                     "PValues"),
                                   stringsAsFactors = FALSE))
      
      assign("metadata", pipeplot('res1[ , qq10(pvalue.GC, pch = 20)]',
                                  filename = paste("QQ", gtxpipe.models[modelid, "model"], group1, "vs", group2, sep = "_"),
                                  title = paste("QQ plot for", gtxpipe.models[modelid, "model"], "contrasting", group1, "vs", group2),
                                  metadata,
                                  number = 5, # *start* at 5 to leave space for 04_summary_results
                                  plotdata = plotdata,
                                  plotpar = list(mar = c(4, 4, 0, 0) + 0.1)),
             pos = parent.frame(n = 4))
      ## Have to use assign(..., pos = ) to update metadata from inside two levels of nested anonymous function
      assign("metadata", pipeplot('res1[ , manhattan(pvalue.GC, SNP, pch = 20, cex = 0.5)]',
                                  ## Note manhattan() must cope with missing pvalue.GC from SNPs in group1 but not in group2
                                  filename = paste("Manhattan", gtxpipe.models[modelid,"model"], group1, "vs", group2, sep = "_"),
                                  title = paste("Manhattan plot for", gtxpipe.models[modelid,"model"], "contrasting", group1, "vs", group2),
                                  metadata,
                                  number = 5, # *start* at 5 to leave space for 04_summary_results
                                  plotdata = plotdata,
                                  plotpar = list(mar = c(4, 4, 0, 0) + 0.1)),
             pos = parent.frame(n = 4))
      ## Have to use assign(..., pos = ) to update metadata from inside two levels of nested anonymous function

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
  gtxpipe.results$lambda <- round(gtxpipe.results$lambda, 4)
  
  ## set row names to something meaningful?
  metadata <- pipetable(gtxpipe.results,
                        "summary_results", "PGx analysis summary",
                        metadata,
                        number = 4) # specifying number because first 4 tables are generated out-of-order

  if (all(c("PC1", "PC2", "PC3") %in% names(ancestrypcs))) {
    ## Hard coding the options here; we want RACE; ETHNIC if available otherwise RACE,
    if (all(c("demo.RACE", "demo.ETHNIC") %in% names(anal1))) {
      adata <- merge(anal1, ancestrypcs, all = FALSE)
      metadata <- pipeplot('with(adata, pcaplot(cbind(PC1, PC2, PC3), paste(as.character(demo.RACE), as.character(demo.ETHNIC), sep = "; ")))', 
                           filename = "Ancestry_PC_RACE_ETHNIC", 
                           title = "Ancestry principal components by race and ethnicity",
                           metadata,
                           plotdata = rbind(snippets["Project", , drop = FALSE],
                             data.frame(value = "demo.RACE; demo.ETHNIC", 
                                        row.names = c("GroupedBy"),
                                        stringsAsFactors = FALSE)))
      ## plotpar not needed since pcaplot calls par within each (sub)screen
    } else if ("demo.RACE" %in% names(anal1)) {
      adata <- merge(anal1, ancestrypcs, all = FALSE)
      metadata <- pipeplot('with(adata, pcaplot(cbind(PC1, PC2, PC3), as.character(demo.RACE)))', 
                           filename = "Ancestry_PC_RACE", 
                           title = "Ancestry principal components by race",
                           metadata,
                           plotdata = rbind(snippets["Project", , drop = FALSE],
                             data.frame(value = "demo.RACE", 
                                        row.names = c("GroupedBy"),
                                        stringsAsFactors = FALSE)))
    }
    if ("demo.COUNTRY" %in% names(anal1)) {
      adata <- merge(anal1, ancestrypcs, all = FALSE)
      metadata <- pipeplot('with(adata, pcaplot(cbind(PC1, PC2, PC3), as.character(demo.COUNTRY)))', 
                           filename = "Ancestry_PC_COUNTRY", 
                           title = "Ancestry principal components by country",
                           metadata,
                           plotdata = rbind(snippets["Project", , drop = FALSE],
                             data.frame(value = "demo.COUNTRY", 
                                        row.names = c("GroupedBy"),
                                        stringsAsFactors = FALSE)))
    }
    if ("demo.REGION" %in% names(anal1)) {
      adata <- merge(anal1, ancestrypcs, all = FALSE)
      metadata <- pipeplot('with(adata, pcaplot(cbind(PC1, PC2, PC3), as.character(demo.REGION)))', 
                           filename = "Ancestry_PC_REGION", 
                           title = "Ancestry principal components by region",
                           metadata,
                           plotdata = rbind(snippets["Project", , drop = FALSE],
                             data.frame(value = "demo.REGION", 
                                        row.names = c("GroupedBy"),
                                        stringsAsFactors = FALSE)))
    }
  }
    
  message('gtxpipe: Writing source display metadata')
  metadata <- metadata[order(as.integer(metadata$number)), ]
  rownames(metadata) <- metadata$number
  metadata$number <- NULL
  ## FIXME, Should warn or truncate if titles >100 characters
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
    suppressMessages(requireNamespace("knitr", quietly = TRUE)) || stop("knitr package not available")
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

  ## gtxpipe.packages should have been cleaned up to match installed.packages by gtxpipe, but
  ## check here in case pipeslave is running in a different environment
  ## Note, any packages not available are silently ignored
  for (package in intersect(getOption('gtxpipe.packages', NULL), rownames(installed.packages()))) {
    if (!require(package, character.only = TRUE, quietly = TRUE)) {
      stop("gtxpipe: fatal error, unable to load package '", package, "'")
    }
  }

  usubjid <- getOption("gtx.usubjid", "USUBJID")
  tmp <- unlist(strsplit(basename(target), ".", fixed = TRUE))
  if (length(tmp) < 2 || !identical(tmp[length(tmp)], "done")) stop ('Bad target "', target, '"')
  job <- paste(tmp[-length(tmp)], collapse = ".") # job is basename with trailing .done stripped off
  rm(tmp)
    
  adataf <- file.path(adir, "analysis-dataset.csv") # analysis data filename
  if (!file.exists(adataf)) stop('Analysis dataset "', adataf, '" does not exist')
  adata0 <- scan(file.path(adir, "analysis-dataset.csv"), sep = "\n", character(0), quiet = TRUE)
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
  ## match cvlist, error if more than 1  FIXME should it be an error if there is no cvlist?
  cvlist <- which(substr(adata0, 1, 10) == '# cvlist: ')
  stopifnot(identical(length(cvlist), 1L))
  cvlist = tokenise.whitespace(substring(adata0[cvlist], 11))
  ## in read.csv should force some settings (stringsAsFactors = TRUE) just in case user options
  ## try to override
  adata <- read.csv(textConnection(adata0[substr(adata0, 1, 1) != "#"]),
                    stringsAsFactors = TRUE)
  adata[[usubjid]] <- as.character(adata[[usubjid]])
  ## apply transformations
  if (nrow(gtxpipe.transformations) > 0) {
    for (idx in 1:nrow(gtxpipe.transformations)) {
      target <- gtxpipe.transformations$targets[idx]
      if (target %in% names(adata)) stop("Transformation overwrites existing variable")
      message(target, " <- ", gtxpipe.transformations$fun[idx]) # debugging message
      adata[[target]] <- eval(parse(text = gtxpipe.transformations$fun[idx]), envir = adata)
    }
  }

  ## ? could sink to .done file, blockassoc should use message() not cat()
  res <- blockassoc(qcall = qcall, data = adata, 
                    minimac = file.path(getOption("gtxpipe.genotypes", "genotypes"), job),
                    threshold.MAF = getOption("gtxpipe.threshold.MAF", 0.01),
                    threshold.Rsq = getOption("gtxpipe.threshold.Rsq", 0.01),
                    threshold.pass = cvlist,
                    ## message specific to intended output file in subsequent write.table
                    message.begin = paste("blockassoc(", file.path(adir, basename(job)), ")", sep = ""))
  ## Require pvalue always, beta and SE if being used for interaction contrasts
  ## Consider writing some comments to this file (blockassoc should return character vector instead of printing messages)
  
  ## Determine which values to write to results file:
  if (getOption("gtxpipe.extended.output", TRUE)) {
    out.columns <- c("SNP", "beta", "SE", "pvalue","Al1","Al2","analysed.Freq1","analysed.Rsq")
  } else {
    out.columns <- c("SNP", "beta", "SE", "pvalue")
  }
  write.table(res[ , out.columns], 
              row.names = FALSE, quote = FALSE, sep = "\t",
              file = gzfile(file.path(adir, paste(job, "out.gz", sep = "."))))
  return(invisible(NULL))
}
        
pipetable <- function(data, filename, title,
                      mdata = data.frame(NULL), number,
                      width = 11.69, height = 8.27) {
  ## FIXME allow option for row names to be suppressed in pdf, requires pass through option to textgrid()

  
  ## number argument is the smallest allowable display number
  number.used <- c(0, as.integer(mdata$number))
  if (missing(number) || number %in% number.used) number <- max(number.used) + 1
  number <- gsub("^ ", "0", format(number, width = getOption("gtxpipe.display.numwidth", 2)))
  number2 <- paste(c(getOption("gtxpipe.display.section", NULL), number), collapse = ".")
  ## Updates to code above should be made here and in pipeplot()

  path <- file.path(getOption("gtxpipe.outputs", "."), paste(number, filename, sep = "_"))
  message('gtxpipe: Writing ', title, ' to "', path, '.[csv|pdf]"')
  dir.create(getOption("gtxpipe.outputs", "."), recursive = TRUE, showWarnings = FALSE) # FIXME throw error if fails

  write.csv(data,
            file = paste(path, "csv", sep = "."),
            row.names = TRUE) # can fail silently if directory does not exist

  pdf(width = width, height = height,
      file = paste(path, "pdf", sep = "."))
  scs <- split.screen(matrix(c(1/11.69, # left margin will be 1" *if* a4 portrait
                               1-1/11.69, # right margin
                               1/8.27, # top margin
                               1-1/8.27 # bottom margin
                               ), nrow = 1))
  screen(scs[1])
  par(family="mono", mar = c(0, 0, 4, 0) + 0.1)
  plot.new()
  plot.window(c(0, 1), c(0, 1))
  textgrid(data)
  mtext(paste("Protocol:", getOption("gtxpipe.protocol", "NA")), side = 3, line = 3, adj = 0)
  mtext("Page 1 of 1", side = 3, line = 3, adj = 1)
  mtext(paste("Population:", "ITT"), side = 3, line = 2, adj = 0) # hard coded to ITT
  mtext(paste("Data as of:", getOption("gtxpipe.date", "NA")), side = 3, line = 2, adj = 1)
  mtext(paste("Source table", number2), side = 3, line = 1, adj = 0.5)
  mtext(title, side = 3, line = 0, adj = 0.5)
  close.screen(scs); dev.off()

  return(rbind(mdata,
               data.frame(Source_File_Name = paste(paste(number, filename, sep = "_"), "pdf", sep = "."), 
                          Display_Category = "PHARMACOGENETIC", 
                          Display_Type = "TABLE", 
                          Display_Number = number2,
                          Title = title,
                          number = number, 
                          stringsAsFactors = FALSE)))
}

pipeplot <- function(plotfun, filename, title,
                     mdata = data.frame(NULL), number,
                     plotdata, plotpar, 
                     width = 8.27, height = 11.69) {

  ## number argument is the smallest allowable display number
  number.used <- c(0, as.integer(mdata$number))
  if (missing(number) || number %in% number.used) number <- max(number.used) + 1
  number <- gsub("^ ", "0", format(number, width = getOption("gtxpipe.display.numwidth", 2)))
  number2 <- paste(c(getOption("gtxpipe.display.section", NULL), number), collapse = ".")
  ## Updates to code above should be made here and in pipetable()
  
  path <- file.path(getOption("gtxpipe.outputs", "."), paste(number, filename, sep = "_"))
  dir.create(getOption("gtxpipe.outputs", "."), recursive = TRUE, showWarnings = FALSE) # FIXME throw error if fails
  
  if (all(capabilities(c("png", "cairo")))) {
    png(type = "cairo", filename = paste(path, "png", sep = "."),
        width = width*300, height = height*300, res = 300)
    message('gtxpipe: Plotting ', title, ' to "', path, '.png"')
  } else if (suppressMessages(requireNamespace("Cairo", quietly = TRUE))) {
    ## This generates an undeclared import warning with R CMD CHECK, not sure why.  FIXME
    Cairo::CairoPNG(file = paste(path, "png", sep = "."),
             width = width*300, height = height*300, res = 300)
    message('gtxpipe: Plotting ', title, ' to "', path, '.png"')
  } else {
    pdf(file = paste(path, "pdf", sep = "."),
        width = width, height = height)
    message('gtxpipe: Plotting ', title, ' to "', path, '.pdf"')
  }

  ## By default, all figures are A4 portrait with the TOP half devoted to plot metadata
  ## Should (?) work even if plotfun subsequently alters e.g. par(mfrow)
  scs <- split.screen(matrix(c(1/8.27, 1/8.27, # left margin
                               1-1/8.27, 1-1/8.27, # right margin
                               0.5+0.25/11.69, 1/11.69, # top margin
                               1-1/11.69, 0.5-0.25/11.69 # bottom margin
                               ), nrow = 2))
  
  screen(scs[1])
  oldpar <- par(family="mono", mar = c(0, 0, 4, 0) + 0.1)
  plot.new()
  plot.window(c(0, 1), c(0, 1))
  if (!missing(plotdata)) {
    textgrid(plotdata)
  } else {
    text(0.5, 0.5, "MISSING PLOT METADATA")
  }
  mtext(paste("Protocol:", getOption("gtxpipe.protocol", "NA")), side = 3, line = 3, adj = 0)
  mtext("Page 1 of 1", side = 3, line = 3, adj = 1)
  mtext(paste("Population:", "ITT"), side = 3, line = 2, adj = 0) # hard coded to ITT
  mtext(paste("Data as of:", getOption("gtxpipe.date", "NA")), side = 3, line = 2, adj = 1)
  mtext(paste("Source figure", number2), side = 3, line = 1, adj = 0.5)
  mtext(title, side = 3, line = 0, adj = 0.5)
  par(oldpar)
  
  screen(scs[2])
  if (!missing(plotpar)) oldpar <- do.call(par, as.list(plotpar)) else oldpar <- par()
  eval.parent(parse(text = plotfun)) # was: eval(..., envir = parent.frame())
  par(oldpar)
  close.screen(scs); dev.off()
  
  return(rbind(mdata,
               data.frame(Source_File_Name = paste(paste(number, filename, sep = "_"), "png", sep = "."), 
                          Display_Category = "PHARMACOGENETIC", 
                          Display_Type = "FIGURE", 
                          Display_Number = number2,
                          Title = title,
                          number = number,
                          stringsAsFactors = FALSE)))
}

pipebed <- function(snplist, flank, outfile) {
  #parse chromosomes and coordinates from snplist
  chr <- vapply(strsplit(snplist, ":"), function(ss) return(ss[1]), character(1))
  pos <- as.integer(vapply(strsplit(snplist, "[:_]"), function(ss) return(ss[2]), character(1)))
  #create dataframe with 4 standard bed columns
  cvbed <- data.frame(chr = chr, s = pos - flank - 1, e = pos + flank,SNP = snplist, stringsAsFactors=FALSE)
  #sort by chromosome and start (since all single base coordinates, no need to sort on end)
  cvbed <- cvbed[order(cvbed$chr,cvbed$s),]
  #re-label columns with standard bed headers
  colnames(cvbed) <- c("#chrom","start","end","SNP")
  #confirmed ok for replicate records in bed file, no need to uniquify for tabix
  suppressWarnings(write.table(cvbed, quote=FALSE, sep = "\t", row.names = FALSE, 
                               file = file.path(outfile),
                               append = FALSE))
}
  

  
