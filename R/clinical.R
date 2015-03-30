# import should detect .txt.gz extensions and handle automatically
# import might be nice to autodetect SAS exported dates

clinical.import <- function(d, pattern = "^[a-zA-Z][a-zA-Z1-9]*\\.txt",
                            usubjid = getOption("gtx.usubjid", "USUBJID"),
                            verbose = TRUE, convert.YN = TRUE, 
                            only) {
  f <- dir(d, pattern = pattern) # list of files
  if (!missing(only)) {
    only <- tokenise.whitespace(only)
    f <- intersect(f, paste(only, ".txt", sep = ""))
  }
  dl <- list() # empty list of data frames
  for (f1 in f) {
    if (verbose) message("Reading ", f1, " ...")
    tmp <- read.table(file.path(d, f1),
                      header = TRUE, sep = "\t", quote = "",
                      na.strings = ".", comment.char = "",
                      stringsAsFactors = TRUE) # be sure to override any global setting
    if (verbose) message("  Read ", nrow(tmp), " rows and ", ncol(tmp), " columns ")
    if (usubjid %in% names(tmp)) {
      tmp[[usubjid]] <- as.character(tmp[[usubjid]])
      ## Other than usubjid, all character columns are factors
      for (colidx in names(which(sapply(tmp, class) == "factor"))) {
        ## Replace empty strings with NA
        ev <- which(tmp[[colidx]] == "")
        if (length(ev) > 0) {
          ## if (verbose) message(length(ev), " empty values in ", colidx, " converted to NA") # is very verbose
          tmp[[colidx]][ev] <- NA
        }
        ## Convert strings that are all Y/N to logicals
        if (convert.YN && all(grepl("^[YN]$", na.omit(tmp[[colidx]])))) {
          ## if (verbose) message(colidx, " Y/N converted to logical") # is very verbose
          tmp[[colidx]] <- tmp[[colidx]] == "Y"
        }
        ## FIXME grepl for "[0-9]{2}[A-Za-z]{3}[0-9]{4}", try as.Date and if works on >99% convert to Date
      }
      dl[[sub("\\.txt$", "", f1)]] <- tmp
    } else {
      warning("Data in file '", f1, "' was not imported because there was no usubjid column '", usubjid, "'")
    }
  }
  return(dl)
}

## pgx.models <- function(model, deps, fun, groups, contrasts) {
##   models <- getOption("pgx.models", NULL)
##   if (is.null(models)) {
##     modelfile <- getOption("pgx.modelfile", "config/pgx.models")
##     message("Reading pgx.models from [ ", modelfile, " ]")
##     models <- tryCatch(suppressWarnings(read.csv(modelfile)),
##                        error = function(e) return(data.frame(model = NULL, deps = NULL, fun = NULL, groups = NULL, contrasts = NULL)))
##   }
##   if (missing(model)) {
##     options(pgx.models = models)
##     return(models)
##   }
##   if (missing(deps) || missing(fun) || missing(groups
  
##   options(pgx.models = models)
##   return(models)
## }





derivation.add <- function(derivations, targets, types, deps, data, fun, aept.list, verbose = TRUE) {
  ## If no existing derivations, create an empty data frame with required columns
  if (missing(derivations)) {
    derivations <- data.frame(targets = character(0), types = character(0),
                              deps = character(0), data = character(0), fun = character(0),
                              stringsAsFactors = FALSE)
  }
  if (missing(targets)) stop("targets is a required argument")
  targets <- as.character(targets)
  if (!missing(types) && !missing(deps) && !missing(data) && !missing(fun)) {
    types <- as.character(types)
    stopifnot(identical(length(targets), length(types)))
    ## parse on whitespace split, check length(typev)==length(targetv) elementwise
    deps <- as.character(deps)
    stopifnot(identical(length(targets), length(deps)))
    data <- as.character(data)
    stopifnot(identical(length(targets), length(data)))
    fun <- as.character(fun)
    stopifnot(identical(length(targets), length(fun)))
    ## All arguments supplied, literally add to existing derivations and return
    message("Adding derivations for targets : ", targets)
    return(rbind(derivations,
                 data.frame(targets = targets, types = types, deps = deps, data = data, fun = fun,
                            stringsAsFactors = FALSE)))
  }

  message("Guessing derivations for targets : ", targets)

  ## Attempt to guess using v1 rules
  targetv <- unlist(strsplit(targets, "\\s+"))
  sapply(targetv, function(t1) {
    t1s <- unlist(strsplit(t1, "\\."))
    if (length(t1s) != 3) stop("Cannot guess derivation for target : ", t1)
    return(t1s)
  })
  stop("Intelligent construction not yet implemented")
}
  
derive1 <- function(datalist, targets, types, deps, data, fun) {
  usubjid <- getOption("gtx.usubjid", "USUBJID")
  stopifnot(all(deps %in% names(datalist)))
  targetv <- tokenise.whitespace(targets)
  typev <- rep(tokenise.whitespace(types), length.out = length(targetv)) # recycle or truncate silently
  output1 <- tryCatch(with(datalist, eval(parse(text = data))),
                      error = function(e) stop(data, " failed with error ", e))
  stopifnot(is.data.frame(output1))
  stopifnot(usubjid %in% names(output1))
  u <- output1[ , usubjid, drop = TRUE]
  uu <- unique(u)
  ##  message('Deriving ', paste(targetv, collapse = ", "), ' for ', length(uu), ' subjects')
  ## TO DO: Elegantly handle situation with 0 subjects
  ##  stopifnot(identical(length(targetv), length(eval(parse(text = type)))))
  foo <- data.frame(uu, 
                    do.call(rbind, lapply(uu,
                                          FUN = function(u1) with(output1[which(u == u1), , drop = FALSE], eval(parse(text = fun))))),
                    stringsAsFactors = FALSE)
  ## try eval(...,envir=...)
  if (ncol(foo) != length(targetv) + 1) {
    stop('derivation function "', fun, '" returned length ', ncol(foo)-1, ' for targets ', paste(targetv, collapse = ", "))
  }
  ##  stopifnot() ???
  ## cat(identical(ncol(foo), length(targetv) + 1), length(targetv) + 1, ncol(foo), names(foo), "\n")
  names(foo) <- c(usubjid, targetv)
  ## coerce to types
  for (idx in 1:length(typev)) {
    if (identical(typev[idx], "factor")) { # special case
      storage.mode(foo[, idx + 1]) <- "character"
      foo[ , idx + 1] <- factor(foo[ , idx + 1])
    } else if (identical(typev[idx], "Date")) { # special case, Date must be like 21MAR1975, which is default for SAS proc export
      storage.mode(foo[, idx + 1]) <- "character"
      foo[ , idx + 1] <- as.Date(foo[ , idx + 1], "%d%b%Y")
    } else {
      ## should check typev[idx] is a valid storage mode
      storage.mode(foo[, idx + 1]) <- typev[idx]
    }
  }
  return(foo)
}

clinical.derive <- function(datalist, derivations, verbose = TRUE, only) {
  usubjid <- getOption("gtx.usubjid", "USUBJID")
  stopifnot(all(c("targets", "types", "deps", "data", "fun") %in% names(derivations)))
  if (missing(only)) only <- derivations$targets
  only <- tokenise.whitespace(only)
  pgx <- subset(datalist[["pop"]], select = usubjid) # check this is robust
  ## if (verbose) cat("Initial N = ", nrow(pgx), "\n", sep = "")
  for (idx in 1:nrow(derivations)) {
    targetv <- tokenise.whitespace(derivations[idx, "targets"])
    if (length(intersect(targetv, only)) > 0) {
      if (verbose) message("Deriving ", paste(targetv, collapse = ", "))
      pgx <- merge(pgx,
                   with(derivations[idx, ], derive1(datalist, targets = targets, types = types, deps = deps, data = data, fun = fun)),
                   all = TRUE)
      ## if (verbose) cat("N = ", nrow(pgx), "\n", sep = "")
    }
  }
  ## warn/error if setdiff(only, names(pgx)) nonempty
  pgx <- subset(pgx, select = intersect(c(usubjid, only), names(pgx)))
  return(pgx)
}



#derivations.standard <- derivation.add(derivations.standard, targets = 'lab.ALT.max lab.ALT.tt3uln lab.ALT.tt5uln lab.ALT.tmax', types = 'i#nteger',#
#               deps = 'lab', data = 'subset(lab, LBTEST == "Alanine Amino Transferase" & ATTYPE == "On-therapy" & !is.na(LBSTRESN) & !is.na#(LBSTRHI))',
#               fun = 'c(safe(max, LBSTRESN/LBSTNRHI), safe(min, LBACTDY[LBSTRESN > 3*LBSTNRHI]), safe(min, LBACTDY[LBSTRESN > 5*LBSTNRHI]),# safe(max, LBACTDY))')

