demographics <- function(object, by, style, digits) UseMethod("demographics")

demographics.numeric <- function(object, by, style, digits) {
  if (missing(style) || is.null(style)) style <- "Mean (SD)" # getOption?
  if (missing(digits) || is.null(digits)) digits <- getOption("demographics.digits", 2)
  ## nt = numeric table of summaries
  nt <- t(vapply(X = by, FUN = function(by1) {
    object1 <- object[which(by1)]
    m <- is.na(object1)
    if (all(m)) return(c(Missing = sum(m),
                         Mean = NA, Median = NA, 
                         SD = NA, Q1 = NA, Q3 = NA,
                         Min = NA, Max = NA))
    ## o1q = object1 quantiles 
    o1q <- unname(quantile(object1[!m], 0:4/4))
    return(c(Missing = sum(m),
             Mean = mean(object1[!m]), Median = o1q[3],
             SD = sd(object1[!m]), Q1 = o1q[2], Q3 = o1q[4],
             Min = o1q[1], Max = o1q[5]))
  },
                 double(8)))
  ## force colnames(nt) <-
  ## FOLLOWING CODE ASSUMES MISSING IS COLUMN 1
  ## ntf = nt formatted
  ntf <- vapply(1:ncol(nt), function(icol) {
    return(formatC(nt[ , icol], digits = digits, format = "f", flag = "-"))
  }, character(nrow(nt)))
  if (!is.array(ntf)) ntf <- matrix(ntf, nrow = 1)
  dimnames(ntf) <- dimnames(nt)
#  ntf <- apply(nt[ , -1, drop = FALSE], 2, format, digits = getOption("demographics.digits", 4),
#                scientific = FALSE)
  ## Add extra columns derved from existing columns
  ntf <- cbind(ntf,
               IQR = paste(ntf[ , "Q1"], ntf[ , "Q3"], sep = "-"),
               Range = paste(ntf[ , "Min"], ntf[ , "Max"], sep = "-"))
  ## dt = demographics table
  dt <- vapply(X = 1:nrow(ntf), FUN = function(irow) {
    style1 <- style
    for (icol in 1:ncol(ntf)) style1 <- sub(colnames(ntf)[icol], ntf[irow, icol], style1, fixed = TRUE)
    return(c(style1, nt[irow, "Missing"]))
  },
               character(2))
  rownames(dt) <- c(style, "Missing, N")
  # if (getOption("demographics.omit0rows", TRUE) && all(dt["Missing, N", ] == "0")) dt <- dt[1, , drop = FALSE]
  colnames(dt) <- names(by)
  return(dt)
}

demographics.factor <- function(object, by, style, digits) {
  if (missing(style) || is.null(style)) style <- "N (pc0%)" # getOption?
  if (missing(digits) || is.null(digits)) digits <- getOption("demographics.digits", 2)
  ## ct = count table of summaries
  ct <- vapply(X = by, FUN = function(by1) {
    object1 <- object[which(by1)]
    m <- is.na(object1)
    return(c(table(object1[!m]), Missing = sum(m)))
  },
               integer(length(levels(object)) + 1))
  ctf <- apply(ct, 2, function(ct1) {
    ## Note that pc must come after pc[012] for the pattern matching to work
    tmp <- data.frame(N = formatC(ct1, digits = 0, format = "d"),
                      pc0 = prettypc(ct1, digits = 0),
                      pc1 = prettypc(ct1, digits = 1),
                      pc2 = prettypc(ct1, digits = 2),
                      pc = prettypc(ct1, digits = digits),
                      stringsAsFactors = FALSE)
    return(sapply(1:nrow(tmp), function(irow) {
      style1 <- style
      for (icol in 1:ncol(tmp)) style1 <- sub(colnames(tmp)[icol], tmp[irow, icol], style1, fixed = TRUE)
      return(style1)
    }))})
  if (nrow(ctf) > 0) rownames(ctf) <- paste(rownames(ct), gsub("pc[0-9]*%*", "%", style), sep = ", ")
  ## remove rows of all zeros
  if (getOption("demographics.omit0rows", TRUE)) ctf <- ctf[!apply(ct == 0, 1, all), , drop = FALSE]
  ## coerce to character?

  return(ctf)
  ## add nice pretty pcs
}

demographics.Surv <- function(object, by, style, digits) {
  if (missing(style) || is.null(style)) style <- "Median (95% CI)" # getOption?
  if (missing(digits) || is.null(digits)) digits <- getOption("demographics.digits", 2)
  ## nt = numeric table of summaries
  nt <- t(vapply(X = by, FUN = function(by1) {
    object1 <- object[which(by1)]
    m <- is.na(object1)
    if (all(m)) return(c(Missing = sum(m),
                         events = NA, median = NA,
                         "0.95LCL" = NA, "0.95UCL" = NA))
    ## o1f = object 1 survival fit summary
    o1f <- tryCatch(summary(survfit(object1[!m] ~ 1))$table[c("events", "median", "0.95LCL", "0.95UCL")],
                    error = function(e) return(c("events" = NA, "median" = NA, "0.95LCL" = NA, "0.95UCL" = NA)))
    return(c(Missing = sum(m), o1f))
  },
                 double(5)))
  colnames(nt) <- c("Missing", "Events", "Median", "95% LCL", "95% UCL")
  ## FOLLOWING CODE ASSUMES MISSING IS COLUMN 1
  ## ntf = nt formatted
  ntf <- vapply(1:ncol(nt), function(icol) {
    return(formatC(nt[ , icol], digits = digits, format = "f", flag = "-"))
  }, character(nrow(nt)))
  if (!is.array(ntf)) ntf <- matrix(ntf, nrow = 1)
  dimnames(ntf) <- dimnames(nt)
  ## Add extra columns derved from existing columns
  ntf <- cbind(ntf,
               "95% CI" = paste(ntf[ , "95% LCL"], ntf[ , "95% UCL"], sep = "-"))
  ## dt = demographics table
  dt <- vapply(X = 1:nrow(ntf), FUN = function(irow) {
    style1 <- style
    for (icol in 1:ncol(ntf)) style1 <- sub(colnames(ntf)[icol], ntf[irow, icol], style1, fixed = TRUE)
    return(c(style1, nt[irow, "Missing"]))
  },
               character(2))
  rownames(dt) <- c(style, "Missing, N")
  colnames(dt) <- names(by)
  return(dt)
}

demographics.logical <- function(object, by, style, digits) {
  ## style and digits not used
  ## tmt = truth and missingness table
  tmt <- vapply(X = by, FUN = function(by1) {
    object1 <- object[which(by1)]
    m <- is.na(object1)
    return(c(N = sum(object1[!m] == TRUE), Missing = sum(m)))
  },
                integer(2))
  ## remove missing row if all zero
  if (getOption("demographics.omit0rows", TRUE)) if (all(tmt[nrow(tmt), ] == 0)) tmt <- tmt[-nrow(tmt), , drop = FALSE]
  ## coerce to character?
  return(tmt)
  ## add nice pretty pcs
}


demographics.default <- function(object, by, style, digits) {
  ## cmt = class and missingness table
  cmt <- vapply(X = by, FUN = function(by1) {
    object1 <- object[which(by1)]
    m <- is.na(object1)
    return(c(Class = class(object), Missing = format(sum(m))))
  },
                character(2))
  return(cmt)
  ## add nice pretty pcs
}
  


demographicsCheckBy <- function(object, by) {
  if (is.list(by) && length(by) > 0) {
    ## Check whether elements of by are named
    if (is.null(names(by))) {
      warning("by argument is unnamed, using \"Unnamed group 1\" etc.")
      names(by) <- paste("Unnamed group", 1:length(by))
    }
    bybad <- is.na(names(by))
    if (any(bybad)) {
      warning("by argument has missing names, using \"Unnamed group 1\" etc.")
      names(by)[which(bybad)] <- paste("Unnamed group", 1:sum(bybad))
    }
    ## Check elements of by are logical vectors with length nrow(object)
    byclass <- lapply(by, class)
    bybad <- byclass != "logical"
    if (any(bybad) > 0) warning("ignoring ", paste("by[[\"", names(by[bybad]), "\"]]", sep = "", collapse = ", "),
                                " because not of class logical")
    by <- by[!bybad]
    bylen <- sapply(by, length)
    bybad <- bylen != nrow(object)
    if (any(bybad) > 0) warning("ignoring ", paste("by[[\"", names(by[bybad]), "\"]]", sep = "", collapse = ", "),
                                " because length not equal to nrow(object)")
    by <- by[!bybad]
    ## Replace any NA elements with FALSE
    ## by <- lapply(by, function(x) return(is.na(x) & x))
  }
  if (identical(length(by), 0L)) by <- list("All subjects" = rep(TRUE, nrow(object))) # failsafe
  if (!is.list(by)) by <- list("All subjects" = rep(TRUE, nrow(object))) # failsafe
  return(by)
}


demographics.data.frame <- function(object, by, style, digits) {
  ## check by argument, modify if necessary
  if (missing(by)) by <- list("All subjects" = rep(TRUE, nrow(object))) # failsafe
  if (is.list(by)) {
    bylist <- by
  } else if (is.character(by)) {
    ## build list, assuming by contains names of logical or factor elements of object
    bybad <- is.na(match(by, names(object)))
    if (any(bybad)) warning("by value(s) ", paste(by[bybad], collapse = ", "),
                            " not in names(object)")
    by <- by[!bybad]
    if (length(by) > 0) {
      bym <- do.call(cbind, lapply(by, function(by1) {
        if (identical(class(object[[by1]]), "logical")) {
          tmp <- matrix(object[[by1]], ncol = 1)
          colnames(tmp) <- by1
          return(tmp)
        } else if (identical(class(object[[by1]]), "factor")) {
          tmp <- vapply(levels(object[[by1]]), function(by1v) return(object[[by1]] == by1v),
                        logical(nrow(object)))
          colnames(tmp) <- paste(by1, levels(object[[by1]]), sep = "=")
          return(tmp)
        } else {
          warning("by argument contains invalid value \"", by1, "\", object[[", by1, "]] is neither logical nor factor")
          return(NULL)
        }
      }))
      ## remove columns from object used to construct by
      object <- object[ , -match(by, names(object)), drop = FALSE] #bugfix when reduced object only has one column left
      ## remake by as a list
      bylist <- lapply(colnames(bym), function(bym1) bym[ , bym1])
      names(bylist) <- colnames(bym)
    } else {
      bylist <- list()
    }
  } else {
    warning("by argument of incorrect class, must be list or character")
    bylist <- list("All subjects" = rep(TRUE, nrow(object))) # failsafe
  }

  bylist <- demographicsCheckBy(object, bylist)

  ## Need something to pass into specialised functions which provide their own defaults
  if (missing(style)) style <- list(NULL)
  if (missing(digits)) digits <- list(NULL) 

  ## demographics table n(N) row
  n <- vapply(X = bylist, FUN = function(by1) return(sum(by1 == TRUE, na.rm = TRUE)), integer(1))
  dtn <- matrix(format(n, scientific = FALSE), nrow = 1)
  rownames(dtn) <- "N"

  ## demographics table all items
  dta <- do.call(rbind,
                 lapply(setdiff(names(object), getOption("gtx.usubjid", "USUBJID")), 
                        function(oname) {
                          ## Note, S3 dispatching automatically dispatches to demographics.numeric for
                          ## objects of class integer, but the style and digits lookups need to handle this
                          ## as a special case
                          style1 <- style[[oname]]
                          if (is.null(style1)) style1 <- style[[class(object[[oname]])]]
                          if (is.null(style1) && class(object[[oname]]) == "integer") style1 <- style[["numeric"]]
                          digits1 <- digits[[oname]]
                          if (is.null(digits1)) digits1 <- digits[[class(object[[oname]])]]
                          if (is.null(digits1) && class(object[[oname]]) == "integer") digits1 <- digits[["numeric"]]
                          dt <- demographics(object[[oname]], by = bylist, style = style1, digits = digits1)
                          if (nrow(dt) > 0) {
                            odesc <- getOption("clinical.descriptors", list(NULL))[[oname]]
                            if (is.null(odesc) || is.na(odesc)) odesc <- oname
                            bstring <- paste(rep.int(".", nchar(odesc)), collapse = "")
                            rownames(dt) <- paste(c(odesc, rep.int(bstring, nrow(dt) - 1)), rownames(dt), sep = ", ")
                          }
                          ## Alternatively support the style of a separate row for the variable name...?
                          return(dt)
                        }))
  dta <- rbind(dtn, dta)
  colnames(dta) <- names(bylist)
  ## if stripping empty columns
  if (getOption("demographics.omit0groups", TRUE)) dta <- dta[ , which(n > 0), drop = FALSE]
  return(dta)
}

prettypc <- function(x, digits = 0) {
  x <- as.integer(x)  
  pc <- 100*x/sum(x, na.rm = TRUE)
  pcr <- round(pc, digits = digits)
  pcf <- formatC(pcr, digits = digits, format = "f", flag = "-") # add width = digits + 4 ?
  pcf[pcr == 0 & pc > 0] <- paste("<", formatC(0.+10^-digits, digits = digits, format = "f", flag = "-"), sep = "") # <1
  pcf[pcr == 100 & pc < 100] <- paste(">", formatC(100.-10^-digits, digits = digits, format = "f", flag = "-"), sep = "") # >99
  return(pcf)
}

#hlines <- function(x) {
#  if (!("demographicstable" %in% class(x))) return (NULL)
#  return(unique(sort(c(-1, 0, 1, nrow(x), which(x$variable[-nrow(x)] != x$variable[-1])))))
#}

  

#n2pc <- function(x, digits = 1, 
#                 excludeval = c("mean", "sd", "min", "max", "median", "0.95LCL", "0.95UCL")) {
#  if (!("demographicstable" %in% class(x))) return (x)
#  excludecol <- match(c("variable", "value"), colnames(x))
#  excluderow <- x$value %in% excludeval
#  tt <- cbind(x[ , excludecol, drop = FALSE],
#              do.call("cbind", lapply(setdiff(1:ncol(x), excludecol), function(ii) cbind(x[ , ii, drop = FALSE], pc = ifelse(!excluderow, paste(round(100*x[ , ii]/x[1 , ii], digits = d#igits), "%", sep = ""), NA)))))
#  class(tt) <- unique(c("demographicstable", class(tt)))
#  return(tt)  
#}



#makesubsets <- function(dd, subsets) {
#  dd <- as.data.frame(dd)
#  if (is.null(subsets)) subsets <- data.frame(All = rep(TRUE, nrow(dd)))
#  if (is.character(subsets) && all(subsets %in% names(dd))) subsets <- subset(dd, select = subsets)
#  subsets <- as.data.frame(subsets)
#  stopifnot(nrow(subsets) == nrow(dd))
#  for (jj in which(sapply(1:ncol(subsets), function(ii) is.factor(subsets[ , ii])))) {
#    for (ll in levels(subsets[ , jj])) subsets[[make.names(ll)]] <- subsets[ , jj] == ll
#  }
#  subsets <- subsets[ , which(sapply(1:ncol(subsets), function(ii) is.logical(subsets[ , ii]))), drop = FALSE]
#  subsets[is.na(subsets)] <- FALSE
#  return(subsets)
#}

#options(list(sort.tables.by.frequency,
#             max.table.length,
#             use.missing.always,
#             numeric.mean/median/etc
#             surv.)


#factornames <- function(dd) {
#  dd <- as.data.frame(dd)
#  names(dd)[sapply(1:ncol(dd), function(ii) is.factor(dd[ , ii]))]
#}
#numericnames <- function(dd) {
#  dd <- as.data.frame(dd)
#  names(dd)[sapply(1:ncol(dd), function(ii) is.numeric(dd[ , ii]))]
#}
#makesubsets <- function(dd, subsets) {
#  dd <- as.data.frame(dd)
#  if (is.null(subsets)) subsets <- data.frame(All = rep(TRUE, nrow(dd)))
#  if (is.character(subsets) && all(subsets %in% names(dd))) subsets <- subset(dd, select = subsets)
#  subsets <- as.data.frame(subsets)
#  stopifnot(nrow(subsets) == nrow(dd))
#  for (jj in which(sapply(1:ncol(subsets), function(ii) is.factor(subsets[ , ii])))) {
#    for (ll in levels(subsets[ , jj])) subsets[[make.names(ll)]] <- subsets[ , jj] == ll
#  }
#  subsets <- subsets[ , which(sapply(1:ncol(subsets), function(ii) is.logical(subsets[ , ii]))), drop = FALSE]
#  subsets[is.na(subsets)] <- FALSE
#  return(subsets)
#}
#factortable <- function(dd, 
#                        subsets = NULL,
#                        select = factornames(dd), 
#                        na.rm = TRUE) {
#  dd <- as.data.frame(dd)
#  subsets <- makesubsets(dd, subsets)
#  stopifnot(nrow(subsets) == nrow(dd))
#  tt <- rbind(data.frame(factor = NA, level = NA),
#              do.call(rbind, lapply(select, function(nn) {
#                ll = levels(dd[[nn]])
#                data.frame(factor = rep(nn, length(ll)), level = ll, 
#                           stringsAsFactors = FALSE)})))
#  if (nrow(tt) > 1) {
#    for (ss in names(subsets)) {
#      xx <- which(subsets[[ss]])
#      tt[[ss]] <- c(length(xx), sapply(2:nrow(tt), function(ii) {
#        sum(dd[[tt$factor[ii]]][xx] == tt$level[ii], na.rm = na.rm)}))
#    }
#  }
#  return(tt)
#}
#srange <- function(xx, na.rm = FALSE, digits = 1) {
#  rr <- round(range(xx, na.rm = na.rm), digits = digits)
#  return(paste(rr[1], rr[2], sep = "-"))
#}
## no nice way to pass extra ags like digits into funs
#numerictable <- function(dd, subsets = data.frame(all = rep(TRUE, nrow(dd))), 
#                         select = numericnames(dd), funs = c("min", "median", "max"), 
#                         na.rm = TRUE) {
#  dd <- as.data.frame(dd)
#  subsets <- makesubsets(dd, subsets)
#  stopifnot(nrow(subsets) == nrow(dd))
#  tt <- expand.grid(fun = funs, variable = select, 
#                    stringsAsFactors = FALSE)[ , c(2, 1)]
#  for (ss in names(subsets)) {
#    xx <- which(subsets[[ss]])
#    tt[[ss]] <- sapply(1:nrow(tt), function(ii) {
#      do.call(tt$fun[ii], list(dd[[tt$variable[ii]]][xx], na.rm = na.rm))})
#  }
#  return(tt)
#}
#n2pc <- function(tt, exclude = 1:2) {
#  tt <- as.data.frame(tt)
#  return(cbind(tt[ , exclude, drop = FALSE],
#               do.call("cbind", lapply(setdiff(1:ncol(tt), exclude), function(ii) cbind(tt[ , ii, drop = FALSE], pc = 100*tt[ , ii]/tt[1 , ii])))))
#}
#forcebind <- function(tt1, tt2) {
#  names(tt1)[1:2] <- 1:2
#  names(tt2)[1:2] <- 1:2
#  tt <- rbind(tt1, tt2)
#  names(tt)[1:2] <- c("variable", "value")
#  return(tt)
#}
