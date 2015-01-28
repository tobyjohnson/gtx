grade <- function(x, threshvals, strict = TRUE) {
  g <- ifelse(!is.na(x), 0, NA)
  threshvals <- sort(threshvals)
  if (strict) {
    for (ii in 1:length(threshvals)) g[x > threshvals[ii]] <- ii
  } else {
    for (ii in 1:length(threshvals)) g[x >= threshvals[ii]] <- ii
  }
  return(g)
}

binarise <- function(is0, is1) {
  stopifnot(identical(length(is0), length(is1)))
  return(ifelse(is0 & !is1, 0, ifelse(is1 & !is0, 1, NA)))
}

pgx.exposure <- function(d, exposure,
                         estart = "exstdtR", eend = "exendtR",
                         subjid = "USUBJID") {
                         
  stopifnot(all(c(subjid) %in% names(d)))
  stopifnot(all(c(subjid, estart, eend) %in% names(exposure)))

  ## Calculate total days of exposure (de) for each subject
  de <- vapply(d[[subjid]], function(ii) {
    mincover(exposure[[estart]][exposure[[subjid]] == ii],
             exposure[[eend]][exposure[[subjid]] == ii])
  }, double(1))

  ## Calculate start of exposure (se) and end of exposure (ee)
  #se <- aggregate(exposure[[estart]], list(exposure[[subjid]]),
  #                FUN = min)
  #ee <- aggregate(exposure[[eend]], list(exposure[[subjid]]),
  #                FUN = max)
  ## Should merge these in...
  
  ## Make data frame of derived variables
  dv <- data.frame(d[[subjid]],
                   exposure = de,
                   stringsAsFactors = FALSE)
  names(dv) <- c(subjid, "exposure")
  return(dv)
}

pgx.trtreat <- function(measure, exposure,
                        mtime, winlen,
                        estart = "exstdtR", eend = "exendtR",
                        subjid = "USUBJID") {
  ## Time in Relation to TREATment
  
  stopifnot(all(c(subjid, mtime) %in% names(measure)))
  stopifnot(all(c(subjid, estart, eend) %in% names(exposure)))
  
  ## Calculate start of exposure (se) and end of exposure (ee)
  ## suppressWarnings allows Inf and -Inf to be returned
  se <- suppressWarnings(aggregate(exposure[[estart]], list(exposure[[subjid]]),
                                   FUN = min, na.rm = TRUE))
  ee <- suppressWarnings(aggregate(exposure[[eend]], list(exposure[[subjid]]),
                                   FUN = max, na.rm = TRUE))
  ## No need for explicit list of ever exposed (ne) subjects

  ## Match to measure by subjid; mne = matched never exposed
  mse <- finitise(se[match(measure[[subjid]], se[ , 1]), 2])
  mee <- finitise(ee[match(measure[[subjid]], ee[ , 1]), 2] + winlen)
  mne <- !(measure[[subjid]] %in% exposure[[subjid]])
  ## Calculate result
  return(ifelse(is.na(measure[[mtime]]), NA,
                ifelse(mne, "No-therapy",
                       ifelse(is.na(mse) | is.na(mee), NA, 
                              ifelse(measure[[mtime]] <= mse, "Pre-therapy",
                                     ifelse(measure[[mtime]] <= mee, "On-therapy",
                                            "Post-therapy"))))))
}
                                          
pgx.endpoints <- function(d, measure, exposure,
                mtime, mvalue, threshval,
                msign = c("greater", "less"), mstrict = c(TRUE, FALSE), 
                estart = "exstdtR", eend = "exendtR",
                dby = NULL, subjid = "USUBJID") {

  stopifnot(all(c(subjid, dby) %in% names(d)))
  stopifnot(all(c(subjid, mtime, mvalue) %in% names(measure)))
  stopifnot(all(c(subjid, estart, eend) %in% names(exposure)))

  if (identical(msign, "greater")) {
    flip <- 1 ## do nothing
    flipname <- "max"
  } else if (identical(msign, "less")) {
    flip <- -1 ## used to flip measure[[mvalue]], v1 and v2
    flipname <- "min"
  } else {
    stop("msign must be either \"greater\" or \"less\"")
  }

  m0 <- measure[!is.na(measure[[mvalue]]), ] # no missing values
  if (identical(mstrict, TRUE)) {
    m1 <- m0[!is.na(m0[[mtime]]) & m0[[mvalue]]*flip > threshval*flip, ] # and [flipped] value > threshval
  } else if (identical(mstrict, FALSE)) {
    m1 <- m0[!is.na(m0[[mtime]]) & m0[[mvalue]]*flip >= threshval*flip, ] # and [flipped] value >= threshval
  } else {
    stop("mstrict must be either TRUE or FALSE")
  }

  ## Calculate maximum value for each subject
  mm <- aggregate(m0[[mvalue]]*flip, list(m0[[subjid]]), FUN = max) # max [flipped] value
  mm <- mm[match(d[[subjid]], mm[ , 1]), 2]*flip
  ## Calculate days of exposure until first [flipped] value >[>=] threshval for each subject
  t1 <- aggregate(m1[[mtime]], list(m1[[subjid]]), FUN = min)
  t1 <- vapply(1:nrow(t1), function(ii) {
    mincover(exposure[[estart]][exposure[[subjid]] == t1[ii, 1]],
             exposure[[eend]][exposure[[subjid]] == t1[ii, 1]],
             trunc.end = t1[ii, 2])
  }, double(1))[match(d[[subjid]], t1[ , 1])]
  ## Calculate days for 'sufficient exposure', by stratum if dby is not NULL
  if (is.null(dby)) {
    dse <- rep(median(t1, na.rm = TRUE), length(d[[subjid]]))
  } else {
    t1by <- aggregate(t1, list(d[[dby]]), FUN = median, na.rm = TRUE)
    dse <- t1by[match(d[[dby]], t1by[ , 1]), 2]
  }
  ## Calculate total days of exposure for each subject
  de <- vapply(d[[subjid]], function(ii) {
    mincover(exposure[[estart]][exposure[[subjid]] == ii],
             exposure[[eend]][exposure[[subjid]] == ii])
  }, double(1))
  ## Make data frame of derived variables
  dv <- data.frame(d[[subjid]],
                   !is.na(t1) | de >= dse,
                   mm,
                   ifelse(!is.na(t1), 1, 0),
                   ifelse(!is.na(t1), t1, de),
                   stringsAsFactors = FALSE)
  
  names(dv) <- c(subjid, paste(mvalue, c("exposed", flipname, "event", "time"), sep = "."))
  return(dv)
  ## binary case/control analysis could be event | exposed
  ## QT analysis is minmax | exposed
  ## TTE is event, time
}

## gradeplot <- function(time, grade, subject) {
##   stopifnot(identical(length(grade), length(time)))
##   stopifnot(identical(length(subject), length(time)))
##   glev <- unique(sort(grade))
##   cumcount <- sapply(glev, function(gmin) {
##     return(sapply(1:length(time), function(idx) {
##       length(unique(subject[which(grade[1:idx] >= gmin)]))
##     }))
##   })
##   plot.new()
##   plot.window(range(time), c(0, max(cumcount)*1.25))
##   axis(1); axis(2, las = 1); box()
##   gcol <- rev(rainbow(length(glev), start = 0, end = 1/3))
##   for (idx in 1:length(glev)) {
##     lines(time, cumcount[ , idx], col = gcol[idx])
##   }
##   title(xlab = "Time", ylab = "Cumulative incidence")
##   legend("top", horiz = TRUE, lwd = 1, col = gcol, 
##        legend = paste("grade ", glev, "+", sep = ""))
## }

