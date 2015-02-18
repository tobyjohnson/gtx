sanitise.whitespace <- function(tt) return(gsub(" $", "", gsub("^ ", "", gsub(" +", " ", tt))))

tokenise.whitespace <- function(tt) {
  tmp <- unlist(strsplit(unname(tt), '\\s+'))
  return(tmp[tmp != ""])
}

text2factor <- function(t) return(as.factor(ifelse(t != "", t, NA)))

## GSK format, but not strinctly CDER/CDISC compliant

## because USUBJID must be 1:1 with subjects across all datasets
## supporting a single submission.  E.g. VEG105192 is randomised study
## and VEG107769 is an open label rollover, given subject would need
## to have same USUBJID in both datasets.  Similarly ANORO priming +
## followup trials...

my.usubjid <- function(studyid, subjid) {
  return(paste(studyid, gsub(" ", "0", format(subjid, digits = 0, width = 7, scientific = FALSE)), sep = "."))
}
