fixNA <- function(tt, na.string = "NA") {
  if (length(dim(tt)) == 1 || length(dim(tt)) == 2) rownames(tt)[is.na(rownames(tt))] <- na.string
  if (length(dim(tt)) == 2) colnames(tt)[is.na(colnames(tt))] <- na.string
  return(tt)
}
