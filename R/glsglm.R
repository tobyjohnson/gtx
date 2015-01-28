gls.approx.logistic <- function(snpdata, leftvar, rightvars) {
  stopifnot(is.data.frame(snpdata$data))
  stopifnot(leftvar %in% names(snpdata$data))
  stopifnot(all(rightvars %in% names(snpdata$data)))
  stopifnot(all(snpdata$data[[leftvar]] %in% 0:1, na.rm = TRUE))
  null <- as.formula(paste(leftvar, paste(rightvars, collapse = "+"), sep = "~"))
  print(null)
  snpdata$data[[paste(leftvar)]]
  p <- 1/(1 + exp(-predict(glm(null, family = "binomial", data = snpdata$data, na.action = na.exclude))))
  snpdata$data$weight <- p*(1-p)
  snpdata$data[[paste(leftvar, "star", sep = "")]] <- (snpdata$data[[leftvar]] - p)/snpdata$data$weight
  return(snpdata)
}
