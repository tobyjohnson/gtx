###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
#' Retrieve data for forest plots
#'
#' @export
retrieve_data <- function(id) {

  rs_data <- phewas.data(rs=id)

  rs_data <- within(rs_data, {
    mlog10p <- pmin(30, -log10(pval))
    category <- ifelse(!is.na(entity), "eQTL", "disease")
  })

  rs_data$category <- as.factor(rs_data$category)
  rs_data <- rs_data[order(rs_data$category, rs_data$pval), ]

  cat <- unique(rs_data$category[!is.na(rs_data$category)])

  sub_disease <- subset(rs_data, category %in% cat[1])
  sub_disease$row_num <- 1:nrow(sub_disease)
  sub_eqtl    <- subset(rs_data, category %in% cat[2])
  dis_res <- nrow(sub_disease)
  eqtl_res <- nrow(sub_eqtl)
  comm <- paste("There are ", dis_res, " disease focussed analyses, and ",  eqtl_res, " eQTL focussed analyses associated with ", id, ".", sep = "")
  print(comm)
  res <- list(disease = sub_disease, eqtl = sub_eqtl)
}

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

#' @export
fplot_disease <- function(res, n = 20) {

  # Now paint a Forest Plot of the 'n' most statistically significant disease associations.
  # Turn off warnings
  options(warn=-1)                 # To turn back on  ...  options(warn=0)
 
  fplot_dis <- res$disease[1:n, c(1, 12, 10, 11)]
  fplot_dis$L95 <- fplot_dis$beta - 1.96 * fplot_dis$se
  fplot_dis$U95 <- fplot_dis$beta + 1.96 * fplot_dis$se
  fplot_dis$analysis <- factor(fplot_dis$analysis, levels = rev(fplot_dis$analysis))

  date = date
  title <- paste('Simple Forest Plot of', myquery, 'Disease Associations :', date, sep=" ")

  fp <- ggplot(data = fplot_dis, aes(y = analysis, x = beta, xmin = L95, xmax = U95)) +
  geom_errorbarh(height = 0.3, col = 'royalblue', alpha = 0.6) +
  geom_point(aes(size = round(-log10(pval), 2)), col = 'royalblue') +
  geom_vline(xintercept = 0, color = 'grey50', linetype = 'dashed') +
  theme_bw() +
  xlab("beta estimate (95% CI)") + 
  ylab(" ") +
  labs(title = title)       
  
  suppressMessages(ggplotly(fp))
}
  

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

#' @export
fplot_eqtl <- function(res, n = 20) {

  # Now paint a Forest Plot of the n most statistically significant eQTL associations.
  # Turn off warnings
  options(warn=-1)                 # To turn back on  ...  options(warn=0)
 
  fplot_eqtl <- res$eqtl[1:n , c(1, 9, 12, 10, 11)]
  fplot_eqtl$L95 <- fplot_eqtl$beta - 1.96 * fplot_eqtl$se
  fplot_eqtl$U95 <- fplot_eqtl$beta + 1.96 * fplot_eqtl$se
  fplot_eqtl$descript <- paste(fplot_eqtl$analysis, fplot_eqtl$entity, sep = '_')
  # reverses the factor level ordering for labels after coord_flip()
  fplot_eqtl$descript <- factor(fplot_eqtl$descript, levels = rev(fplot_eqtl$descript))

  date = date
  title <- paste('Simple Forest Plot of', myquery, 'eQTL Associations :', date, sep=" ")

  fp <- ggplot(data = fplot_eqtl, aes(y = descript, x = beta, xmin = L95, xmax = U95)) +
    geom_errorbarh(height = 0.3, col = 'darkgreen', alpha = 0.6) +
    geom_point(aes(size = round(-log10(pval), 2)), col = 'darkgreen') +
    geom_vline(xintercept = 0, color = 'grey50', linetype = 'dashed') +
    theme_bw() +
    xlab("beta estimate (95% CI)") + 
    ylab(" ") +
    labs(title = title)   
  
  suppressMessages(ggplotly(fp))
}
   

 



###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
#' @export
fplot_search <- function(res, term = 'asthma') {

  # Now paint a Forest Plot of the 10 most statistically significant disease associations.
  # Turn off warnings
  options(warn=-1)                 # To turn back on  ...  options(warn=0)
 
  idx <- grep(term, res$disease[ ,'description'], ignore.case = TRUE, value=F)
  # implement an error check here 
  if(length(idx) == 0) {
    print( paste("'", term, "'", " was not found in this result.  Please try an alternative term.", sep = ""))
    return(0)
  }
  search_disease <- res$disease[idx, ]
  search_disease <- search_disease[order(search_disease$pval), ]

  fplot_search_disease <- search_disease[ , c(1, 2, 3, 12, 10, 11)]
  fplot_search_disease$L95 <- fplot_search_disease$beta - 1.96 * fplot_search_disease$se
  fplot_search_disease$U95 <- fplot_search_disease$beta + 1.96 * fplot_search_disease$se
  date = date
  title <- paste('Simple Forest Plot of', myquery, term, 'Disease Associations', sep=" ")

  fplot_search_disease$description <- factor(fplot_search_disease$description, levels = fplot_search_disease$description)

  search_fp <- ggplot(data = fplot_search_disease, aes(y = description, x = beta, xmin = L95, xmax = U95)) +
    geom_errorbarh(height = 0.3, col = 'orangered', alpha = 0.5) +
    geom_point(aes(size = round(-log10(pval), 2)), col = 'orangered') +
    geom_vline(xintercept = 0, color = 'grey50', linetype = 'dashed', alpha = 0.5) +
    scale_y_discrete(labels = function(x) str_wrap(fplot_search_disease$description, width = 60)) +
    labs(title = title, x = "beta estimate (95% CI)", y = " ", caption = date) +
    scale_size(name = "-log10(pvalue)") +  
    theme(axis.title = element_text(size = 11, face = "bold"), axis.title.x = element_text(colour = "black", size = 10, face = "bold"), axis.text.y = element_text(size = 10), plot.caption = element_text(size = 9, face = "italic")) +
    theme_minimal()
 
  return(search_fp)
  
  #suppressMessages(ggplotly(search_fp))
}
   


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
