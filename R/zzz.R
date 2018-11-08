#' .onAttach - Report GTX package info
#' 
#' To time & version stamp notebooks, this internal function will parse and report the github SHA tag.
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @import stringr
#' @import glue
#' @import futile.logger
#' @import dplyr
#' @import purrr
.onAttach <- function(libname, pkgname){
  description <- packageDescription("gtx")
  
  params <- c("GithubRepo", "GithubRef", "GithubSHA1", "Packaged")
  
  if(!all(params %in% names(description))){
    flog.warn("GTX was not installed from github, no SHA version available.")
  }
  else {
    version_info <-
      params %>% 
      map_chr(~pluck(description, .)) %>% 
      glue_collapse(., sep = " | ")
    
    flog.info(glue("GTX package info: {version_info}"))
  }
}

# gtx_dirs <- 
#   .libPaths() %>% 
#   map(~dir(path = ., pattern = "gtx", full.names = TRUE)) %>% 
#   compact() 
# 
# if(length(gtx_dirs) > 1){
#   flog.warn(glue("Multiple gtx installations detected, defaulting to: {pluck(gtx_dirs, 1)}"))
# }
# 
# default_description_file <- 
#   pluck(gtx_dirs, 1) %>% 
#   list.files(path = ., pattern = "DESCRIPTION", full.names = TRUE)
# 
# description <- suppressMessages(suppressWarnings(
#   read_delim(file = default_description_file, 
#              delim = ":", 
#              col_names = c("key", "value"), 
#              trim_ws = TRUE)))
# 
# if(description %>% filter(key == "Package") %>% pluck("value", 1) != "gtx"){
#   flog.warn(glue("The DESCRIPTION file ({default_description_file}) \"Package\" doesn't match \"gtx\". "))
# }
# 
# if(description %>% filter(str_detect(key, "Github")) %>% nrow == 0){
#   flog.warn(glue("GTX package was not instaled from github, unable to determine version/SHA."))
# }
# 
# version_info <- 
#   description %>% 
#   filter(str_detect(key, "Github") | key == "Packaged") %>% 
#   filter(key != "GithubUsername") %>% 
#   pluck("value") %>% 
#   glue_collapse(sep = " | ")