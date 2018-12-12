#' Validate module input has been changed
#' 
#' \strong{validate_module_input - validate module input}
#' In the gwas_results modules, many inputs are set to "DOUBLE_CLICK_HERE_TO_CHANGE".
#' This functions validates that the input was changed from the default. 
#' It is important to wrap this functionality in a function to keep it hidden from users.
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @param ... A single or multiple comma seperated inputs to validate. Inputs can be strings or vectors. 
#' @return Returns nothing. Will fail if encounters a problem. 
#' @examples
#' validate_module_input(hgncid)
#' validate_module_input(chrom, pos_start, pos_end)
#' @import stringr
#' @import glue
#' @import futile.logger
#' @importFrom rlang eval_tidy
validate_module_input <- function(...){
  if(missing(...)){
    futile.logger::flog.error("validate_module_input | Missing input.")
    stop();
  }
  futile.logger::flog.debug("validate_module_input | Setup.")
  dots_quos <- dplyr::quos(...)
  
  validation_fxn <- function(x){
    x_name <- quo_name(x)
    if(any(stringr::str_detect(rlang::eval_tidy(x), "DOUBLE_CLICK_HERE_TO_CHANGE")) == TRUE){
      futile.logger::flog.error(glue::glue("validate_module_input | \`{x_name}\` needs to be changed from the default \"DOUBLE_CLICK_HERE_TO_CHANGE\"."))
      stop(call. = FALSE);
    }
  }
 
  futile.logger::flog.debug("validate_module_input | Evaluate inputs.")
  purrr::map(dots_quos, validation_fxn)  
  
  futile.logger::flog.debug("validate_module_input | Done.")
}

#' config_db
#' 
#' \strong{config_db - Read default project config database name}
#' Read in a default config file and identify the named database to use. 
#' This is useful in different projects where you want the code to point to different databses.
#' The config file is csv, 2 cols (key, value), and config_db uses key=database
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @param file [Default = ~/.gtx_config.csv] Path to the config file to read in
#' @return String of the database name
#' @examples
#' gtxconnection(use_database = config_db())
#' gtxconnection(use_database = config_db(file = "/some/path/config.csv"))
#' @import futile.logger
#' @import glue
#' @import purrr
#' @import readr

config_db <- function(file = "~/.gtx_config.csv"){
  # Check the file exists
  if(!file.exists(file)){
    futile.logger::flog.error(glue::glue("config_db | The input file doesn't exist:{file}"))
    stop()
  }
  
  # Safely read in the config file
  safely_read_csv <- purrr::safely(readr::read_csv)
  exec <- safely_read_csv(file = file,
                          col_types = cols("key"   = col_character(), 
                                           "value" = col_character()))
  
  if (!is.null(exec$error)){
    futile.logger::flog.error(glue::glue("config_db | Unable to use the config file: {file} because: {exec$error}"))
    stop()
  }
  
  mandatory_cols <- c("key", "value")
  if(any(purrr::map_lgl(mandatory_cols, ~ .x %in% names(exec$result)) == FALSE)){
    futile.logger::flog.error("config_db | config file is missing \"key\" and/or \"value\" col header.")
    stop()
  }
  
  # Select the key-value pair where key=database
  ret <- exec$result %>% dplyr::filter(key == "database") %>% dplyr::select(value)
  
  if(nrow(ret) != 1){
    futile.logger::flog.error("config_db | key \"database\" does not have exactly 1 row/value.")
  }
  ret <- dplyr::pull(ret)
  
  # Validate the database name doesn't have "DROP" in it for SQL protection
  if(stringr::str_detect(ret, stringr::coll("DROP", ignore_case = TRUE))){
    futile.logger::flog.error("config_db | database value has \"drop\" in it, refusing to use it.")
    stop()
  }
  
  # Validate the database name has valid characters in it
  if(!stringr::str_detect(ret, stringr::regex("^[A-Za-z0-9_]+$"))){
    illegal <- str_replace_all(ret, stringr::regex("[A-Za-z0-9_]+"), "")
    if(!stringr::str_detect(illegal, stringr::regex(".+"))){
      futile.logger::flog.error(glue::glue("config_db | database value contains no characters, value:{ret}"))
      stop()   
    }
    else {
      futile.logger::flog.error(glue::glue("config_db | database value contains illegal characters: {illegal} : value:{ret}"))
      stop() 
    }
  }
  
  return(ret)
}