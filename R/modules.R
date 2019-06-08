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
    gtx_error("validate_module_input | Missing input.")
    stop();
  }
  gtx_debug("validate_module_input | Setup.")
  dots_quos <- dplyr::quos(...)
  
  validation_fxn <- function(x){
    x_name <- quo_name(x)
    if(any(stringr::str_detect(rlang::eval_tidy(x), 
                               "DOUBLE_CLICK_HERE_TO_CHANGE")) == TRUE){
      gtx_error("validate_module_input | \`{x_name}\` needs to be changed \\
                from the default \"DOUBLE_CLICK_HERE_TO_CHANGE\".")
      stop(call. = FALSE);
    }
  }
 
  gtx_debug("validate_module_input | Evaluate inputs.")
  purrr::map(dots_quos, validation_fxn)  
  
  gtx_debug("validate_module_input | Done.")
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
    gtx_error("config_db | The input file doesn't exist:{file}")
    stop()
  }
  
  # Safely read in the config file
  safely_read_csv <- purrr::safely(readr::read_csv)
  exec <- safely_read_csv(file = file,
                          col_types = cols("key"   = col_character(), 
                                           "value" = col_character()))
  
  if (!is.null(exec$error)){
    gtx_error("config_db | Unable to use the config file: {file} because: {exec$error}")
    stop()
  }
  
  mandatory_cols <- c("key", "value")
  if(any(purrr::map_lgl(mandatory_cols, ~ .x %in% names(exec$result)) == FALSE)){
    gtx_error("config_db | config file is missing \"key\" and/or \"value\" col header.")
    stop()
  }
  
  # Select the key-value pair where key=database
  ret <- exec$result %>% dplyr::filter(key == "database") %>% dplyr::select(value)
  
  if(nrow(ret) != 1){
    gtx_error("config_db | key \"database\" does not have exactly 1 row/value.")
  }
  ret <- dplyr::pull(ret)
  
  # Validate the database name doesn't have "DROP" in it for SQL protection
  if(stringr::str_detect(ret, stringr::coll("DROP", ignore_case = TRUE))){
    gtx_error("config_db | database value has \"drop\" in it, refusing to use it.")
    stop()
  }
  
  # Validate the database name has valid characters in it
  if(!stringr::str_detect(ret, stringr::regex("^[A-Za-z0-9_]+$"))){
    illegal <- str_replace_all(ret, stringr::regex("[A-Za-z0-9_]+"), "")
    if(!stringr::str_detect(illegal, stringr::regex(".+"))){
      gtx_error("config_db | database value contains no characters, value:{ret}")
      stop()   
    }
    else {
      gtx_error("config_db | database value contains illegal \\
                characters: {illegal} : value:{ret}")
      stop() 
    }
  }
  
  return(ret)
}

#' add_ssh_known_host
#'
#'  This function adds a url to ssh known hosts. This allows us to enable users 
#'  to run some code (a module) to automatically do git updates without the user 
#'  needing to manually uthenticate that they want to add the url as a known host (e.g. github.com)
#'  when you do a git query (fetch/pull/etc).
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @param known_host A host to search for and add as known host. e.g. "github.com"
#' @import futile.logger
#' @import glue
#' @import dplyr
#' @import stringr
add_ssh_known_host <- function(known_host){
  
  if(missing(known_host)){
    gtx_error("add_ssh_known_host | Must specify known_host = X.")
    stop()
  }
  
  home_dir <- Sys.getenv(x = "HOME", unset = NA)
  if(is.na(home_dir)){
    gtx_error("add_ssh_known_host | Cannot determine home directory.")
    stop()
  }
  
  known_hosts_file <- glue::glue("{home_dir}/.ssh/known_hosts")
  gtx_debug("add_ssh_known_host | ssh known_hosts file determined: {known_hosts_file}")  
  if(file.exists(known_hosts_file)){
    search_file <- system2(command = "grep", 
                           args = c(known_host, known_hosts_file), 
                           stdout = TRUE)
    
    if(!is.null(purrr::pluck(search_file, attr_getter("status")))){
      gtx_error("add_ssh_known_host | Unable to properly grep the known_hosts file.")  
    }
    
    if(!is_bare_character(search_file)){
      if(str_detect(string = search_file, pattern = known_host)){
        gtx_debug("add_ssh_known_host | Found previous url in known host, skipping add.")  
        return()
      }
    }
  }
  
  exec <- system2(command = "ssh-keyscan", 
                  args   = known_host, 
                  stdout = glue::glue("{home_dir}/tmp_keyscan"),
                  stderr = FALSE)
  
  if(!is.null(purrr::pluck(exec, attr_getter("status")))){
    gtx_error("add_ssh_known_host | Unable to ssh-keyscan.")
    stop()
  }
  
  exec <- system2(command = "ssh-keygen", 
                  args   = c("-lf", glue::glue("{home_dir}/tmp_keyscan")),
                  stdout = FALSE,
                  stderr = FALSE)
  
  if(!is.null(purrr::pluck(exec, attr_getter("status")))){
    gtx_error("add_ssh_known_host | Unable to ssh-keyscan.")  
    stop()
  }
  
  exec <- file.append(known_hosts_file, glue::glue("{home_dir}/tmp_keyscan"))
  if(exec == FALSE){
    gtx_error("add_ssh_known_host | Unable to append keyscan to known_hosts.")  
    stop()
  }
  gtx_debug("add_ssh_known_host | completed adding: {known_host}")
}

#' config_tmp_write_db
#' 
#' \strong{config_tmp_write_db - Read default project config database name to write temporary tables to.}
#' Read in a default config file and identify the named database to use to write temporary tables to. 
#' This is useful in different projects where you want the code to point to different databses.
#' The config file is csv, 2 cols (key, value), and config_tmp_write_db uses key=tmp_write_db
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @param file [Default = ~/.gtx_config.csv] Path to the config file to read in
#' @return String of the database name
#' @examples
#' example_tbl <- impala_copy_to(df = your_df, dest = validate_impala(), database = config_tmp_write_db())
#' @import futile.logger
#' @import glue
#' @import purrr
#' @import readr
config_tmp_write_db <- function(file = "~/.gtx_config.csv"){
  # Check the file exists
  if(!file.exists(file)){
    gtx_error("config_db | The input file doesn't exist:{file}")
    stop()
  }
  
  # Safely read in the config file
  safely_read_csv <- purrr::safely(readr::read_csv)
  exec <- safely_read_csv(file = file,
                          col_types = cols("key"   = col_character(), 
                                           "value" = col_character()))
  
  if (!is.null(exec$error)){
    gtx_error("config_db | Unable to use the config file: {file} because: {exec$error}")
    stop()
  }
  
  mandatory_cols <- c("key", "value")
  if(any(purrr::map_lgl(mandatory_cols, ~ .x %in% names(exec$result)) == FALSE)){
    gtx_error("config_db | config file is missing \"key\" and/or \"value\" col header.")
    stop()
  }
  
  # Select the key-value pair where key=tmp_write_db
  ret <- 
    exec$result %>% 
    dplyr::filter(key == "tmp_write_db") %>% 
    dplyr::select(value)
  
  if(nrow(ret) != 1){
    gtx_error("config_db | key \"tmp_write_db\" does not have exactly 1 row/value.")
  }
  ret <- dplyr::pull(ret)
  
  # Validate the database name doesn't have "DROP" in it for SQL protection
  if(stringr::str_detect(ret, stringr::coll("DROP", ignore_case = TRUE))){
    gtx_error("config_db | tmp_write_db value has \"drop\" in it, refusing to use it.")
    stop()
  }
  
  # Validate the database name has valid characters in it
  if(!stringr::str_detect(ret, stringr::regex("^[A-Za-z0-9_]+$"))){
    illegal <- str_replace_all(ret, stringr::regex("[A-Za-z0-9_]+"), "")
    if(!stringr::str_detect(illegal, stringr::regex(".+"))){
      gtx_error("config_db | tmp_write_db value contains no characters, value:{ret}")
      stop()   
    }
    else {
      gtx_error("config_db | tmp_write_db value contains illegal \\
                characters: {illegal} : value:{ret}")
      stop() 
    }
  }
  
  return(ret)
}