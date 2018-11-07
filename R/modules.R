#' Validate module input has been changed
#' 
#' \strong{validate_module_input - validate module input}
#' In the gwas_results modules, many inputs are set to "DOUBLE_CLICK_HERE_TO_CHANGE".
#' This functions validates that the input was changed from the default. 
#' It is important to wrap this functionality in a function to keep it hidden from users.
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @param .data A single/vector of inputs to validate. 
#' @return component \code{error} as TRUE or FALSE.
#' @examples 
#' inputs <- validate_module_input(hgncid)
#' inputs <- validate_module_input(c(chrom, pos))

#' @import stringr
#' @import futile.logger
validate_module_input <- function(.data = NULL){
  if(is.null(.data)){
    flog.error("validate_mInput | no input specified. Example: valid = validate_module_input(hgncid) ")
    stop;
  }
  flog.debug("validate_mInput | checking input strings")
  ret = list(error = any(str_detect(.data, "DOUBLE_CLICK_HERE_TO_CHANGE") == TRUE))
  
  flog.debug("validate_mInput | Done, returning results.")
  return(ret)
}