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
    flog.error("validate_module_input | Missing input.")
    stop();
  }
  flog.debug("validate_module_input | Setup.")
  dots_quos <- quos(...)
  
  validation_fxn <- function(x){
    x_name <- quo_name(x)
    if(any(str_detect(rlang::eval_tidy(x), "DOUBLE_CLICK_HERE_TO_CHANGE")) == TRUE){
      flog.error(glue("validate_module_input | \`{x_name}\` needs to be changed from the default \"DOUBLE_CLICK_HERE_TO_CHANGE\"."))
      stop(call. = FALSE);
    }
  }
 
  flog.debug("validate_module_input | Evaluate inputs.")
  map(dots_quos, validation_fxn)  
  
  flog.debug("validate_module_input | Done.")
}
