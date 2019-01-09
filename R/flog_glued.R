#' flog_error
#' 
#' Simple wrapper around futile.logger::flog.error(glue::glue(.msg)).
#' For full control over flog or glue opts use full explicit code. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @param .msg
#' @import glue
#' @import futile.logger
flog_error <- function(.msg){
  futile.logger::flog.error(glue::glue(.msg, .envir = parent.frame()))
}

#' flog_debug
#' 
#' Simple wrapper around futile.logger::flog.debug(glue::glue(.msg)).
#' For full control over flog or glue opts use full explicit code. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @param .msg
#' @import glue
#' @import futile.logger
flog_debug <- function(.msg){
  futile.logger::flog.debug(glue::glue(.msg, .envir = parent.frame()))
}

#' flog_info
#' 
#' Simple wrapper around futile.logger::flog.info(glue::glue(.msg)).
#' For full control over flog or glue opts use full explicit code. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @param .msg
#' @import glue
#' @import futile.logger
flog_info <- function(.msg){
  futile.logger::flog.info(glue::glue(.msg, .envir = parent.frame()))
}

#' flog_warn
#' 
#' Simple wrapper around futile.logger::flog.warn(glue::glue(.msg)).
#' For full control over flog or glue opts use full explicit code. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @param .msg
#' @import glue
#' @import futile.logger
flog_warn <- function(.msg){
  futile.logger::flog.warn(glue::glue(.msg, .envir = parent.frame()))
}