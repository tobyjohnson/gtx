#' gtx_error
#' 
#' Simple wrapper around futile.logger::flog.error(glue::glue(.msg)).
#' For full control over flog or glue opts use full explicit code. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @family gtx_log
#' @param .msg msg to glue and flog
#' @import glue
#' @import futile.logger
gtx_error <- function(.msg){
  futile.logger::flog.error(glue::glue(.msg, .envir = parent.frame()))
}

#' gtx_debug
#' 
#' Simple wrapper around futile.logger::flog.debug(glue::glue(.msg)).
#' For full control over flog or glue opts use full explicit code. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @family gtx_log
#' @param .msg msg to glue and flog
#' @import glue
#' @import futile.logger
gtx_debug <- function(.msg){
  futile.logger::flog.debug(glue::glue(.msg, .envir = parent.frame()))
}

#' gtx_info
#' 
#' Simple wrapper around futile.logger::flog.info(glue::glue(.msg)).
#' For full control over flog or glue opts use full explicit code. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @family gtx_log
#' @param .msg msg to glue and flog
#' @import glue
#' @import futile.logger
gtx_info <- function(.msg){
  futile.logger::flog.info(glue::glue(.msg, .envir = parent.frame()))
}

#' gtx_warn
#' 
#' Simple wrapper around futile.logger::flog.warn(glue::glue(.msg)).
#' For full control over flog or glue opts use full explicit code. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @family gtx_log
#' @param .msg msg to glue and flog
#' @import glue
#' @import futile.logger
gtx_warn <- function(.msg){
  futile.logger::flog.warn(glue::glue(.msg, .envir = parent.frame()))
}

#' gtx_fatal
#' 
#' Simple wrapper around futile.logger::flog.fatal(glue::glue(.msg)).
#' For full control over flog or glue opts use full explicit code. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @family gtx_log
#' @param .msg msg to glue and flog
#' @import glue
#' @import futile.logger
gtx_fatal <- function(.msg){
  futile.logger::flog.fatal(glue::glue(.msg, .envir = parent.frame()))
}