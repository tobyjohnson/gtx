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

#' gtx_fatal_stop
#' 
#' Wrapper around futile.logger::flog.fatal(glue::glue(.msg)), which
#' also throws a true fatal error with the same message, using 
#' stop().  Useful if calling within a function body where we want to 
#' stop execution, and also ensures the stop() message goes to
#' the console (even if futile.logger messages have been redirected).
#' 
#' @author Toby Johnson \email{Toby.x.Johnson@gsk.com}
#' @export
#' @family gtx_log
#' @param .msg msg to glue and flog
#' @import glue
#' @import futile.logger
gtx_fatal_stop <- function(.msg){
  # evaluate once, just in case there are side effects:
  msg_eval <- glue::glue(.msg, .envir = parent.frame())
  futile.logger::flog.fatal(msg_eval)
  stop(msg_eval, call. = FALSE)
  # changing call. = TRUE would just make the stop message
  # say that the error occurred in gtx_fatal2().  We want
  # to inspect the parent.frame() to find where the error
  # *actually* occurred...  FIXME
}
