#' Validate and/or establish a spark connection
#' 
#' \strong{validate_sc() - Validate spark connection}
#' This will validate a spark connection. If it does not 
#' exists, it will try to safely establish the the connection 
#' and return it. Note, if \code{\link{validate_sc}} has to 
#' create the Spark connection it will use the default settings. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @param sc Spark connection to validate. Defaults to checking options(gtx.sc).
#' @param flog_level [Default = getOption("gtx.flog_level", "WARN")] \code{\link{futile.logger}} \code{\link{flog.threshold}} [INFO, WARN, ERROR]
#' @examples 
#' sc <- validate_sc(sc)
#' sc <- validate_sc()
#' @return Spark connection
#' @import sparklyr
#' @import futile.logger
#' @import glue
validate_sc <- function(sc         = getOption("gtx.sc", NULL),
                        flog_level = getOption("gtx.flog_level", "WARN")){
  if(flog_level != "INFO" & flog_level != "WARN" & flog_level != "ERROR"){
    flog.error(glue("flog_level defined as: {flog_level}. Needs to be [INFO, WARN, or ERROR]."))
  }
  flog.threshold(flog_level)

  # Check we have a spark connection
  flog.info("Validating Spark connection.")
  if(is.null(sc) | !any(str_detect(class(sc), "spark_connection"))){ 
    flog.info("Spark connection is not valid, will try to establish a connection.")
    # Try to establish a spark connection 
    sc_config <- spark_config()
    sc_config$spark.port.maxRetries <- 60
    
    safe_spark_connect <- safely(spark_connect)
    safe_sc <- safe_spark_connect(master     = "yarn-client",
                                  spark_home = "/opt/cloudera/parcels/SPARK2/lib/spark2",
                                  version    = "2.1",
                                  config     = sc_config)
    # If there was no error, use the connection
    if(is.null(safe_sc$error)){ 
      sc <- safe_sc$result 
      # Record that we made a connection in aba.query_locus
      abaQ_sc_created = TRUE
      flog.info("Spark connection established.")
    }
    # Otherwise advise user to manually create and pass the spark connection
    else if (!is.null(safe_sc$error)){
      flog.fatal(glue('Unable to make the spark connection. Try initiating the spark connection manually and passing the connection to through option sc. To build a spark connection run:\n
                      \tsc <- spark_connect(master = {double_quote("yarn-client")}, \n\t\tspark_home = {double_quote("/opt/cloudera/parcels/SPARK2/lib/spark2")}, \n\t\tversion = {double_quote("2.1")})\n
                      After establishing the spark connection (above), you can also set the gtx options to use this connection by default running:\n
                      \t\toptions(gtx.sc = sc)'))
      stop()
    }
  }
  else if(any(str_detect(class(sc), "spark_connection"))){
    flog.info("Spark connection is valid.")
  }
  else{
    flog.fatal("Unsure how to process sc.")
    stop()
  }
  return(sc)
}