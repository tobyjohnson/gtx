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
#' @examples 
#' sc <- validate_sc(sc)
#' sc <- validate_sc()
#' @return Spark connection
#' @import purrr
#' @import sparklyr
#' @import futile.logger
#' @import glue
validate_sc <- function(sc = getOption("gtx.sc", NULL)){
  # Check we have a spark connection
  flog.debug("tidy_connections::validate_sc | Validating Spark connection.")
  if(is.null(sc) | !any(str_detect(class(sc), "spark_connection"))){ 
    flog.debug("tidy_connections::validate_sc | Spark connection is not valid, will try to establish a connection.")
    # Try to establish a spark connection 
    sc_config <- spark_config()
    sc_config$spark.port.maxRetries <- 60
    
    safe_spark_connect <- purrr::safely(spark_connect)
    safe_sc <- safe_spark_connect(master     = "yarn-client",
                                  spark_home = "/opt/cloudera/parcels/SPARK2/lib/spark2",
                                  version    = "2.1",
                                  config     = sc_config)
    # If there was no error, use the connection
    if(is.null(safe_sc$error)){ 
      sc <- safe_sc$result 
      # Record that we made a connection
      attr(sc, "internal_conn") <- TRUE
      flog.debug("tidy_connections::validate_sc | Spark connection established.")
    }
    # Otherwise advise user to manually create and pass the spark connection
    else if (!is.null(safe_sc$error)){
      flog.error(glue('tidy_connections::validate_sc | Unable to make the spark connection. Try initiating the spark connection manually and passing the connection to through option sc. To build a spark connection run:\n
                      \tsc <- spark_connect(master = {double_quote("yarn-client")}, \n\t\tspark_home = {double_quote("/opt/cloudera/parcels/SPARK2/lib/spark2")}, \n\t\tversion = {double_quote("2.1")})\n
                      After establishing the spark connection (above), you can also set the gtx options to use this connection by default running:\n
                      \t\toptions(gtx.sc = sc)'))
      stop()
    }
  }
  else if(any(str_detect(class(sc), "spark_connection"))){
    flog.debug("tidy_connections::validate_sc | Spark connection is valid.")
  }
  else{
    flog.error("tidy_connections::validate_sc | Unsure how to process sc.")
    stop()
  }
  return(sc)
  }

#' Validate and/or establish an impala odbc connection
#' 
#' \strong{validate_impala() - Validate impala connection}
#' This will validate an implyr impala connection. If it does not 
#' exists, it will try to safely establish the the connection 
#' and return it. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @param impala impala connection to validate. Defaults to checking options(gtx.impala).
#' @examples 
#' impala <- validate_impala()
#' @return implyr impala connection
#' @import purrr
#' @import implyr
#' @import odbc
#' @import futile.logger
validate_impala <- function(impala = getOption("gtx.impala", NULL)){
  # Check we have a spark connection
  flog.debug("tidy_connections::validate_impala | Validating impala connection.")
  if(is.null(impala) | !any(str_detect(class(impala), "src_impala"))){ 
    flog.debug("tidy_connections::validate_impala | impala connection is not valid, will try to establish a connection.")
    # Try to establish an impala connection 
    safe_connect <- purrr::safely(implyr::src_impala)
    safe_con     <- safe_connect(drv = odbc::odbc(), dsn = "impaladsn")
    
    # If there was no error, use the connection
    if(is.null(safe_con$error)){ 
      impala <- safe_con$result 
      # Record that we made a connection
      flog.debug("tidy_connections::validate_impala | impala connection established.")
      attr(impala, "internal_conn") <- TRUE
    }
    # Otherwise advise user to manually create and pass the impala connection
    else if (!is.null(safe_con$error)){
      flog.error('tidy_connections::validate_impala | Unable to make the impala connection.
                 FIRST, make sure you have used kinit.
                 Second, try to create+trouble shoot your own connection using: 
                   impala = implyr::src_impala(drv = odbc::odbc(), dsn = "impaladsn")')
      stop()
    }
  }
  else if(any(str_detect(class(impala), "src_impala"))){
    flog.debug("tidy_connections::validate_impala | impala connection is valid.")
  }
  else{
    flog.error("tidy_connections::validate_impala | Unsure how to process impala Please check input and read documentation.")
    stop()
  }
  return(impala)
}
############################################
#' Copy data to user's home database for joins with impala tables.
#' 
#' \strong{impala_copy_to() - copy data to RDIP}
#' This will copy a dataframe to the home directory of the user 
#' and return the 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @param df Data to copy to RDIP
#' @param dest Impala implyr connection
#' @return table reference to the \code{.data} in RDIP
#' @import dplyr
#' @import glue
#' @import implyr
#' @import readr
#' @import whoami
#' @import futile.logger
############################################
impala_copy_to <- function(df, dest = getOption("gtx.impala", NULL)){
  
  # Validate input has rows
  if(nrow(df) == 0){ flog.error("tidy_connections::impala_copy_to | The data has no rows.") }
  
  # Validate impala connection
  impala = validate_impala(impala = dest)
  
  # Determine user name for tables
  user_name <- whoami::username()
  if(is.null(user_name)){ 
    flog.error("tidy_connections::impala_copy_to | Unable to determine user name.") 
    return()
  }
  ############################################
  # If we have a small data frame, copy the data directly to RDIP
  if(nrow(df) < 250){
    safely_dbExecute <- purrr::safely(implyr::dbExecute)
    
    sql_statement <- glue("DROP TABLE IF EXISTS {`user_name`}.tmp_data4join PURGE")
    
    exec <- safely_dbExecute(impala, sql_statement)
    if (!is.null(exec$error)){
      flog.error(glue("tidy_connections::impala_copy_to | unable to remove {user_name}.tmp_data4join because:\n{exec$error}"))
      return()
    }
    
    input_tbl <- 
      copy_to(
        dest = impala, 
        df = df,
        name = glue("{`user_name`}.tmp_data4join"), 
        temporary = FALSE)
  }
  ############################################
  # If we have a big table, need to more complicated copy
  else {
    input_tbl <- big_copy_to(df = df, dest = impala)
  }
  return(input_tbl)
}

############################################
#' Copy large data to RDIP tables.
#' 
#' \strong{big_copy_to() - copy data to RDIP}
#' This will copy a dataframe to the home directory of the user 
#' and return the table reference to the data in RDIP. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @param df Data to copy to RDIP
#' @param dest Impala implyr connection
#' @param chrom_as_string Convert "chrom" col's to character instead of integers
#' @import glue
#' @import implyr
#' @import dplyr
#' @import futile.logger
#' @import whoami
#' @import purrr
############################################
big_copy_to <- function(df, dest = getOption("gtx.impala", NULL), chrom_as_string = TRUE){
  # Data 
  if(any(str_detect(names(df), "chrom")) & chrom_as_string == TRUE){
    if(typeof(df$chrom) != "character"){
      flog.warn(glue("column \"chrom\" detected, forcing it to be type character. 
                     Set chrom_as_string to FALSE if you don't want this."))
      df <-
        df %>% 
        mutate(chrom = as.character(chrom))
    }
  }
  # Determine user name for tables
  user_name <- whoami::username()
  if(is.null(user_name)){ 
    flog.error("Unable to determine user name which is needed for copying data to RDIP.") 
    return()
  }
  
  ############################################
  # Write the data 4 join to a tmp file on the edge node
  flog.debug("tidy_connections::big_copy_to | create tmp dir")
  system("mkdir -p ~/tmp")
  
  flog.debug("tidy_connections::big_copy_to | write data to tmp csv")
  write_csv(df, path = "~/tmp/tmp_data4join.csv")
  
  # Move the data to hdfs
  cmd <- glue("hdfs dfs -mkdir -p /user/{user_name}/staging/")
  flog.debug(glue("tidy_connections::big_copy_to | create staging dir on hdfs"))
  system(cmd)
  
  cmd <- glue("hdfs dfs -put -f ~/tmp/tmp_data4join.csv /user/{user_name}/staging/")
  flog.debug(glue("tidy_connections::big_copy_to | move data to hdfs"))
  system(cmd)
  
  cmd <- "rm ~/tmp/tmp_data4join.csv"
  flog.debug(glue("tidy_connections::big_copy_to | remove tmp file on edge node"))
  system(cmd)
  ############################################
  # Make a list for conversion of R class to SQL class
  sql_types <- tibble(r_class   = c("character", "integer", "double", "logical"),
                      sql_class = c("STRING",    "INT",     "DOUBLE", "BOOLEAN"))
  
  # Determine proper SQL col types
  col_info <- tibble(col_names = names(df), col_class = map_chr(df, typeof))
  col_info <- 
    left_join(col_info, 
              sql_types, 
              by = c("col_class" = "r_class"))
  # mutate(sql_create_tbl = glue("{col_names} {sql_class}"))
  if(any(is.na(col_info$sql_class))){
    bad_cols <- 
      col_info %>% 
      filter(is.na(sql_class)) %>% 
      pull(col_names) %>% 
      glue_collapse(., sep = ",")
    
    flog.error(glue("Could not match {bad_cols}"))
    return()
  }
  ############################################
  # Drop pre-existing tmp table
  safely_dbExecute <- purrr::safely(implyr::dbExecute)
  
  sql_statement <- glue("DROP TABLE IF EXISTS {`user_name`}.tmp_data4join PURGE")
  
  exec <- safely_dbExecute(dest, sql_statement)
  if (!is.null(exec$error)){
    flog.error(glue("tidy_connections::big_copy_to | unable to remove {user_name}.tmp_data4join because:\n{exec$error}"))
    return()
  }
  ############################################
  # Create new table based on correct col types
  sql_statement <- 
    glue(
      glue("CREATE TABLE {`user_name`}.tmp_data4join ("), 
      glue_collapse(glue("`{col_info$col_names}` {col_info$sql_class}"), sep = ", ", last = ""),
      glue(") 
           ROW FORMAT DELIMITED FIELDS TERMINATED BY \",\"
           STORED AS TEXTFILE
           TBLPROPERTIES(\"skip.header.line.count\"=\"1\")"))
  
  exec <- safely_dbExecute(dest, sql_statement)
  if (!is.null(exec$error)){
    flog.error(glue("tidy_connections::big_copy_to | unable to create table: {user_name}.tmp_data4join because:\n{exec$error}"))
    return()
  }
  ############################################
  # Load data from tmp HDFS file into table
  sql_statement <- 
    glue("LOAD DATA INPATH '/user/{`user_name`}/staging/tmp_data4join.csv' 
         INTO TABLE {`user_name`}.`tmp_data4join`")
  
  exec <- safely_dbExecute(dest, sql_statement)
  if (!is.null(exec$error)){
    flog.error(glue("tidy_connections::big_copy_to | Unable to load tmp_data4join.csv into table: {user_name}.tmp_data4join because:\n{exec$error}"))
    return()
  }
  ############################################
  # Return table to data in RDIP
  input_tbl <- tbl(dest, glue("{user_name}.tmp_data4join"))
  return(input_tbl)
}
#############################################
#' Close connections they were internally created
#' 
#' \strong{close_int_conn() - Close internally created connections}
#' This will close implyr or sparklyr connections that were internally created.
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @param conn Connection to check
#' @examples 
#' close_int_conn(conn = impala)
#' @import sparklyr
#' @import implyr
#' @import odbc
#' @import futile.logger
close_int_conn <- function(conn){
  internal_conn_atrr <- pluck(conn, attr_getter("internal_conn"))
  if(!is.null(internal_conn_atrr)){
    if(internal_conn_atrr ==  TRUE){
      flog.debug("tidy_connections::close_int_conn | internal connection detected.")
      # Check for impala connection
      if(any(str_detect(class(conn), "src_impala"))){
        implyr::dbDisconnect(conn)
        flog.debug("tidy_connections::close_int_conn | impala connection closed.")
      }
      # Check for sparklyr connection
      else if(any(str_detect(class(sc), "spark_connection"))){
        sparklyr::spark_disconnect(conn)
        flog.debug("tidy_connections::close_int_conn | spark connection closed.")
      }
    }
  }
  else {
    flog.debug("tidy_connections::close_int_conn | internal connection not detected.")
  }
}