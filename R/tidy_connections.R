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

validate_sc <- function(sc = getOption("gtx.sc", NULL), spark_version = "2.2.0", app_name = "sparklyr", sc_config=spark_config()){

  # Check we have a spark connection
  futile.logger::flog.debug("tidy_connections::validate_sc | Validating Spark connection.")
  if(is.null(sc) | !any(stringr::str_detect(class(sc), "spark_connection"))){ 
    futile.logger::flog.debug("tidy_connections::validate_sc | Spark connection is not valid, will try to establish a connection.")
    # Try to establish a spark connection 
    sc_config <- sparklyr::spark_config()
    sc_config$spark.driver.memory                <- '16G'
    sc_config$spark.executor.memory              <- '16G'
    sc_config$spark.yarn.executor.memoryOverhead <- '8G'
    sc_config$spark.port.maxRetries              <- 60
    sc_config$spark.rpc.message.maxSize          <- 512    # This works best for uploading data via copy_to
    
    safe_spark_connect <- purrr::safely(sparklyr::spark_connect)
    safe_sc <- safe_spark_connect(master     = "yarn-client",
                                  spark_home = "/opt/cloudera/parcels/SPARK2/lib/spark2",
                                  version    = spark_version,
                                  config     = sc_config)
    # If there was no error, use the connection
    if(is.null(safe_sc$error)){ 
      sc <- safe_sc$result 
      # Record that we made a connection
      attr(sc, "internal_conn") <- TRUE
      futile.logger::flog.debug("tidy_connections::validate_sc | Spark connection established.")
    }
    # Otherwise advise user to manually create and pass the spark connection
    else if (!is.null(safe_sc$error)){
      futile.logger::flog.error(glue::glue('tidy_connections::validate_sc | Unable to make the spark connection. Try initiating the spark connection manually and passing the connection to through option sc. To build a spark connection run:\n
                      \tsc <- spark_connect(master = {double_quote("yarn-client")}, \n\t\tspark_home = {double_quote("/opt/cloudera/parcels/SPARK2/lib/spark2")}, \n\t\tversion = {double_quote("2.1")})\n
                      After establishing the spark connection (above), you can also set the gtx options to use this connection by default running:\n
                      \t\toptions(gtx.sc = sc)'))
      stop()
    }
  }
  else if(any(stringr::str_detect(class(sc), "spark_connection"))){
    futile.logger::flog.debug("tidy_connections::validate_sc | Spark connection is valid.")
    # Check if spark version and app_name matches those requested. Otherwise, disconnect and call function again
    if (grepl(spark_version, sc$spark_home)){
	    futile.logger::flog.debug("Spark version detected is different from current version. Restarting Spark connection")
	    sparklyr::spark_disconnect(sc)
	    sc = validate_sc(spark_version=spark_version, app_name=app_name)
    }
    if (sc$app_name != app_name){
	    futile.logger::flog.debug("App name requested is different from current value. Restarting Spark connection")
      sparklyr::spark_disconnect(sc)
	    sc = validate_sc(spark_version=spark_version, app_name=app_name)
    }
  }
  else{
    futile.logger::flog.error("tidy_connections::validate_sc | Unsure how to process sc.")
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
  # Check we have an implyr connection
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
      futile.logger::flog.debug("tidy_connections::validate_impala | impala connection established.")
      attr(impala, "internal_conn") <- TRUE
    }
    # Otherwise advise user to manually create and pass the impala connection
    else if (!is.null(safe_con$error)){
      futile.logger::flog.error('tidy_connections::validate_impala | Unable to make the impala connection.
                 FIRST, make sure you have used kinit.
                 Second, try to create+trouble shoot your own connection using: 
                   impala = implyr::src_impala(drv = odbc::odbc(), dsn = "impaladsn")')
      stop()
    }
  }
  else if(any(stringr::str_detect(class(impala), "src_impala"))){
    futile.logger::flog.debug("tidy_connections::validate_impala | impala connection is valid.")
  }
  else{
    futile.logger::flog.error("tidy_connections::validate_impala | Unsure how to process impala Please check input and read documentation.")
    stop()
  }
  return(impala)
}

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
#' @param database Defaults to user's database: \code{whoami}.
#' @param table_name [Default = "tmp_data4join"] Specify table name, overrides random_name = TRUE
#' @param random_name [Default = FALSE] TRUE = ~random number name for the table. 
#' @return table reference to the \code{df} in RDIP
#' @import dplyr
#' @import glue
#' @import implyr
#' @import readr
#' @import futile.logger
impala_copy_to <- function(df, dest = getOption("gtx.impala", NULL), 
                           database = NULL, table_name = "tmp_data4join", random_name = FALSE ){
  safely_dbExecute  <- purrr::safely(implyr::dbExecute)
  safely_dbGetQuery <- purrr::safely(implyr::dbGetQuery)
  # Validate input has rows
  if(nrow(df) == 0){ futile.logger::flog.error("tidy_connections::impala_copy_to | The data has no rows.") }
  
  # Validate impala connection
  impala = validate_impala(impala = dest)
  
  # Determine database for tables
  if(is.null(database)){
    database <- whoami()
    if(is.null(database)){ 
      futile.logger::flog.error("tidy_connections::impala_copy_to | \\
                                Unable to determine user name for database.") 
      stop()
    }
  }

  # If old tmp table, clean it up. 
  if(database == whoami() & table_name == "tmp_data4join" & random_name == FALSE){
    sql_statement <- glue::glue('SHOW TABLES IN {`database`} LIKE "{table_name}"')
    gtx_debug("tidy_connections::impala_copy_to | {sql_statement}")
    exec <- safely_dbGetQuery(impala, sql_statement)
    
    if (is.null(exec$error) & nrow(exec$result) >= 1){
      futile.logger::flog.debug(glue::glue("tidy_connections::impala_copy_to | \\
                                           dropping {database}.{table_name}"))
    }
    
    sql_statement <- glue::glue("INVALIDATE METADATA {database}.tmp_data4join")
    gtx_debug("tidy_connections::impala_copy_to | {sql_statement}")
    exec <- safely_dbExecute(impala, sql_statement)
    
    if (!is.null(exec$error)){
      futile.logger::flog.error(glue::glue("tidy_connections::impala_copy_to | \\
                                           unable to remove {database}.tmp_data4join because:\n{exec$error}"))
      stop()
    }
    
    sql_statement <- glue::glue("DROP TABLE IF EXISTS {database}.tmp_data4join PURGE")
    gtx_debug("tidy_connections::impala_copy_to | {sql_statement}")
    exec <- safely_dbExecute(impala, sql_statement)
    
    if (!is.null(exec$error)){
      futile.logger::flog.error(glue::glue("tidy_connections::impala_copy_to | \\
                                           unable to remove {database}.tmp_data4join because:\n{exec$error}"))
      stop()
    }
  }
  # If a new table, determine if the table already exists
  else if(table_name != "tmp_data4join" & random_name == FALSE){
    sql_statement <- glue::glue('SHOW TABLES IN {`database`} LIKE "{table_name}"')
    gtx_debug("tidy_connections::impala_copy_to | {sql_statement}")
    exec <- safely_dbGetQuery(impala, sql_statement)
    
    if (!is.null(exec$error)){ 
      gtx_error("tidy_connections::impala_copy_to | \\
                 failed to determine if {database}.{table_name} exist because:\n{exec$error}")
      stop()
    }
    # If we have rows here, we found tables that have conflicting names
    else if (nrow(exec$result) >=1){
      gtx_error("tidy_connections::impala_copy_to | \\
                {database}.{table_name} already exists.")
      stop();
    }
  }
  else if(random_name == TRUE){
    time_stamp <- as.integer(Sys.time())
    table_name <- glue::glue("tmp_{time_stamp}")
    
    sql_statement <- glue::glue('SHOW TABLES IN {`database`} LIKE "{table_name}"')
    gtx_debug("tidy_connections::impala_copy_to | {sql_statement}")
    exec <- safely_dbGetQuery(impala, sql_statement)
    
    if (!is.null(exec$error)){ 
      gtx_error("tidy_connections::impala_copy_to | \\
                 unable to query if random_name table exists because:\n{exec$error}")
      stop();
    }
    # If we have rows here, we found tables that have conflicting names
    else if (nrow(exec$result) >=1){
      gtx_error("tidy_connections::impala_copy_to | {database}.{table_name} already exists.")
      stop();
    }
  }
  else {
    gtx_error("tidy_connections::impala_copy_to | Unable to determine appropriate table_name.")
    stop();
  }

  # If we have a small data frame, copy the data directly to RDIP
  if(nrow(df) < 50 & ncol(df) < 10){
    input_tbl <- 
      dplyr::copy_to(
        dest = impala, 
        df   = df,
        name = glue::glue("{`database`}.{`table_name`}"), 
        temporary = FALSE)
  }

  # If we have a big table, need to more complicated copy
  else {
    input_tbl <- big_copy_to(df = df, dest = impala, database = database, table_name = table_name)
  }
  
  # Attach "tmp" status to delete
  attr(input_tbl, "tmp_impala_tbl") <- TRUE
  
  # Return final table
  return(input_tbl)
}

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
#' @param database Name of the database
#' @param table_name Name of the table within the database. 
#' @param chrom_as_string [Default = TRUE] TRUE = Override "chrom" cols as character instead of integers
#' @param compute_stats [Default = FALSE] TRUE = SQL execute COMPUTE STATS
#' @import readr
#' @import glue
#' @import implyr
#' @import dplyr
#' @import futile.logger
#' @import purrr
big_copy_to <- function(df, dest = getOption("gtx.impala", NULL), 
                        chrom_as_string = TRUE, database = NULL, 
                        table_name = NULL, compute_stats = FALSE){
  safely_dbExecute  <- purrr::safely(implyr::dbExecute)
  safely_dbGetQuery <- purrr::safely(implyr::dbGetQuery)
  safely_system     <- purrr::safely(system)
  if(is.null(database)){
    gtx_error("tidy_connections::big_copy_to | database is NULL.")
  }
  if(is.null(table_name)){
    gtx_error("tidy_connections::big_copy_to | table_name is NULL.")
  }
  
  # Data 
  if(any(stringr::str_detect(names(df), "chrom")) & chrom_as_string == TRUE){
    if(typeof(df$chrom) != "character"){
      futile.logger::flog.warn(glue::glue("column \"chrom\" detected, forcing it to be type character. 
                     Set chrom_as_string to FALSE if you don't want this."))
      df <-
        df %>% 
        dplyr::mutate(chrom = as.character(chrom))
    }
  }
  # Determine user name for tables
  user_name <- whoami()
  if(is.null(user_name)){ 
    futile.logger::flog.error("Unable to determine user name which is needed for copying data to RDIP.") 
    stop();
  }
  
  # --
  # Write the data 4 join to a tmp file on the edge node
  gtx_debug("tidy_connections::big_copy_to | create tmp dir")
  exec <- safely_system("mkdir -p ~/tmp")
  if(!is.null(exec$error)){
    gtx_error("tidy_connections::big_copy_to | Unable to create ~/tmp because: {exec$error}")
    stop();
  }
  # Export data to tmp.csv in user's home directory
  gtx_debug("tidy_connections::big_copy_to | write data to tmp csv")
  safely_write_csv <- purrr::safely(readr::write_csv)
  exec <- safely_write_csv(df, path = glue::glue("~/tmp/{table_name}.csv"))
  if(!is.null(exec$error)){
    gtx_error("tidy_connections::big_copy_to | Unable to export data to: ~/tmp/{table_name}.csv because: {exec$error}")
    stop();
  }
  
  # Make sure there is a user staging direction on HDFS
  gtx_debug("tidy_connections::big_copy_to | create staging dir on hdfs")
  exec <- safely_system(glue::glue("hdfs dfs -mkdir -p /user/{user_name}/staging/"))
  if(!is.null(exec$error)){
    gtx_error("tidy_connections::big_copy_to | \\
              Unable to create: /user/{user_name}/staging/ on HDFS because: {exec$error}")
    stop();
  }
  
  # Move data to HDFS
  gtx_debug("tidy_connections::big_copy_to | move data to hdfs")
  exec <- safely_system(glue::glue("hdfs dfs -put -f ~/tmp/{table_name}.csv /user/{user_name}/staging/"))
  if(!is.null(exec$error)){
    gtx_error("tidy_connections::big_copy_to | Unable to put data on HDFS because: {exec$error}")
    stop();
  }
  
  # Remove tmp data in user's home directory
  gtx_debug("tidy_connections::big_copy_to | remove tmp file on edge node")
  exec <- safely_system(glue::glue("rm ~/tmp/{table_name}.csv"))
  if(!is.null(exec$error)){
    gtx_error("tidy_connections::big_copy_to | Unable to remove tmp data in ~/tmp because: {exec$error}")
    stop();
  }
  # Data is now ready to be read into a table
  
  # Make a list for conversion of R class to SQL class
  sql_types <- 
    dplyr::tibble(r_class   = c("character", "integer", "double", "logical"),
                  sql_class = c("STRING",    "INT",     "DOUBLE", "BOOLEAN"))
  
  # Determine proper SQL col types
  col_info <- dplyr::tibble(col_names = names(df), col_class = purrr::map_chr(df, typeof))
  col_info <- 
    dplyr::left_join(col_info, 
                     sql_types, 
                     by = c("col_class" = "r_class"))

  if(any(is.na(col_info$sql_class))){
    bad_cols <- 
      col_info %>% 
      dplyr::filter(is.na(sql_class)) %>% 
      dplyr::pull(col_names) %>% 
      glue::glue_collapse(., sep = ",")
    
    futile.logger::flog.error(glue::glue("Could not match {bad_cols}"))
    stop()
  }

  # Create new table based on correct col types
  gtx_debug("tidy_connections::big_copy_to | create new table: {database}.{table_name}")
  sql_statement <- 
    glue::glue(
      glue::glue("CREATE TABLE {`database`}.{`table_name`} ("), 
      glue_collapse(glue::glue("`{col_info$col_names}` {col_info$sql_class}"), sep = ", ", last = ""),
      glue::glue(") \\
           ROW FORMAT DELIMITED FIELDS TERMINATED BY \",\" \\
           STORED AS TEXTFILE \\
           TBLPROPERTIES(\"skip.header.line.count\"=\"1\")"))
  
  gtx_debug("tidy_connections::impala_copy_to | {sql_statement}")
  exec <- safely_dbExecute(dest, sql_statement)
  if (!is.null(exec$error)){
    gtx_error("tidy_connections::big_copy_to | unable to create table: {database}.{table_name} because:\n{exec$error}")
    stop()
  }

  # Load data from tmp HDFS file into table
  gtx_debug("tidy_connections::big_copy_to | Load data into new table: {database}.{table_name}")
  sql_statement <- 
    glue::glue("LOAD DATA INPATH '/user/{`user_name`}/staging/{table_name}.csv' \\
                INTO TABLE {`database`}.{`table_name`}")
  gtx_debug("tidy_connections::impala_copy_to | {sql_statement}")
  exec <- safely_dbExecute(dest, sql_statement)
  if (!is.null(exec$error)){
    futile.logger::flog.error(glue::glue("tidy_connections::big_copy_to | Unable to load tmp_data4join.csv into table: {user_name}.{table_name} because:\n{exec$error}"))
    stop()
  }

  # COMPUTE STATS on the new table
  if(isTRUE(compute_stats)){
    gtx_debug("tidy_connections::big_copy_to | COMPUTE STATS on: {database}.{table_name}")
    sql_statement <- glue("COMPUTE STATS {`database`}.{`table_name`}")
    gtx_debug("tidy_connections::impala_copy_to | {sql_statement}")
    exec <- safely_dbExecute(dest, sql_statement)
    if (!is.null(exec$error)){
      gtx_error("tidy_connections::big_copy_to | Unable to COMPUTE STATS on: \\
                {database}.{table_name} because:\n{exec$error}")
      stop()
    }
  }
  
  # Return table to data in RDIP
  input_tbl <- dplyr::tbl(dest, glue::glue("{database}.{table_name}"))
  return(input_tbl)
}

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
  internal_conn_atrr <- purrr::pluck(conn, attr_getter("internal_conn"))
  
  if(!is.null(internal_conn_atrr)){
    if(internal_conn_atrr ==  TRUE){
      futile.logger::flog.debug("tidy_connections::close_int_conn | internal connection detected.")
      # Check for impala connection
      if(any(stringr::str_detect(class(conn), "src_impala"))){
        implyr::dbDisconnect(conn)
        futile.logger::flog.debug("tidy_connections::close_int_conn | impala connection closed.")
      }
      # Check for sparklyr connection
      else if(any(stringr::str_detect(class(sc), "spark_connection"))){
        sparklyr::spark_disconnect(conn)
        futile.logger::flog.debug("tidy_connections::close_int_conn | spark connection closed.")
      }
    }
  }
  else {
    futile.logger::flog.debug("tidy_connections::close_int_conn | internal connection not detected.")
  }
}

#' whoami - Determine user name/id for working with CDH. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @import futile.logger
#' @import glue
#' @import dplyr
#' @import stringr
whoami <- function(){
  names <-dplyr::tibble(env = c("USER", "HADOOP_USER_NAME"))
  
  names <- 
    names %>% 
    dplyr::mutate(getenv = purrr::map_chr(env, ~Sys.getenv(., unset = NA)))
  
  futile.logger::flog.debug("whoami | Determined System environmental variables to be:")
  futile.logger::flog.debug(glue::glue("whoami | {names}"))
  
  ret <- 
    names %>% 
    dplyr::filter(!stringr::str_detect(getenv, "cdsw") & !is.na(getenv)) %>% 
    dplyr::distinct(getenv) %>% 
    dplyr::pull(getenv) %>% 
    stringr::str_to_lower()
  
  futile.logger::flog.debug(glue::glue("whoami | Determined username to be: {ret}"))
  
  if(!stringr::str_detect(ret, "\\w+")){
    futile.logger::flog.error("whoami | username appears empty. Email kbs14104@gsk.com to report the bug.")
  }
  
  return(ret)
}

# TODO: Add debug comments
#' Drop (tmp) impala_copy_to tables.
#' 
#' \strong{drop_impala_copy() - drop tmp tables}
#' This function will easily drop \code{tbl} referenecs that were created by \code{impala_copy_to}.
#' This function will *not* drop other tables. 
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}
#' @export
#' @param .table Table reference to drop
#' @import futile.logger
#' @import purrr
#' @import glue
#' @import dplyr
#' @import stringr
drop_impala_copy <- function(.table = NULL){
  if(is.null(.table)){
    futile.logger::flog.error("tidy_connections::drop_tmp_impala_tbl | no .table specified.")
    return();
  } 
  
  table_path <- purrr::pluck(.table, "ops") %>% purrr::pluck("x", 1)
  if(!stringr::str_detect(table_path, stringr::regex("\\w+\\.\\w+"))){
    futile.logger::flog.error(glue::glue("tidy_connections::drop_tmp_impala_tbl | {table_path} path doesn't match database.table convention."))
  }
  
  # Remove the table IF the table is an impala table
  if("tbl_impala" %in% class(.table)){ 
    tmp_impala_tbl_value <- attr(.table, "tmp_impala_tbl", exact = TRUE)
    # And the table is marked with attribute tmp_impala_tbl == TRUE
    if(!is.null(tmp_impala_tbl_value)){
      if(tmp_impala_tbl_value == TRUE){
        safely_dbExecute <- purrr::safely(implyr::dbExecute)
        sql_statement <- glue::glue("DROP TABLE IF EXISTS {`table_path`} PURGE")
        
        exec <- safely_dbExecute(impala, sql_statement)
        if (!is.null(exec$error)){
          futile.logger::flog.error(glue::glue("tidy_connections::drop_tmp_impala_tbl | unable to remove {`table_path`} because:\n{exec$error}"))
          stop();
        }
        else {
          futile.logger::flog.debug(glue::glue("tidy_connections::drop_tmp_impala_tbl | Successfully removed {`table_path`}"))
        }
      }
      else{
        futile.logger::flog.warn("tidy_connections::drop_tmp_impala_tbl | This is not a tbl flagged as tmp by impala_copy_to. Skipping.")
      }
    }
    else{
      futile.logger::flog.warn("tidy_connections::drop_tmp_impala_tbl | This is not a tbl flagged as tmp by impala_copy_to. Skipping.")
    }
  }
  else{
    futile.logger::flog.warn("tidy_connections::drop_tmp_impala_tbl | This is not an impala tbl. Skipping.")
  }
}
