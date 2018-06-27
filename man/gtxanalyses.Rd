\name{gtxanalyses}
\alias{gtxanalyses}
\alias{gtxwhat}
\title{Find GWAS analyses}
\description{
  Find GWAS analyses matching specified criteria.
}
\usage{
gtxanalyses(analysis,
            analysis_not, 
            phenotype_contains,
            description_contains,
            ncase_ge,
            ncohort_ge,
            has_access_only = FALSE, 
            dbc = getOption("gtx.dbConnection", NULL))
gtxwhat(analysis1,
        analysis,
        analysis_not, 
        description_contains,
        phenotype_contains,
        ncase_ge,
        ncohort_ge,
        tablename)
}
\arguments{
  \item{analysis}{Key value(s) for GWAS analysis/es required}
  \item{analysis_not}{Key value(s) for GWAS analysis/es not required}
  \item{phenotype_contains}{Word(s) to search for in phenotype metadata}
  \item{description_contains}{Word(s) to search for in description metadata}
  \item{ncase_ge}{Numeric value that number of cases must be
    greater-or-equal than}
  \item{ncohort_ge}{Numeric value that number in cohort must be
    greater-or-equal than}
  \item{has_access_only}{Logical whether to only return results in
    databases user has access to}
  \item{dbc}{Database connection}
  \item{analysis1}{Key value for exactly one GWAS analysis required}
  \item{tablename}{Table name for SQL query}
}
\details{
  Note, in the near future, a better high level interface will be
  provided to find GWAS analyses of interest.  The \code{gtxanalyses()}
  and \code{gtxwhat()} functions were originally intended as internal
  functions, to abstract part of the functionality required within the
  \code{phewas()} family of functions.  However, in the absence of a
  better interface, these functions are currently exported for users to
  use.
}
\value{
  \code{gtxanalyses} returns a data frame containing selected columns
  from \code{TABLE analyses}.  \code{gtxwhat} returns a string that can
  be used as (part of) the \code{WHERE} clause in an SQL statement.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
		  