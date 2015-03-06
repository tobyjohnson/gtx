\name{pipeline}
\alias{pipeline}
\title{(Pharmaco)genetic analysis and reporting pipeline}
\description{
  Convenience function to convert infinite values to missing.
}
\usage{
pipeline(config = "config.R", slave)
}
\arguments{
  \item{config}{name of a configuration file}
  \item{slave}{used in recursive calling}
}
\description{
  This is a work in progress.
  
  Runs a (to be) complete (pharmaco)genetic analysis and reporting
  pipeline.  The premise is that clinical trial data are available as
  SAS datasets exported in plain text format using SAS proc export, and
  genotype data are available (currently assumed to be MACH/minimac
  dosage files).  Given option settings that
  describe in a particular way the (pharmaco)genetic analyses to be
  conducted, a single \code{pipeline()} function call will conduct the
  analyses and generate a report and source displays (plots and tables).

  The function is controlled by options and variables set when the
  configuation file is \code{source}d.

  Essential variables to set include pgx.models, pgx.groups,
  pgx.transformations, pgx.eigenvec, clinical.derivations

  Derivations are (here) by definition rules for converting clinical
  data (which may not be one-row-per-subject) to one-row-per-subject
  variables, and which can be applied independently over subjects.
  I.e. the derived variable for the i-th subject is independent of data
  for other subjects.

  Transformations are (here) by definition rules for converting one
  one-row-per-subject variable to another, which may depend on the data
  over all subjects in a given analysis dataset.  E.g. rank or quantile
  transformations, converting text to ordered factors.

  
  
  
  Optional variables to set include option(clinical.usubjid),
  genotype.path, clinical.path,

  TO do, allow control of make command.
  Candidate gene lists by analysis.  Better interface for specifying the
  options.   Currently deps have to be specified manually.  The deps
  have two purposes, firstly to make data loading and derivations more
  efficient by only loading/deriving the data needed, and secondly to
  compute analysis datasets and subsets of subjects with nonmissing data
  for each analysis.
  
  WARNING: The master/slave arangement means BAD THINGS MAY HAPPEN if
  you upgrade the gtx package while \code{pipeline()} is running.

  The reasons that options are specified in a configuration file instead of as
  function arguments

  Analysis datasets are stored as csv files with two special comment
  rows.  Hence R classes such as Surv and ordered factors are not
  preserved.  (This is a problem?)  Can write models like
  coxph(Surv(SRVMO,SRVCFLCD)~...) and
  clm(factor(BR,c("CR","PR","SD","PD"))~...)
  but this is tedious and inefficient.
  Consider applying transformations within slave process?
}
\value{
  Returns an invisible NULL.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
