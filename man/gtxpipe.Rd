\name{gtxpipe}
\alias{gtxpipe}
\title{Pipeline for routine genetic association analysis and reporting}
\description{
  An implementation of a pipeline that simplifies and standardizes the
  analysis and report generation for routine genetic association projects.
}
\usage{
gtxpipe(gtxpipe.models = getOption("gtxpipe.models"),
        gtxpipe.groups = getOption("gtxpipe.groups", data.frame(group = 'ITT', deps = 'pop.PNITT', fun = 'pop.PNITT', stringsAsFactors = FALSE)),
        gtxpipe.derivations = getOption("gtxpipe.derivations", {data(derivations.standard.IDSL); derivations.standard.IDSL}),
        gtxpipe.transformations = getOption("gtxpipe.transformations", data.frame(NULL)),
        gtxpipe.eigenvec,
        stop.before.make = FALSE)
}
\arguments{
  \item{gtxpipe.models}{A data frame defining association models to be
    fitted}
  \item{gtxpipe.groups}{A data frame defining (sub)groups of
    individuals in which to fit the models}
  \item{gtxpipe.derivations}{A data frame defining methods to derive
    analysis variables from underlying clinical data}
  \item{gtxpipe.transformations}{A data frame defining transformations
    required for the analysis variables}
  \item{gtxpipe.eigenvec}{A filename containing eigenvectors or other
    covariates to adjust for in all association models}
  \item{stop.before.make}{Logical whether to stop before actually
    fitting association models}
}
\details{
  The pipeline implemented by \code{gtxpipe} takes as input some
  clinical data and some genotype data, performs association analyses,
  and outputs tables, figures, and documents summarising the results.

  The association analyses that are conducted are controlled by function
  arguments passed to \code{gtxpipe}, and options with names that begin
  \sQuote{\code{gtx.}} or \sQuote{\code{gtxpipe.}}.  The intention is
  that settings that are likely to vary on a project by project basis
  are controlled by function arguments, and settings that are likely to
  be constant over multiple projects are controlled by options.

  The input clinical data must be a set of plain text files inside a single
  directory, which is specifed by \code{options(gtxpipe.clinical)} (default a
  subdirectory \sQuote{clinical} of the current working directory).  [In future
  the intent is to support multiple directories for multiple clinical
  studies to be analysed together.]  These are expected to have
  \code{.txt} extenstions and to be 
  SAS datasets exported in plain text format using SAS \code{proc
  export}, but other plain text files in the same format may work.  The
  files are read using \code{gtx::\link{clinical.import}}.

  The input genotype data must be a set of minimac dosage and
  information files, inside a single directory specifed by
  \code{options(gtxpipe.genotypes)} (default a subdirectory
  \sQuote{genotypes} of the current working directory).  [In future the
  intent is to support other genotype data formats].

  Analyses are conducted inside a directory specifed by
  \code{options(gtxpipe.analyses)} (default a subdirectory
  \sQuote{analyses} of the current working directory), and the pipeline
  outputs are written inside a directory specifed by
  \code{options(gtxpipe.outputs)} (default a subdirectory
  \sQuote{outputs} of the current working directory).  Both directories
  are created if they do not already exist.

  \bold{Models} are defined (for the purposes of \code{gtxpipe}) as
  regression models on derived and transformed clinical data
  (derivations and transformations are defined below).  Genetic
  variables are not explicitly included in the model definition; these
  are automatically added by \code{gtxpipe}.  The \code{gtxpipe.models}
  argument must be a dataframe with the following columns (with the
  following information in each row): \dQuote{model} (a name used for
  organising and reporting results); \dQuote{deps} (a string with space
  delimited names of variables the model depends on); \dQuote{fun} (a
  string with an R language statement of the model); \dQuote{groups} (a
  string with space delimited subject groups [as defined in Groups
  below] in which to evaluate the model); \dQuote{contrasts} (a string
  with space delimited group contrasts of interest); \dQuote{cvlist} (a
  string with space delimited identifiers for candidate genetic variants
  of interest).  [In future the intent is to provide an interface so
  that the information in \dQuote{deps} is determined automatically from
  a specification of \dQuote{fun}.]
    
  \bold{Derivations} are defined (for the purposes of \code{gtxpipe}) as
  methods that: (i) convert clinical data (which may not be
  one-row-per-subject) to analysis variables (which must be
  one-row-per-subject); and (ii) can be applied independently over
  subjects.  That is, a derived analysis variable for the i-th subject
  must be a scalara that depends \emph{only} on the clinical data for
  the i-th subject.  An example of a derivation is to compute the
  highest grade of adverse event experienced by each subject, using the
  clinical data for all adverse events of a given type.  The
  \code{gtxpipe.derivations} argument must be a dataframe suitable for
  passing as the \code{derivations} argument to the function
  \code{\link{clinical.derive}}.  Note that all analysis variables must
  be derived (even if the derivation is a simple extraction from a
  clinical dataset).  See the examples for how to write simple
  derivations.  For efficiency, it is possible to specify the derivation
  of multiple variables in a single row of \code{gtxpipe.derivations},
  as long as the \dQuote{deps} and \dQuote{data} part of the derivation
  are the same for all the variables.

  \bold{Transformations} are defined (for the purposes of \code{gtxpipe}) as
  methods that convert one or more analysis variables (which must be
  one-row-per-subject) to another analysis variable (also
  one-row-per-subject).  Transformations may act independently over
  subjects (e.g. a log transform), but also may act such that the
  transformed value for the i-th subject depends on the values for other
  subjects (e.g. a rank or quantile transform).  Thus in general the
  action of the transformation may depend on which set or subset of
  subjects it is applied to.  Transformations may also be used to
  convert data types e.g. to create a variable of the \code{Surv} class
  from a time variable and a censoring indicator variable.  The
  \code{gtxpipe.transformations} argument must be a dataframe with the
  following columns (with the following information in each row):
  \dQuote{targets} (a name [or space delimited names] for the derived
  variable[s]); \dQuote{deps} (a string with space delimited
  names of variables the derivation depends on); \dQuote{fun}
  (a string with an R language statement of the transformation).  

  \bold{Groups} are defined using R language statements that depend on
  derived variables. FIXME: Can they depend on transformations?  The
  \code{gtxpipe.groups} argument must be a dataframe with the following
  columns (with the following information in each row): \dQuote{group}
  (a name for the group, as referred to in Model groups and contrasts);
  \dQuote{deps} (a string with space delimited names of variables the
  derivation depends on); \dQuote{fun} (a string with an R language
  statement that evaluates logical for group membership).  \bold{Group
  contrasts} are specified in Models using the names of two groups,
  separated by a \sQuote{/}.  For this reason, and because group names
  are used as directory names, all group names must be simple alphanumeric.

  \bold{Descriptors} is a name for the concept of storing a table of long
  descriptive names for individual variables and substituting those
  names in output displays.  FIXME Document more.

  make command needs to be set.  FIXME should be a sensible default.
  
  FIXME list:  Candidate variants should be analysed even if failing MAF
  or Rsq filters.  Better interface for specifying the arguments (helper
  functions to build models, derivations etc in the above format.
  Currently deps have to be specified manually but this should be
  automatic for most cases).

  FIXME:  Analysis datasets are stored as csv files with two special comment
  rows.  Hence R classes such as Surv and ordered factors are not
  preserved.  (This is a problem.)  Can write models like
  coxph(Surv(SRVMO,SRVCFLCD)~...) and
  clm(factor(BR,c("CR","PR","SD","PD"))~...)
  but this is tedious and inefficient.
  Consider applying transformations within slave process?
  
  Note regarding deps: The deps have two purposes, firstly to make data
  loading and derivations more efficient by only loading/deriving the
  data needed, and secondly to compute analysis datasets and subsets of
  subjects with nonmissing data for each analysis.

  WARNING: The master/slave arrangement means BAD THINGS MAY HAPPEN if
  you upgrade the gtx package while \code{pipeline()} is running.

} 
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
