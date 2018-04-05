\name{gtxregion}
\alias{gtxregion}
\alias{gtxwhere}
\title{Interface to define a genomic region}
\description{
  Unified interface to define a genomic region.
}
\usage{
gtxregion(chrom, pos_start, pos_end, 
          hgncid, ensemblid, rs, pos, surround = 500000, 
          dbc = getOption("gtx.dbConnection", NULL))
gtxwhere(chrom, 
         pos, pos_ge, pos_le, 
         pos_end_ge, pos_start_le, 
         pos_start_ge, pos_end_le, 
         rs, hgncid, ensemblid,
         tablename)
}
\arguments{
  \item{chrom}{Character specifying chromosome}
  \item{pos_start}{Integer start position on chromosome}
  \item{pos_end}{Integer end position on chromosome}
  \item{pos}{Integer position on chromosome}
  \item{hgncid}{HGNC gene identifier}
  \item{ensemblid}{ENSEMBL gene identifier}
  \item{rs}{dbSNP rs identifier}
  \item{surround}{Distance around entity to include in region}
  \item{pos_ge}{Position greater-or-equal required}
  \item{pos_le}{Position less-or-equal required}
  \item{pos_end_ge}{End position greater-or-equal required}
  \item{pos_start_le}{Start position less-or-equal required}
  \item{pos_start_ge}{Start position greater-or-equal required}
  \item{pos_end_le}{End position less-or-equal required}
  \item{tablename}{Database table name}
  \item{dbc}{Database connection}
}
\details{

  The \code{gtxregion()} function provides a unified interface for other
  functions to define a genomic region (or potentially for a user to
  invoke directly).  For any valid combination of its optional
  arguments, it returns genomic coordinates (chromosome, start and end
  positions) as described below, using the database connection
  \code{dbc} to resolve any queries (such as the coordinates of a named
  gene).

  When accessing this functionality indirectly via higher level
  functions (such as \code{\link{regionplot}()} and
  \code{\link{coloc}()}), the functionality should be almost completely
  intuitive for most users, and if necessary can be learned by example
  from the manual pages and vignettes for those higher level functions.
  It suffices to add that the optional arguments are used according to a
  priority order, which is exactly the order of arguments in the
  function definition.  For example if \code{chrom}, \code{pos_start},
  \code{pos_end} and \code{hgnc} are all provided, \code{hgnc} has lower
  priority and is ignored.  Similarly if \code{hgnc} and \code{pos} are
  provided, \code{pos} has lower priority and is ignored.

  It is an intended design feature that \code{pos} and \code{rs} are
  lowest in the priority order.  When used in conjunction with higher
  priority arguments such as \code{hgnc}, a \code{pos} or \code{rs}
  argument can be used \emph{without} affecting the genomic region
  specified, which then allows a function that wraps \code{gtxregion()}
  to use \code{pos} or \code{rs} for secondary purposes, such as to
  highlight a specific position or variant in a visual display.  Thus,
  \code{regionplot(..., pos = 1234567, surround = 500000)} selects a
  500kb region around position 1234567 and visually highlights any
  variant present at position 1234567, and \code{regionplot(..., hgnc =
  'ABC123', surround = 10000, pos = 1234567)} selects a 10kb region
  around the ABC123 gene and visually highlights any variant present at
  position 1234567.
  
  The remainder of this manual page is more technical documentation,
  intended for programmers writing new high level functions that will
  work alongside \code{\link{regionplot}()} and \code{\link{coloc}()},
  and should be read in combination with the source code.
  
  The \code{gtxregion()} function resolves its arguments to genomic
  coordinates as follows:
  
  If the arguments \code{chrom}, \code{pos_start} and \code{pos_end} are
  all provided, these are checked for validity and used to directly
  specify the return value.

  Otherwise, if the argument \code{hgnc} is provided, \code{TABLE genes}
  is queried (using \code{dbc} and \code{gtxwhere}) and a region
  spanning the gene(s) plus \code{surround}ing distance is returned.
  Otherwise, if the argument \code{ensg} (integer) is provided,
  \code{TABLE genes} is similarly queried.

  Otherwise, if the arguments \code{chrom} and \code{pos} are both
  provided, these are checked for validity and used plus
  \code{surround}ing distance to directly specify the return value.

  Otherwise, if the argument \code{rs} is provided, \code{TABLE sites}
  (\code{sites_by_rs}) is queried (using \code{gtxwhere}) and a region
  plus \code{surround}ing distance is returned.

  The methods just described are implemented using \code{if ... else if
  ... else if ...} logic, so for example if a \code{hgnc} argument is
  provided then any \code{ensg} argument is ignored, etc.

  The \code{gtxwhere} function provides a standardized and sanitized way
  to dynamically construct part of a SQL WHERE statement.  This is best
  illustrated by the examples below.  When more than one argument value
  is given, either as multiple values for a single argument, or for more
  that one argument, the following logic seems most useful: Multiple
  values for a single argument are combined or-wise, and multiple
  arguments are combined and-wise.

  To use \code{gtxwhere} to select chromosome segments (such as genes or
  other entities, recombination rate segments, etc) that wholly \emph{or
  partially} overlap a query region, use \code{pos_end_ge=query_start}
  and \code{pos_start_le=query_end}.  To select \emph{only} chromosome
  segments that wholly overlap, instead use
  \code{pos_start_ge=query_start} and \code{pos_end_le=query_end}.

  For identifiers that are represented (for efficiency) as integers in
  database tables but as strings in \dQuote{user space}, \code{gtxwhere}
  is the layer at which string-to-integer checking and conversion should
  occur.
}
\value{
  \code{gtxregion} returns a named list with elements \sQuote{chrom}
  (character), \sQuote{pos_start} (integer) and \sQuote{pos_end}
  (integer).

  \code{gtxwhere} returns a character string suitable for inclusion
  after the WHERE clause in a SQL statement.
}
\examples{
\dontrun{
  gtxregion(chrom = 1, pos_start = 109616403, pos_end = 109623689)
  # dies without an open ODBC connection
}
gtxwhere(rs = 'rs599839')
gtxwhere(chrom = 1, pos = c(109616403, 109623689))
gtxwhere(chrom = 1, pos_end_ge = 109616403, pos_start_le = 109623689)
\dontrun{
gtxwhere()
}
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
