\name{gtx.version}
\alias{gtx.version}
\docType{data}
\title{Detailed package version information}
\description{
Detailed package information from git repository log immediately prior
to package building.
}
\usage{data(gtx.version)}
\format{
  A character vector with three rows, which are the SHA1 sum, the author
  name and email, and the date, of the last commit in the git repository
  from which this package was built.

  A custom script is used to update the value of
  \code{gtx.version} each time the package is built.

  The SHA1 sum can be used to obtain the package source code from the
  git repository at \url{https://github.com/tobyjohnson/gtx}.  Therefore,
  printing the contents of \code{gtx.version} into analysis log files is
  an important part of documenting and ensuring reproducibility of the
  analysis.
}
\examples{
data(gtx.version)
print(gtx.version)
}
\keyword{datasets}
