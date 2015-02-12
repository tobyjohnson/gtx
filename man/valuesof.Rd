\name{valuesof}
\alias{valuesof}
\title{Convert character vector to scalar}
\description{
  Convenience function to summarize a character vector (or argument that
  can be coerced), by converting to a character scalar containing the
  unique values present.
}
\usage{
valuesof(x, sep = ",", na.convert = NA)
}
\arguments{
  \item{x}{a character vector of allele frequencies}
  \item{sep}{character used to separate values of \code{x} in the return
    value}
  \item{na.convert}{character used to represent NA values in \code{x}.
    After conversion, NA values are removed}
}
\value{
  A character vector of length 1.
}
\examples{
  valuesof(c(1:5, NA))
  valuesof(c(1:5, NA), na.convert = "NA")
  valuesof(factor(c("Drug A", "Drug A", "Drug A")), sep = "/")
  valuesof(factor(c("Drug A", "Drug B", "Drug A", "Placebo")), sep = "/")
  valuesof(matrix(rbinom(100, 2, 0.5), ncol = 10))
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
