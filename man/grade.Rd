\name{grade}
\alias{grade}
\title{Convert numeric values to ordinal grades}
\description{
  Convenience function to convert numeric values to ordinal grades at
  user defined threshold values.
}
\usage{
grade(x, threshvals, strict = TRUE)
}
\arguments{
  \item{x}{a vector of numeric values}
  \item{threshvals}{a vector of threshold values}
  \item{strict}{Apply strict inequalities}
}
\value{
  The vector of ordinal grades.
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
