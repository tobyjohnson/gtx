\name{multimatch}
\alias{multimatch}
\title{Match with multiple matching possible.}
\description{
  Function to represent results of matching when each query argument may match
  multiple elements in the target.
}
\usage{
multimatch(query, target, values, sep = ",", use.unique = TRUE)
}
\arguments{
  \item{query}{Vector of values to be matched.}
  \item{target}{Vector of values that \code{query} is to be matched against.}
  \item{values}{Vector of values to be returned for matches in \code{target}.}
  \item{sep}{Character to separate values for multiple matches.}
  \item{use.unique}{Logical flag indicating whether to return unique set
    of values only.}
}
\value{
  A character vector of length \code{length(query)}.
}
\examples{
bmidata <- data.frame(subject = c("A001", "A002", "A003", "A003"),
                      year = c(2001, 2001, 2001, 2005),
                      bmi = c(21.3, 29.7, 25.5, 22.3))
multimatch(c("A002", "A003"), bmidata$subject, bmidata$bmi)
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
