\name{zise}
\alias{zise}
\title{Normal quantile transform data}
\description{
  Transforms data so that it is normally distributed, using only the
  ranks of the data.
}
\usage{
zise(x, only, by)
}
\arguments{
  \item{x}{a vector of values to be transformed}
  \item{only}{optional, logical vector of same length as x}
  \item{by}{optional, factor of same length as x}
}
\details{
  The transformation is sometimes known as a normal quantile
  transform, an inverse normal transformation, or as a quantile-quantile
  normalisation of \code{x}.

  If the argument \code{only} is provided, only the corresponding values
  of \code{x} are transformed, and the remainder are returned as NA.

  If the argument \code{by} is provided, it is coerced to a factor if
  necessary, and the transformation is applied separately 
  by subsets of \code{x} corresponding to each level of \code{by}.
}
\value{
  The vector of transformed values.
}
\examples{
  x <- rbeta(500, shape1 = 0.6, shape2 = 1)
  oldpar <- par(mfrow = c(1, 2))
  hist(x, breaks = 20, col = "red"); box()
  hist(zise(x), breaks = 20, col = "yellow"); box()

  par(mfrow = c(1, 3))
  x <- data.frame(gender = c(rep("Male", 100), rep("Female", 100)),
                  BMI = exp(c(rnorm(100, 3.3, 0.1), rnorm(100, 3.2, 0.2))))
  boxplot(BMI ~ gender, data = x) # not normal
  boxplot(zise(BMI) ~ gender, data = x) # normal overall but not within genders
  boxplot(zise(BMI, by = gender) ~ gender, data = x) # normal overall and within genders

  par(oldpar)
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
