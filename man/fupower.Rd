\name{fupower}
\alias{fupower}
\title{Power for followup studies}
\description{
  Provides a simple calculation of power for followup studies using only
  the discovery P-value and sample size, and followup sample size(s).
}
\usage{
fupower(p0, n0, n1, beta0, invlog = FALSE, alpha = 0.05,
        plty = c("solid", "dashed", "dotdash"), elty = "dotted",
        pcol = "red", ecol = "blue")
}
\arguments{
  \item{p0}{Discovery P-value}
  \item{n0}{Discovery sample size}
  \item{n1}{Vector of followup sample size(s)}
  \item{beta0}{Discovery effect size (optional)}
  \item{invlog}{Whether to label x-axis using inverse log transformed
    values}
  \item{alpha}{Alpha level for one-sided followup test}
  \item{plty}{Power line type, recycled to length of \code{n1}}
  \item{elty}{Effect line type}
  \item{pcol}{Power line colour, recycled to length of \code{n1}}
  \item{ecol}{Effect line colour}
}
\details{
  See the examples.
} 
\examples{
par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(2, 1, 0), mfrow = c(3, 2))
fupower(.05, 1000, c(1000, 300), beta0 = log(1.5), invlog = TRUE)
fupower(5e-4, 1000, c(1000, 300), beta0 = log(1.5), invlog = TRUE)
fupower(5e-8, 1000, c(1000, 300), beta0 = log(1.5), invlog = TRUE)
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
