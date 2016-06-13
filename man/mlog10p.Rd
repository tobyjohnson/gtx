\name{mlog10p}
\alias{mlog10p}
\title{Negative log10 P values}
\description{
  Compute minus log10 P values from effect sizes and standard errors.
}
\usage{
mlog10p(beta, se)
}
\arguments{
  \item{beta}{a vector of effect sizes}
  \item{se}{a vector of standard errors}
}
\value{
  The vector of minus log10 P values.
}
\details{
  The results are the same as \code{-log10(pnorm(-abs(beta)/se)*2)}, but
  are calculated to machine precision even when the P values are very
  small (such that the intermediate step of calculating the P values
  would give zero to machine precision).
}
\examples{
beta <- c(1:4, (1:10)*5) * rep(c(-1, 1), 7)
se <- rep(1, length(beta))
pval1 <- pnorm(-abs(beta)/se)*2
pval2 <- pchisq(beta^2/se^2, df = 1, lower.tail = FALSE)
cbind(z = beta/se, pval = pval1, -log10(pval1), -log10(pval2), mlog10p(beta, se))
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
