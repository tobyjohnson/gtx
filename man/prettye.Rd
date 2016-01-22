\name{prettye}
\alias{prettye}
\title{Convert scientific notation to pretty expressions}
\description{
  Converts scientific notation numbers to expressions suitable for
  pretty printing, by converting an e character to an expression for
  times-ten-to-the-power-of.
}
\usage{
prettye(x)
}
\arguments{
  \item{x}{a character vector}
}
\details{
  The argument should be a \emph{character} vector (and is immediately
  coerced using \code{as.character}).  The return value is
  constructed using characters before the \sQuote{e} as the mantissa and
  characters after the \sQuote{e} as the exponent.  Therefore, any
  control of the number of digits to be displayed should be performed
  using e.g. \code{\link{format}} before calling \code{prettye}.
}
\value{
  An expression vector.  (Generated using sapply with individual
  expression return values seems to be converted into an expression with
  multiple elements.)
}
\examples{
p <- 1/(1:10)*10^-(1:10)
plot(1:length(p), -log10(p),
     type = "p", pch = 19, 
     ylim = c(0, max(-log10(p))), las = 1, 
     xlab = "index", ylab = expression(-log[10](P)))
text(1:length(p), -log10(p), prettye(format(p, digits = 3)), pos = 1)
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
