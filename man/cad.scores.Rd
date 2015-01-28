\name{cad.scores}
\alias{cad.scores}
\docType{data}
\title{Genetic risk scores for coronary artery disease risk}
\description{
  Risk score parameterised using GWAS meta-analysis results published by
  the C4D and CARDIoGRAM consortia (C4D Consortium Nature Genetics 2011,
  Schunkert et al. Nature Genetics 2011).  The score is for risk of
  prevalent coronary artery disease (CAD).
}
\usage{data(cad.scores)}
\format{
  A data frame suitable for use with other functions in this package,
  see \link{gtx.params}.

  The score is called CAD2011, and has coefficients in ln(odds) for CAD
  per coded allele.  For details of how results from the two papers were
  combined, see manuscript by T. Johnson.
}
\examples{
data(cad.scores)
cad.scores$MAF <- pmin(cad.scores$coded.freq, 1-cad.scores$coded.freq)
cad.scores$OR <- exp(abs(cad.scores$coef))
plot(cad.scores$MAF, cad.scores$OR,
     xlim = c(0, 0.5), ylim = c(1, max(cad.scores$OR)),
     xlab = "Minor allele frequency",
     ylab = "Odds ratio effect", las = 1)
text(cad.scores$MAF, cad.scores$OR, cad.scores$name, pos = 1, cex = 0.5)
}
\source{
  For the publications from which these data were
  extracted see \url{http://dx.doi.org/10.1038/ng.782} and
  \url{http://dx.doi.org/10.1038/ng.784}.
}
\keyword{datasets}