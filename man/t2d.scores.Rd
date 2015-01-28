\name{t2d.scores}
\alias{t2d.scores}
\docType{data}
\title{Genetic risk scores for type 2 diabetes risk.}
\description{
  Risk score parameterised using GWAS meta-analysis results published by
  the DIAGRAM+ consortium (Voight et al. Nature Genetics 2011).  The
  score is for risk of prevalent type 2 diabetes (T2D).
}
\usage{
data(t2d.scores)
}
\format{
  A data frame suitable for use with other functions in this package,
  see \link{gtx.params}.

  The score is called T2D2010, and has coefficients in ln(odds) for T2D per coded
  allele.
  
  In future releases, additional risk scores used in earlier papers may be added for
  historical interest (WHII and Lin papers).
}
\source{
  For the publication by Voight et al. from which these data were extracted see \url{http://dx.doi.org/10.1038/ng.609}.
}
\keyword{datasets}