\name{plantation}
\alias{plantation}
\title{Plantation (regular forest) plot}
\description{
  Provides a simple forest plot.  Allows an arbitrary functional
  relationship between the plotting scale and the numerical values
  displayed, useful hazard ratios or odds ratios where plotting on a log
  scale is preferred.
}
\usage{
plantation(ntext, beta, ci.lo, ci.hi, se, alpha = 0.05,
                       FUN = I, pvals = TRUE, meta = TRUE,
                       xzero = 0, xlim, xticks, xlab = "Effect",
                       digits = 2,
                       groups = list("Fixed effect meta-analysis" = 1:length(beta)))
}
\arguments{
  \item{ntext}{Names for each item}
  \item{beta}{Effect sizes on the plotting scale}
  \item{ci.lo}{Lower CI limits}
  \item{ci.hi}{Upper CI limits}
  \item{se}{Standard errors}
  \item{alpha}{Alpha level for CIs}
  \item{FUN}{Function to transform from plotting scale to printing scale}
  \item{pvals}{Whether to write P-values}
  \item{meta}{Calculate meta-analysis}
  \item{xzero}{Value of X at which to }
  \item{xlim}{X-axis limits on plotting scale}
  \item{xticks}{X-axis ticks on plotting scale}
  \item{xlab}{X-axis label}
  \item{digits}{digits}
  \item{groups}{List of groups for meta-analyses}
}
\details{
  The user can supply \emph{either} standard errors (SEs) or upper and
  lower limits of confidence intervals (CIs).  Whichever is missing is
  calculated from the other, assuming the CIs are symmetric and at the
  level \code{alpha}.

  A key concept is the separation of the printing scale (which is the
  scale that effect sizes are conventionally reported on), and the
  plotting scale (which should typically be the same as the analysis
  scale; the scale on which normal theory for the estimator best holds).
  The most common examples are the printing scale being hazard ratio
  (HR) or odds ratio (OR), and the plotting scale being (natural) log HR
  or log OR.  A key feature of \code{plantation} is that it
  allows an arbitrary function \code{FUN} to map the plotting scale onto
  the printing scale.
  
  See the examples and vignette.
} 
\examples{
plantation(paste("Study ", 1:4, " (drug ", c("A", "A", "B", "B"), ")", sep = ""),
           c(0.41, 0.06, 0.39, 0.50),
           se = c(0.12, 0.17, 0.18, 0.25),
           FUN = exp, digits = 2,
           xlim = log(c(0.8, 2.5)), xticks = log(c(0.8, 1, 1.5, 2, 2.5)),
           xlab = "HR", 
           groups = list("Fixed effect all drug A" = 1:2,
             "Fixed effect all drug B" = 3:4,
             "Fixed effect all studies" = 1:4))
}
\author{
  Toby Johnson \email{Toby.x.Johnson@gsk.com}
}
