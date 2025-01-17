% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/padog.R
\name{fit_ebayes}
\alias{fit_ebayes}
\title{Perform eBayes analysis from limma.}
\usage{
fit_ebayes(esetm, contrasts, mod, pdata, rna_seq = FALSE)
}
\arguments{
\item{contrasts}{comparison to make. For \code{\link{padog}} is \code{'d-c'}.}

\item{mod}{design matrix}

\item{pdata}{data.frame with columns \code{lib.size} and \code{norm.factors} needed if \code{rna_seq} is \code{TRUE}.}

\item{rna_seq}{is the analysis on RNA seq data? Default is \code{FALSE}. If \code{TRUE} must supply \code{pdata}.}
}
\value{
result from call to limma \code{eBayes}.
}
\description{
Generates contrast matrix then runs eBayes analysis from limma.
}
\keyword{internal}
