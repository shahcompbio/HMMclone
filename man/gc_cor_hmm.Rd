% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gc.R
\name{gc_cor_hmm}
\alias{gc_cor_hmm}
\title{GC bias correction using HMMcopy approach}
\usage{
gc_cor_hmm(counts, gc, span1 = 0.03, span2 = 0.3, valid = NULL, bin_ids = NULL)
}
\arguments{
\item{counts}{Numeric vector of read counts}

\item{gc}{Numeric vector of GC content values (0-1) corresponding to counts}

\item{span1}{Span parameter for first loess regression (default: 0.03)}

\item{span2}{Span parameter for second loess regression (default: 0.3)}

\item{valid}{Logical vector indicating valid bins (currently unused)}

\item{bin_ids}{Character vector of bin identifiers (currently unused)}
}
\value{
Numeric vector of GC-corrected read counts
}
\description{
Corrects GC bias in read count data using a two-stage loess regression
approach similar to the HMMcopy package.
}
\details{
The correction uses a two-stage approach:
\enumerate{
\item Fit initial loess regression of counts vs GC content
\item Fit second loess regression to smooth the correction curve
\item Divide original counts by predicted values to get corrected counts
}

This approach is adapted from the HMMcopy package methodology.
}
\examples{
# Generate example data
gc <- runif(1000, 0.3, 0.7)
counts <- rpois(1000, lambda = 100 * (1 + 2 * (gc - 0.5)^2))  # GC bias

# Apply GC correction
corrected <- gc_cor_hmm(counts, gc)

}
\seealso{
\code{\link{gc_cor_modal}} for modal regression approach
}
