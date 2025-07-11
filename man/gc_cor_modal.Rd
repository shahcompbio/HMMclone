% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gc.R
\name{gc_cor_modal}
\alias{gc_cor_modal}
\title{GC bias correction using modal regression}
\usage{
gc_cor_modal(
  counts,
  gc,
  valid = rep(TRUE, length(counts)),
  bin_ids = names(counts),
  lowess_frac = 0.2,
  n_sample_data = NULL,
  degree = 4,
  knots = c(0.38),
  q = c(0.1, 0.9),
  g = c(0.1, 0.9),
  model = "poly",
  results = c("counts", "default", "full")
)
}
\arguments{
\item{counts}{Numeric vector of read counts}

\item{gc}{Numeric vector of GC content values (0-1) corresponding to counts}

\item{valid}{Logical vector indicating which bins to use (default: all TRUE)}

\item{bin_ids}{Character vector of bin identifiers (default: sequence numbers)}

\item{lowess_frac}{Fraction for lowess smoothing (default: 0.2)}

\item{n_sample_data}{Number of data points to sample for quantile regression (default: NULL, use all)}

\item{degree}{Degree for spline fitting (default: 4, currently unused)}

\item{knots}{Knot positions for spline fitting (default: 0.38, currently unused)}

\item{q}{Quantile range for regression c(lower, upper) (default: c(0.1, 0.9))}

\item{g}{GC content range for integration c(lower, upper) (default: c(0.1, 0.9))}

\item{model}{Model type (default: "poly", currently unused)}

\item{results}{Return format: "counts", "default", or "full" (default: "counts")}
}
\value{
Depends on \code{results} parameter:
\itemize{
\item "counts": Numeric vector of GC-corrected read counts
\item "default": Data frame with correction details (excluding quantile curves)
\item "full": Complete data frame with all quantile regression results
}
}
\description{
Corrects GC bias using modal regression approach with quantile regression
to find the modal relationship between read counts and GC content.
}
\details{
This method uses quantile regression to model the relationship between read counts
and GC content across multiple quantiles. The modal quantile (representing the
most common relationship) is identified and used for correction.

The approach:
\enumerate{
\item Fits polynomial quantile regression across multiple quantiles
\item Integrates under each quantile curve
\item Finds the modal quantile with minimum distance between adjacent curves
\item Uses the modal curve for GC bias correction
}

This implementation is based on methods from Matthew Zatzman (scatools).
}
\examples{
# Generate example data with GC bias
gc <- runif(1000, 0.3, 0.7)
counts <- rpois(1000, lambda = 100 * (1 + 2 * (gc - 0.5)^2))

# Apply modal GC correction
corrected <- gc_cor_modal(counts, gc)

# Get full results
full_results <- gc_cor_modal(counts, gc, results = "full")

}
\seealso{
\code{\link{gc_cor_hmm}} for HMMcopy-style correction
}
