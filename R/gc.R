# Different way similar to HMMcopy
# Code adapated from `HMMcopy` package
#' GC bias correction using HMMcopy approach
#'
#' Corrects GC bias in read count data using a two-stage loess regression
#' approach similar to the HMMcopy package.
#'
#' @param counts Numeric vector of read counts
#' @param gc Numeric vector of GC content values (0-1) corresponding to counts
#' @param span1 Span parameter for first loess regression (default: 0.03)
#' @param span2 Span parameter for second loess regression (default: 0.3)
#' @param valid Logical vector indicating valid bins (currently unused)
#' @param bin_ids Character vector of bin identifiers (currently unused)
#'
#' @return Numeric vector of GC-corrected read counts
#'
#' @details
#' The correction uses a two-stage approach:
#' 1. Fit initial loess regression of counts vs GC content
#' 2. Fit second loess regression to smooth the correction curve
#' 3. Divide original counts by predicted values to get corrected counts
#'
#' This approach is adapted from the HMMcopy package methodology.
#'
#' @examples
#' # Generate example data
#' gc <- runif(1000, 0.3, 0.7)
#' counts <- rpois(1000, lambda = 100 * (1 + 2 * (gc - 0.5)^2))  # GC bias
#'
#' # Apply GC correction
#' corrected <- gc_cor_hmm(counts, gc)
#'
#' @seealso \code{\link{gc_cor_modal}} for modal regression approach
#'
#' @export
gc_cor_hmm <- function(counts, gc, span1 = 0.03, span2 = 0.3, valid = NULL, bin_ids = NULL) {
  counts <- as.vector(counts)
  rough <- loess(counts ~ gc, span = span1)
  idx <- seq(0, 1, by = 0.001)
  final <- loess(predict(rough, idx) ~ idx, span = span2)
  counts_cor <- counts / predict(final, gc)
  return(counts_cor)
}

#R implementation of modal regression GC correction from Matthew Zatzman.
# See scatools
#' GC bias correction using modal regression
#'
#' Corrects GC bias using modal regression approach with quantile regression
#' to find the modal relationship between read counts and GC content.
#'
#' @param counts Numeric vector of read counts
#' @param gc Numeric vector of GC content values (0-1) corresponding to counts
#' @param valid Logical vector indicating which bins to use (default: all TRUE)
#' @param bin_ids Character vector of bin identifiers (default: sequence numbers)
#' @param lowess_frac Fraction for lowess smoothing (default: 0.2)
#' @param n_sample_data Number of data points to sample for quantile regression (default: NULL, use all)
#' @param degree Degree for spline fitting (default: 4, currently unused)
#' @param knots Knot positions for spline fitting (default: 0.38, currently unused)
#' @param q Quantile range for regression c(lower, upper) (default: c(0.1, 0.9))
#' @param g GC content range for integration c(lower, upper) (default: c(0.1, 0.9))
#' @param model Model type (default: "poly", currently unused)
#' @param results Return format: "counts", "default", or "full" (default: "counts")
#'
#' @return Depends on \code{results} parameter:
#' \itemize{
#'   \item "counts": Numeric vector of GC-corrected read counts
#'   \item "default": Data frame with correction details (excluding quantile curves)
#'   \item "full": Complete data frame with all quantile regression results
#' }
#'
#' @details
#' This method uses quantile regression to model the relationship between read counts
#' and GC content across multiple quantiles. The modal quantile (representing the
#' most common relationship) is identified and used for correction.
#'
#' The approach:
#' 1. Fits polynomial quantile regression across multiple quantiles
#' 2. Integrates under each quantile curve
#' 3. Finds the modal quantile with minimum distance between adjacent curves
#' 4. Uses the modal curve for GC bias correction
#'
#' This implementation is based on methods from Matthew Zatzman (scatools).
#'
#' @examples
#' # Generate example data with GC bias
#' gc <- runif(1000, 0.3, 0.7)
#' counts <- rpois(1000, lambda = 100 * (1 + 2 * (gc - 0.5)^2))
#'
#' # Apply modal GC correction
#' corrected <- gc_cor_modal(counts, gc)
#'
#' # Get full results
#' full_results <- gc_cor_modal(counts, gc, results = "full")
#'
#' @seealso \code{\link{gc_cor_hmm}} for HMMcopy-style correction
#'
#' @export
gc_cor_modal <- function(counts,
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
                         results = c("counts", "default", "full")) {
  if (length(counts) != length(gc)) {
    stop("Length of counts and gc vectors are not identical")
  }

  if (length(q) != 2 | (q[1] > q[2])) {
    stop("Must provide lower and upper quantile")
  }

  if (is.null(bin_ids)) {
    bin_ids <- seq_len(length(counts))
  }

  results <- match.arg(results)

  # Into a dataframe pre-filtering
  df_regression_all <- data.frame(reads = counts, gc = gc, bin_ids = bin_ids, valid = valid)

  # Filter the data for non-zero gc and read counts
  df_regression <- df_regression_all[df_regression_all$valid == TRUE, ]

  # 2nd order polynomial quantile regression
  quantiles <- seq(q[1], q[2], by = 0.01)
  quantile_names <- paste0("q", quantiles * 100)

  if (nrow(df_regression) < 10) {
    df_regression <- .gc_warn_error(df_regression, quantile_names)
  } else {
    # Fit second order polynormial quantile regression
    if (is.null(n_sample_data)){
      poly2_quantile_model <- quantreg::rq(reads ~ gc + I(gc^2),
                                           tau = quantiles, data = df_regression, method = "fn")
    } else{
      poly2_quantile_model <- quantreg::rq(reads ~ gc + I(gc^2),
                                           tau = quantiles, data = dplyr::sample_n(df_regression, n_sample_data), method = "fn")
    }
    #poly2_quantile_model <- quantreg::rq(reads ~ splines2::bSpline(gc, df = degree, intercept = T), tau = quantiles, data = df_regression, method = "fn")


    # Fit to our data
    poly_quantile_predict <- predict(object = poly2_quantile_model, newdata = df_regression)

    colnames(poly_quantile_predict) <- quantile_names # rename columns

    # Bind to our data
    df_regression <- dplyr::bind_cols(df_regression, poly_quantile_predict)

    poly2_quantile_params <- as.data.frame(poly2_quantile_model$coefficients)
    colnames(poly2_quantile_params) <- quantile_names # rename columns

    # integration and mode selection
    gc_min <- quantile(df_regression$gc, g[1])
    gc_max <- quantile(df_regression$gc, g[2])

    true_min <- min(df_regression$gc)
    true_max <- max(df_regression$gc)

    poly2_quantile_integration <- c(0, apply(X = poly2_quantile_params, MARGIN = 2, FUN = function(params) {
      poly2 <- polynom::polynomial(params)
      integ <- polynom::integral(poly2)
      integrand <- predict(integ, gc_max) - predict(integ, gc_min)
    }))

    # find the modal quantile
    distances <- diff(poly2_quantile_integration)

    df_dist <- data.frame(quantiles = quantiles, quantile_names = quantile_names, distances = distances)
    dist_max <- quantile(df_dist$distances, 0.95)
    df_dist_filter <- df_dist[df_dist$distances < dist_max, ]

    # Error catch when df_dist_filter is empty
    if (nrow(df_dist_filter) == 0) {
      stop("Not enough data points...")
    } else {
      df_dist_filter$lowess <- lowess(y = df_dist_filter$distances, x = df_dist_filter$quantiles, f = lowess_frac, delta = 0)$y

      modal_quantile <- df_dist_filter[which.min(df_dist_filter$lowess), 2]

      # add values to table
      df_regression["modal_quantile"] <- modal_quantile
      df_regression["modal_curve"] <- df_regression[modal_quantile]
      df_regression["modal_corrected"] <- df_regression["reads"] / df_regression[modal_quantile]
    }
  }



  # Below zeroes are NA? Not sure how to handle exactly
  df_regression[(df_regression$modal_corrected < 0 | is.na(df_regression$modal_corrected)), "modal_corrected"] <- NA

  # Merge back the missing bins to ensure dimensions stay consistent
  df_regression <- dplyr::left_join(df_regression_all, df_regression, by = c("reads", "gc", "bin_ids", "valid"))

  # Do we want to reassign NAs as zeros?
  # df_regression[is.na(df_regression$modal_corrected), "modal_corrected"] <- 0

  if (results == "full") {
    return(df_regression)
  } else if (results == "default") {
    # Only return selected quantile curve
    return(df_regression[, which(!colnames(df_regression) %in% quantile_names)])
  } else if (results == "counts") {
    # only return the corrected counts
    return(as.vector(df_regression$modal_corrected))
  }
}

# Helper function for gc_cor_modal when insufficient data
.gc_warn_error <- function(df_regression, quantile_names) {
  warning("Insufficient data for GC correction. Returning original data with NA for corrected values.")
  df_regression$modal_quantile <- NA
  df_regression$modal_curve <- NA
  df_regression$modal_corrected <- NA
  return(df_regression)
}
