#' Fit Hidden Markov Model for copy number calling
#'
#' Fits a Hidden Markov Model to copy number data using the Viterbi algorithm
#' to find the most likely sequence of copy number states.
#'
#' @param copy Numeric vector of copy number values (e.g., GC-corrected read counts)
#' @param selftransitionprob Probability of staying in the same copy number state (default: 0.99)
#' @param maxCN Maximum copy number state to consider (default: 15)
#' @param sd_value Standard deviation for the normal emission model (default: 0.2)
#' @param usecpp Logical, whether to use C++ implementation of Viterbi algorithm (default: TRUE)
#'
#' @return Integer vector of inferred copy number states (0-based indexing)
#'
#' @details
#' The function models copy number states as hidden states in an HMM where:
#' - States represent integer copy numbers from 0 to maxCN
#' - Emission probabilities follow a normal distribution centered on the true copy number
#' - Transition probabilities favor staying in the same state (selftransitionprob)
#' - The Viterbi algorithm finds the most likely state sequence
#'
#' @examples
#' # Generate some example copy number data
#' copy_data <- c(rep(2, 50), rep(4, 30), rep(2, 20))
#' copy_data <- copy_data + rnorm(length(copy_data), 0, 0.2)
#'
#' # Fit HMM
#' states <- fitHMM(copy_data, selftransitionprob = 0.99, maxCN = 6)
#'
#' @export
fitHMM <- function(copy,
                   selftransitionprob = 0.99,
                   maxCN = 15,
                   sd_value = 0.2,
                   usecpp = TRUE){

  states <- seq(0, maxCN, 1)

  l <- outer(copy, states, FUN = function(x, y) dnorm(x, mean = y, sd = sd_value, log = TRUE))

  transition_prob <- matrix(
    (1 - selftransitionprob) / (length(states) - 1),
    length(states), length(states)
  )
  diag(transition_prob) <- selftransitionprob

  if (usecpp == TRUE){
    res <- viterbi(l,
                   log(transition_prob),
                   observations = seq_len(length(copy)))
  } else{
    res <- viterbiR(l,
                    log(transition_prob),
                    observations = seq_len(length(copy)))
  }

  return(res)
}

HMMcloneperchr <- function(dat,
                     selftransitionprob = 0.99,
                     maxCN = 15,
                     sd_value = 0.2,
                     usecpp = TRUE,
                     pb = NULL){

  if (!is.null(pb)) {
    pb$tick()$print()
  }

  states <- c()
  for (mychr in unique(dat$chr)) {
    hmmresults <- fitHMM(
      copy = dplyr::filter(dat, chr == mychr)$copy,
      selftransitionprob = selftransitionprob,
      maxCN = maxCN,
      sd_value = sd_value,
      usecpp = usecpp
    )
    states <- c(states, hmmresults)
  }

  dat$state <- states
  return(dat)
}

#' Copy number calling for multiple clones using Hidden Markov Models
#'
#' Main function for calling copy number states across multiple clones from
#' single-cell whole genome sequencing data. Processes each clone and chromosome
#' separately using Hidden Markov Models.
#'
#' @param dat Data frame with columns: clone_id, chr, start, end, copy.
#'   Additional columns are preserved in output.
#' @param selftransitionprob Probability of staying in the same copy number state (default: 0.99)
#' @param maxCN Maximum copy number state to consider (default: 15)
#' @param sd_value Standard deviation for emission model when clone_coverage is NULL (default: 0.2)
#' @param clone_coverage Optional data frame with clone-specific parameters.
#'   Should have columns: clone_id and either 'coverage' or 'n_cells'
#' @param usecpp Logical, whether to use C++ implementation (default: TRUE)
#' @param ncores Number of cores for parallel processing (default: 1)
#' @param progressbar Logical, whether to show progress bar (default: TRUE)
#'
#' @return Data frame with original data plus 'state' column containing
#'   inferred copy number states for each genomic bin
#'
#' @details
#' The function:
#' 1. Sorts data by clone, chromosome, and position
#' 2. Applies clone-specific standard deviations if clone_coverage provided
#' 3. Processes each clone separately (optionally in parallel)
#' 4. Fits HMM chromosome-wise to avoid unrealistic inter-chromosomal transitions
#' 5. Fills in states for bins excluded from inference (keep = FALSE)
#'
#' If a 'keep' column is present, only bins with keep = TRUE are used for
#' inference. States for excluded bins are filled based on neighboring states.
#'
#' @examples
#' # Create example data
#' library(data.table)
#' dat <- data.frame(
#'   clone_id = rep("clone1", 100),
#'   chr = rep("1", 100),
#'   start = seq(1, 1000000, 10000),
#'   end = seq(10000, 1000000, 10000),
#'   copy = c(rep(2, 50), rep(4, 30), rep(2, 20)) + rnorm(100, 0, 0.2)
#' )
#'
#' # Call copy numbers
#' results <- HMMclone(dat, selftransitionprob = 0.99, maxCN = 6)
#'
#' # With clone coverage information
#' clone_cov <- data.frame(clone_id = "clone1", coverage = 0.5)
#' results <- HMMclone(dat, clone_coverage = clone_cov)
#'
#' @export
HMMclone <- function(dat,
                     selftransitionprob = 0.99,
                     maxCN = 15,
                     sd_value = 0.2,
                     clone_coverage = NULL,
                     usecpp = TRUE,
                     ncores = 1,
                     progressbar = TRUE){

  # Make sure dataframe is in chromosome position order
  dat <- data.table::as.data.table(dat) %>%
  .[order(clone_id, chr, start)]

  if (!is.null(clone_coverage)){
    if ("n_cells" %in% names(clone_coverage)){
      clone_coverage$sd <- predict(variance_models$model_cells, clone_coverage)
      sd_clone <- clone_coverage$sd
      names(sd_clone) <- clone_coverage$clone_id
    }
    if ("coverage" %in% names(clone_coverage)){
      clone_coverage$sd <- predict(variance_models$model_coverage, clone_coverage)
      sd_clone <- clone_coverage$sd
      names(sd_clone) <- clone_coverage$clone_id
    }
    if (!all(c("coverage", "n_cells") %in% names(clone_coverage))){
      error("columns with names coverage or n_cells not present in clone_coverage data.frame!")
    }
  } else{
    #if clone_coverage is null assign sd to be the same for each clone
    clone_names <- unique(dat$clone_id)
    sd_clone <- rep(sd_value, length(clone_names))
    names(sd_clone) <- clone_names
  }

  if (!("keep" %in% colnames(dat))){
    dat$keep <- TRUE
  }

  dat_for_inference <- dplyr::filter(dat, keep == TRUE)

  if (progressbar == TRUE) {
    pb <- dplyr::progress_estimated(length(unique(dat_for_inference$clone_id)), min_time = 1)
  } else {
    pb <- NULL
  }

  if (ncores > 1) {
    dat_state <- data.table::rbindlist(parallel::mclapply(unique(dat_for_inference$clone_id),
                                                         function(clone) {
                                                           HMMcloneperchr(dplyr::filter(dat_for_inference, clone_id == clone),
                                                                              selftransitionprob = selftransitionprob,
                                                                              sd_value = sd_clone[[clone]],
                                                                              maxCN = maxCN,
                                                                              usecpp = usecpp
                                                           )
                                                         },
                                                         mc.cores = ncores
    )) %>%
      .[order(clone_id, chr, start)]
  } else {
    dat_state <- data.table::rbindlist(lapply(
      unique(dat_for_inference$clone_id),
      function(clone) {
        HMMcloneperchr(dplyr::filter(dat_for_inference, clone_id == clone),
                       selftransitionprob = selftransitionprob,
                       sd_value = sd_clone[[clone]],
                       maxCN = maxCN,
                       usecpp = usecpp
        )
      }
    )) %>%
      .[order(clone_id, chr, start)]
  }

  #fill in bins not used in inference (keep = FALSE) based on neihbouring states
  dat_state <- dplyr::full_join(dat_state,
                                dat,
                                by = setdiff(colnames(dat), "state")) %>%
    sort_by_chr() %>%
    tidyr::fill(., "state", .direction = "updown")

  return(dat_state)
}
