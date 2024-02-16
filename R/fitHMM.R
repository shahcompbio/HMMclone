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
