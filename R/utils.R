#' Sort genomic data by chromosome order
#'
#' Sorts a data frame by clone_id, chromosome (in proper genomic order), and start position.
#' Chromosomes are ordered as 1, 2, ..., 22, X, Y.
#'
#' @param dat Data frame containing genomic data with columns: clone_id, chr, start
#'
#' @return Data frame sorted by clone_id, chromosome order, and start position
#'
#' @details
#' This function ensures proper genomic ordering by converting chromosome names
#' to a standardized order (autosomes 1-22 followed by sex chromosomes X, Y).
#' This is important for genomic analyses where chromosome order matters.
#'
#' @examples
#' # Create example data with mixed chromosome order
#' dat <- data.frame(
#'   clone_id = rep("clone1", 6),
#'   chr = c("X", "2", "1", "Y", "10", "3"),
#'   start = c(1000, 2000, 3000, 4000, 5000, 6000)
#' )
#'
#' # Sort by proper chromosome order
#' sorted_dat <- sort_by_chr(dat)
#'
#' @export
sort_by_chr <- function(dat){
  dfchr <- data.frame(chr = c(paste0(1:22), "X", "Y"), idx = 1:24)
  dat <- dplyr::inner_join(dat, dfchr, by = "chr") %>%
    dplyr::arrange(clone_id, idx, start) %>%
    dplyr::select(-idx)
  return(dat)
}
