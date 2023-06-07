#' @export
sort_by_chr <- function(dat){
  dfchr <- data.frame(chr = c(paste0(1:22), "X", "Y"), idx = 1:24)
  dat <- dplyr::inner_join(dat, dfchr, by = "chr") %>%
    dplyr::arrange(clone_id, idx, start) %>%
    dplyr::select(-idx)
  return(dat)
}
