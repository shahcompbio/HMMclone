
#' Count the number of reads per bin from a file
#'
#' @param bamfile Absolute path to a bamfile
#' @param binsize binsize to use, default = "10000", use quotes!
#' @param genome Reference genome (hg19 or hg38), default = "hg19"
#' @param minMapq filter reads with mapq values below this, default = 20
#' @param isPaired
#' @param isProperPair
#' @param isUnmappedQuery
#' @param hasUnmappedMate
#' @param isMinusStrand
#' @param isMateMinusStrand
#' @param isFirstMateRead
#' @param isSecondMateRead
#' @param isSecondaryAlignment
#' @param isNotPassingQualityControls
#' @param isDuplicate
#'
#' @return Data.frame with read counts and gc and mappability per bin
#' @export
#'
#' @examples
count_reads <- function(bamfile,
                        binsize = "10000",
                        genome = "hg19",
                        minMapq = 20,
                        isPaired=TRUE,
                        isProperPair=NA,
                        isUnmappedQuery=FALSE,
                        hasUnmappedMate=NA,
                        isMinusStrand=NA,
                        isMateMinusStrand=NA,
                        isFirstMateRead=NA,
                        isSecondMateRead=NA,
                        isSecondaryAlignment=NA,
                        isNotPassingQualityControls=FALSE,
                        isDuplicate=FALSE){

  #bamfile <- "/juno/work/shah/users/william1/projects/gbm/pilot/results/bams/P1-1.bam"
  df <- references[[genome]][[paste0(binsize)]]

  gr <- GenomicRanges::GRanges(seqnames=df$chr,
                ranges=IRanges::IRanges(start=df$start, end=df$end))

  flag <- Rsamtools::scanBamFlag(isPaired=isPaired,
                      isProperPair=isProperPair,
                      isUnmappedQuery=isUnmappedQuery,
                      hasUnmappedMate=hasUnmappedMate,
                      isMinusStrand=isMinusStrand,
                      isMateMinusStrand=isMateMinusStrand,
                      isFirstMateRead=isFirstMateRead,
                      isSecondMateRead=isSecondMateRead,
                      isSecondaryAlignment=isSecondaryAlignment,
                      isNotPassingQualityControls=isNotPassingQualityControls,
                      isDuplicate=isDuplicate)

  param <- Rsamtools::ScanBamParam(flag=flag, what=c('rname', 'pos', 'mapq', 'cigar'))
  reads <- Rsamtools::scanBam(bamfile, param=param)
  reads <- reads[[1L]]

  #from QDNAseq, filter for mapq
  hasMapq <- any(is.finite(reads[['mapq']]))
  if (hasMapq) {
    keep <- which(reads[['mapq']] >= minMapq)
    reads <- lapply(reads, FUN=function(x) x[keep])
  }

  #filter for canonical chromosomes
  chroms <- unique(df$chr)
  keep <- which(reads[['rname']] %in% chroms)
  reads <- lapply(reads, FUN=function(x) x[keep])

  reads$length <- GenomicAlignments::cigarWidthAlongReferenceSpace(reads$cigar)

  ref_seq_length <- df %>%
    dplyr::group_by(chr) %>%
    dplyr::summarize(start = min(start), end = max(end)) %>%
    dplyr::summarise(x = sum(end)) %>%
    dplyr::pull(x)

  cov <- round(sum(reads$length) / ref_seq_length, 3)
  message(paste0("Genome wide coverage: ", cov))
  message(paste0("Total number of reads (millions): ", round(length(reads$rname) / 1e6, 3)))

  gr_bam <- GenomicRanges::GRanges(seqnames=reads$rname,
                    ranges=IRanges::IRanges(start=reads$pos,
                                   end=reads$pos))

  count_reads <- GenomicRanges::countOverlaps(gr, gr_bam)
  df$reads <- count_reads
  return(df)
}
