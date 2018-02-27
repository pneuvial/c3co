#' Aligns Multiple Copy-Number Data Objects
#'
#' Expands multiple copy-number data.frame:s to have the exact same set of
#' (chr, pos) loci.  If new rows are injected, the other non-(chr, pos) fields
#' are populated with missing values.
#'
#' @param dat A list of PSCN data.frame:s with required columns \code{chr} and
#' \code{pos}.
#'
#' @return A list of the same length as \code{dat} where all data.frame:s have
#' the exact same set of fields \code{chr} and \code{pos} (and in the same
#' order).
#'
#' @export
alignLoci <- function(dat, ...) {
  stopifnot(is.list(dat))
  ## Nothing to do
  if (length(dat) <= 1L) return(dat)

  ## Expand to union of all (chr, pos):s
  chromosomes <- lapply(dat, FUN = function(df) unique(df$chr))
  chromosomes <- sort(unique(unlist(chromosomes, use.names = FALSE)), na.last = TRUE) 
  stopifnot(!anyNA(chromosomes))
  
  res <- vector("list", length = length(dat))
  for (cc in seq_along(chromosomes)) {
    chr <- chromosomes[cc]
    stopifnot(!is.na(chr))
    dat_cc <- lapply(dat, FUN = function(df) df[df$chr == chr, ])
    pos <- lapply(dat_cc, FUN = function(df) unique(df$pos))
    pos <- sort(unique(unlist(pos, use.names = FALSE)))

    dat_cc <- lapply(dat_cc, FUN = function(df) {
      idxs <- match(pos, table = df$pos)
      stopifnot(length(idxs) == length(pos))
      df <- df[idxs, ]
      stopifnot(nrow(df) == length(pos),
                all(df$chr == chr, na.rm = TRUE),
                all(df$pos == pos, na.rm = TRUE))
      ## Make sure to populate with non-missing (chromosome, x) loci
      df$chr <- chr
      df$pos <- pos
      stopifnot(!anyNA(df$chr), !anyNA(df$pos))
      stopifnot(nrow(df) == length(pos),
                all(df$chr == chr, na.rm = FALSE),
                all(df$pos == pos, na.rm = FALSE))
      df
    })

    if (cc == 1L) {
      res <- dat_cc
    } else {
      for (ii in seq_along(dat)) {
        res[[ii]] <- rbind(res[[ii]], dat_cc[[ii]])
      }
    }
  }

  ## Sanity checks
  ns <- unlist(lapply(res, FUN = nrow))
  stopifnot(all(ns == ns[1]))
  lapply(res, FUN = function(df) {
    stopifnot(!anyNA(df$chr), !anyNA(df$pos))
    chrs <- sort(unique(df$chr), na.last = TRUE)
    stopifnot(!anyNA(chrs), length(chrs) == length(chromosomes),
              all(chrs == chromosomes), !anyNA(df$pos))
  })
  
  res
}

