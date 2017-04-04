#' @importFrom utils str
mstr <- function(...) {
  message(paste(capture.output(str(...)), collapse = "\n"))
}
