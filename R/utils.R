#' @importFrom utils capture.output str
mstr <- function(...) {
  message(paste(capture.output(str(...)), collapse = "\n"))
}
