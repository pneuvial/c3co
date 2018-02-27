#' @importFrom utils capture.output str
mstr <- function(...) {
  message(paste(capture.output(str(...)), collapse = "\n"))
}

mprint <- function(...) {
  message(paste(capture.output(print(...)), collapse = "\n"))
}

mprintf <- function(...) {
  message(sprintf(...))
}
