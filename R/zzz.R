# gtxversion() output when attaching package
.onAttach <- function(...) {
  packageStartupMessage(gtxversion())
}