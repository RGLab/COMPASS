## Copied from flowWorkspace
.isBoolGate <- function (x, y) {
  if (require(flowWorkspace)) {
    return(class(getGate(x, y)) == "booleanFilter")
  }
}
