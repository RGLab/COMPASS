.needsFlowWorkspace <- function() {
  stop("This function requires 'flowWorkspace' -- try biocLite('flowWorkspace') to install.",
    call.=FALSE)
}

## Copied from flowWorkspace
.isBoolGate <- function (x, y) {
  if (requireNamespace("flowWorkspace",quietly=TRUE)) {
    return(class(flowWorkspace::getGate(x, y)) == "booleanFilter")
  } else {
    .needsFlowWorkspace()
  }
}
