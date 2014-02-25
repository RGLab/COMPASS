.needsFlowWorkspace <- function() {
  stop("This function requires 'flowWorkspace' -- try biocLite('flowWorkspace') to install.",
    call.=FALSE)
}

## Copied from flowWorkspace
.isBoolGate <- function (x, y) {
  if (require(flowWorkspace)) {
    return(class(flowWorkspce::getGate(x, y)) == "booleanFilter")
  } else {
    .needsFlowWorkspace()
  }
}
