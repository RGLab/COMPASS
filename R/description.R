##' Get and Set the Description for the Shiny Application
##' 
##' This is used for setting an informative description used in the Shiny
##' application.
##' 
##' Information about the \code{COMPASS} results will be auto-generated.
##' 
##' @param x A \code{COMPASS} fit.
##' @param value A set of paragraphs describing the experiment, as a character
##'   vector.
##' @export
COMPASSDescription <- function(x) {
  x$description
}

##' @rdname COMPASSDescription
##' @export
"COMPASSDescription<-" <- function(x, value) {
  x$description <- value
  return(x)
}
