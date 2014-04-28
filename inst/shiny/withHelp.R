## Place a help icon -- the style may need to be tweaked to get it where you
## really want it on the page
helpIcon <- function(inputId, style=NULL) {
  icon <- icon("question-circle")
  icon[[2]]$attribs$id <- inputId
  if (length(style)) {
    icon[[2]]$attribs$style <- style
  }
  icon
}

## Functions that operate like the regular Shiny inputs, but also generate
## a cute little 'help' icon in the top right corner.
selectizeInputWithHelp <- function(inputId, ..., options = NULL) {
  si <- selectizeInput(inputId, ..., options = options)
  icon <- list(
    name="i",
    attribs=list(
      class="fa fa-question-circle help-icon",
      style="float: right; margin-top: 3px;",
      id=paste0(inputId, "-help")
    ),
    children=list()
  )
  class(icon) <- "shiny.tag"
  output <- c(
    list(icon),
    si
  )
  class(output) <- c("shiny.tag.list", "list")
  output
}

numericInputWithHelp <- function(inputId, label, value, min = NA, max = NA, step = NA) {
  ni <- numericInput(inputId, label, value, min, max, step)
  icon <- list(
    name="i",
    attribs=list(
      class="fa fa-question-circle help-icon",
      style="float: right; margin-top: 3px;",
      id=paste0(inputId, "-help")
    ),
    children=list()
  )
  class(icon) <- "shiny.tag"
  output <- c(
    list(icon),
    ni
  )
  class(output) <- c("shiny.tag.list", "list")
  output
}

checkboxInputWithHelp <- function(inputId, label, value = FALSE) {
  ni <- checkboxInput(inputId, label, value)
  icon <- list(
    name="i",
    attribs=list(
      class="fa fa-question-circle help-icon",
      style="float: right; margin-top: 3px;",
      id=paste0(inputId, "-help")
    ),
    children=list()
  )
  class(icon) <- "shiny.tag"
  ni[[3]][[3]] <- icon
  ni
}
