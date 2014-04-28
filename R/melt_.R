##' Make a 'Wide' data set 'Long'
##'
##' Inspired by \code{reshape2:::melt}, we melt \code{data.frame}s and
##' \code{matrix}s. This function is built for speed.
##'
##' If items to be stacked are not of the same internal type, they will be
##' promoted in the order \code{logical} > \code{integer} > \code{numeric} >
##' \code{character}.
##'
##' @param data The \code{data.frame} to melt.
##' @param ... Arguments passed to other methods.
melt_ <- function(data, ...) {
  UseMethod("melt_")
}

##' @rdname melt_
##' @param id.vars Vector of id variables. Can be integer (variable position)
##'  or string (variable name). If blank, we use all variables not in \code{measure.vars}.
##' @param measure.vars Vector of measured variables. Can be integer (variable position)
##'  or string (variable name). If blank, we use all variables not in \code{id.vars}.
##' @param variable.name Name of variable used to store measured variable names.
##' @param value.name Name of variable used to store values.
##' @export
melt_.data.frame <- function(data, id.vars, measure.vars, variable.name="variable", ..., value.name="value") {

  "%nin%" <- function(x, y) {
    return( !(x %in% y) )
  }

  any_na <- function(x) {
    return( any(is.na(x)) )
  }

  ## figure out which variables belong to id.vars, measure.vars,
  if( missing(measure.vars) ) {
    if( missing(id.vars) ) {
      ## assume that we are going to melt everything
      id.vars <- integer(0)
    }
    if( is.character(id.vars) ) {
      id.vars <- match( id.vars, names(data) )
    }
    measure.vars <- which( 1:length(data) %nin% id.vars )
  }

  if( missing(id.vars) ) {
    if( missing(measure.vars) ) {
      stop("if 'id.vars' is missing, you must supply 'measure.vars'")
    }
    if( is.character(measure.vars) ) {
      measure.vars <- match( measure.vars, names(data) )
    }
    id.vars <- which( 1:length(data) %nin% measure.vars )
  }

  if (is.character(id.vars)) {
    id.vars <- match(id.vars, names(data))
  }

  if (is.character(measure.vars)) {
    measure.vars <- match(measure.vars, names(data))
  }

  if (is.null(id.vars)) {
    id.vars <- integer(0)
  }

  if (any_na(id.vars)) {
    stop("Failed to match all of 'id.vars' to variable names in 'data'")
  }

  if (any_na(measure.vars)) {
    stop("Failed to match all of 'measure.vars' to variable names in 'data'")
  }

  if (any(id.vars < 1) || any(id.vars > length(data))) {
    stop("one or more of the 'id.vars' indexes beyond column range of data")
  }

  if (any(measure.vars < 1) || any(measure.vars > length(data))) {
    stop("one or more of the 'measure.vars' indexes beyond column range of data")
  }

  return( .Call( C_melt_dataframe,
    data,
    as.integer(id.vars-1L),
    as.integer(measure.vars-1L),
    variable.name,
    value.name
  ) )

}

##' @rdname melt_
##' @export
melt_.matrix <- function( data, ... ) {
  return( .Call(C_melt_matrix, data) )
}
