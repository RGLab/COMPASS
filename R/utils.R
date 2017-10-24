stopf <- function(...) stop(gettextf(...), call.=FALSE)
warnf <- function(...) warning(gettextf(...), call.=FALSE)

clamp <- function(x, from, to) {
  x[ x < from ] <- from
  x[ x > to ] <- to
  x
}

##' Markers
##'
##' Returns the markers associated with an experiment.
##'
##' @param object An \R object.
##' @export
markers <- function(object) {
  UseMethod("markers")
}

##' @export
markers.COMPASSContainer <- function(object) {
  colnames(object$data[[1]])
}

##' @export
markers.COMPASSResult <- function(object) {
  cats <- categories(object, FALSE)
  colnames(cats)
}

##' Categories
##'
##' Returns the categories matrix in a \code{COMPASSResult} object.
##'
##' @param x A \code{COMPASSResult} object.
##' @param counts Boolean; if \code{TRUE} we return the counts (degree of
##'   functionality) as well.
##' @export
categories <- function(x, counts) {
  if (!inherits(x, "COMPASSResult")) {
    stop("'x' is not a 'COMPASSResult'", call.=FALSE)
  }
  cats <- x$data$categories
  if (!counts) {
    cats <- cats[, colnames(cats) != "Counts", drop=FALSE]
  }
  return(cats)
}

##' Metadata
##'
##' Functions for getting and setting the metadata associated with an object.
##'
##' @param x An \R object.
##' @param value An \R object appropriate for storing metadata in object \code{x};
##'   typically a \code{data.frame}.
##' @export
metadata <- function(x) UseMethod("metadata")

##' @export
#' @rdname metadata
metadata.COMPASSContainer <- function(x) {
  x$meta
}

##' @export
#' @rdname metadata
metadata.COMPASSResult <- function(x) {
  x$data$meta
}

##' @export
#' @rdname metadata
`metadata<-` <- function(x, value) UseMethod("metadata<-")

##' @export
#' @rdname metadata
`metadata<-.COMPASSContainer` <- function(x, value) {
  x$meta <- check_meta(x$data, x$counts, value, x$individual_id, x$sample_id)
  x
}

##' COMPASSResult Accessors
##'
##' These functions can be used for accessing data within a \code{COMPASSResult}.
##'
##' @param x A \code{COMPASSResult} object.
##' @name COMPASSResult-accessors
NULL

##' The gamma array associated with a \code{COMPASS} model fit.
##' @rdname COMPASSResult-accessors
Gamma <- function(x) {
  return(x$fit$gamma)
}

##' @rdname COMPASSResult-accessors
MeanGamma <- function(x) {
  return(x$fit$mean_gamma)
}


#' Flag COMPASS boolean populations
#'
#' Returns a boolean vector indexing cell populations in \code{cellpops} that match
#' the pattern for boolean combinations of \code{markers}.
#'
#' @details
#' If markers A, B, C, D make up the population names in \code{cellpops} and they
#'the names match the pattern e.g. "A+B-C+D+,Count" (typical of exports from some gating tools),
#'then \code{markers} should be a vector of markers in the same order they appear in \code{cellpops}.
#'
#' @param cellpops \code{vector} of character names of cell populations.
#' @param markers \code{vector} of character names of markers in the order they appear in the population names.
#'
#' @return A boolean vector indexing \code{cellpops} with \code{TRUE} for populations matchin
#' the pattern.
#'
#' @export
#' @seealso translate_marker_names
#' @examples
#' #Generate some population names
#' markers = LETTERS[1:4]
#' pos = c("+","-")
#' popnames = apply(expand.grid(pos,pos,pos,pos),1,
#'             function(x)paste(paste(paste(markers,x,sep=""),
#'             collapse=""),",Count",sep=""))
#' popnames = sample(c(popnames,paste(paste(markers,sample(c("+","-"),
#'              length(markers),replace=TRUE),sep=""),",Count",sep="")))
#' popnames[select_compass_pops(popnames,LETTERS[1:4])]
select_compass_pops = function(cellpops,markers){
  pattern = paste0(paste0(paste0(markers,"[+-]"),collapse=""),",Count$",sep="")
  grepl(pattern,cellpops)
}

#' Translate marker names to format use by COMPASS
#'
#' Translate boolean population names from format exported by common
#' software tools to a format used by COMPASS.
#'
#' @param cellpops \code{character} vector of cell population names.
#'
#' @return \code{character} vector of cell population names used by COMPASS
#' @export
#' @seealso select_compass_pops
#' @examples
#' #Generate marker names
#' markers = LETTERS[1:4]
#' pos = c("+","-")
#' popnames = apply(expand.grid(pos,pos,pos,pos),1,
#'               function(x) paste(paste(paste(markers,x,sep=""),
#'               collapse=""),",Count",sep=""))
#' popnames = sample(c(popnames,
#'            paste(paste(markers,sample(c("+","-"),
#'            length(markers),replace=TRUE),sep=""),
#'            ",Count",sep="")))
#' popnames = popnames[select_compass_pops(popnames,LETTERS[1:4])]
#' #Translate
#' translate_marker_names(popnames)
translate_marker_names = function(cellpops){
  gsub(",Count$","",gsub("&$","",gsub("(\\w*)\\+","\\1&",gsub("(\\w*)-","!\\1&",cellpops))))
}


#' @name COMPASSfitToCountsTable
#' @title Extract a table of counts from a COMPASSResult object
#' @description Returns a table of counts and parent counts for each cell subset in a COMPASS fit.
#' @param x \code{COMPASSResult}
#' @param idcol unquote variable name in the metadata for the subject id.
#' @param stimName the name of the stimulation
#' @param drop numeric vector indicating the columns in the metadata to drop from the output. Usually sample-specific columns rather than subject specific columns.
#' @param parent character name of the parent population for this model fit. e.g. "CD4"
#' @importFrom dplyr left_join
#' @importFrom tidyr gather 
#' @importFrom tidyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr bind_rows
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom dplyr enquo
#' @importFrom rlang !!
#' @export
COMPASSfitToCountsTable = function(x,idcol=NULL,parent=NULL,drop=NULL, stimName=NULL){
    stopifnot(inherits(x,"COMPASSResult"))
    ns = x$data$n_s
    ns = data.frame(id = rownames(ns),ns,check.names = FALSE)
    nu = x$data$n_u
    nu = data.frame(id = rownames(nu), nu,check.names = FALSE)
    cs = x$data$counts_s
    cs = data.frame(id = names(cs), ParentCount = cs,check.names = FALSE)
    cu = x$data$counts_u
    cu = data.frame(id = names(cu),ParentCount = cu,check.names = FALSE)
    ns = gather(ns,population,Count,-1L)
    nu = gather(nu,population,Count,-1L)
    stim = left_join(ns,cs,by = "id")
    unstim = left_join(nu,cu,by = "id")
    meta = x$data$meta
    stim = stim %>% mutate(Stim = stimName)
    unstim = unstim %>% mutate(Stim = "UNS")
    data = bind_rows(stim,unstim) %>% mutate(id = as.character(id))
    if(is.null(match.call()$parent)){
      stop("argument 'parent' is missing. 'parent' is a character string describing the cell population parent subset: .e.g. \"cd4\"",call. =FALSE)
    }
    idcol = enquo(idcol)
    parent = enquo(parent)
    meta = meta %>% rename(id = !!idcol) %>% select(-!!drop) %>% mutate(id = as.character(id))
    data = left_join(data,meta, by ="id")
    return(data)
}
