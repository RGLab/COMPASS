##' Extract Thresholded Counts from a GatingSet
##' 
##' This function extracts thresholded intensities for children of a
##' node \code{node}, as specified through the \code{map} argument.
##' 
##' \code{map} should be a list mapping node names, as specified in the
##' gating hierarchy of the gating set, to channel names, as specified
##' in either the \code{desc} or \code{name} columns of the parameters
##' of the associated \code{flowFrame}s in the \code{GatingSet}.
##' 
##' @param gs A \code{GatingSet} or \code{GatingSetList}.
##' @param node The name, or path, of a single node in a 
##'   \code{GatingSet} / \code{GatingSetList}.
##' @param map A \code{list}, mapping node names to markers.
##' @examples
##' \dontrun{
##'   GetThresholdedIntensities(gs, "4+")
##' }
##' @return A \code{list} with two components:
##' \item{\code{data:}}{A \code{list} of thresholded intensity measures.}
##' \item{\code{counts:}}{A named vector of total cell counts at the node \code{node}.}
##' @export
GetThresholdedIntensities <- function(
  gs,
  node,
  map
) {
  
  node_names <- names(map)
  channel_names <- unlist( unname(map) )
  
  if (!(inherits(gs, "GatingSet") || inherits(gs, "GatingSetList"))) {
    stop("'gs' must be either a 'GatingSet' or a 'GatingSetList'")
  }
  
  ## Turn the gs into a list of gs's
  if (inherits(gs, "GatingSet")) {
    gslist <- vector("list", length(gs))
    for (i in 1:length(gs)) {
      gslist[[i]] <- gs[[i]]
    }
    names(gslist) <- sampleNames(gs)
  } else {
    stop("Internal Error: not implemented for GatingSetList's yet.", call.=FALSE)
  }
  
  ## Make 'node' act more like a regular expression if it isn't one already
  n <- nchar(node)
  if (!substring(node, 1, 1) == "/") node <- paste0("/", node)
  if (!substring(node, n, n) == "$") node <- paste0(node, "$")
  node <- gsub("(?<!\\\\)\\+", "\\\\+", node, perl=TRUE)
  
  paths <- getNodes(gslist[[1]])
  path <- paths[ grepl(node, paths, fixed = FALSE) ]
  
  if (length(path) > 1) {
    stop(sprintf("The node expression %s is not unique.", node))
  }
  
  if (length(path) == 0) {
    stop(sprintf("The node expression %s doesn't identify any nodes.", 
      node))
  }
  
  # Extract the parent node name from the full population name
  node_name <- gsub(".*/", "", path)
  message(sprintf("Fetching data from children of '%s'.", path))
  
  # extract all the counts
  message("Extracting cell counts")
  counts <- unlist(lapply(gslist, function(x) {
    getTotal(x, path)
  }))
  
  # Get the children of that parent and filter out boolean gates Test if
  # children exist, and test if non-empty set returned.
  message("Looking for expected children nodes")
  child.nodes <- getChildren(gs[[1]], node_name)
  child.nodes <- basename(child.nodes)
  if (length(child.nodes) == 0) {
    stop(sprintf("Population %s has no children! Choose a different parent population.", 
      node_name))
  }
  
  child.nodes <- child.nodes[ !sapply(child.nodes, function(x) 
    flowWorkspace:::.isBoolGate(gs[[1]], x))]
  
  if (length(child.nodes) == 0) {
    stop(sprintf("All the children of %s are boolean gates. Choose a population with non-boolean child gates.", 
      node_name))
  }
  
  message("We will extract data from the following nodes:\n\n",
    paste( strwrap( paste(child.nodes, collapse=", "), 60), collapse="\n" ))
  
  ## This next block of code tries to figure out which channel names we will be using
  
  ## Try to guess whether we should be pulling names from the 'desc'
  ## column or the 'name' column of the flowSets
  dat <- getData(gs, use.exprs=FALSE)
  ff <- get( objects(dat@frames)[[1]], envir=dat@frames )
  params <- parameters(ff)@data
  
  ## First, check for a perfect match using a basic regex
  .check_match <- function(x, vec) {
    all( sapply(channel_names, function(x) {
      any( x == na.omit(vec) )
    }) )
  }
  
  matched_flag <- FALSE ## did we match everything in 'map' to the channel names?
  column_to_use <- NULL ## are we going to use 'desc' or 'name'?
  indices <- NULL ## what indices are we using to pull from the flowFrame?
  
  if (!all(params$name == colnames( exprs(ff) ))) {
    stop("Internal Error: expected 'params$desc' and 'colnames( exprs(ff) )' to ",
      "be identical but they are not!", call.=FALSE)
  }
  
  if (.check_match(channel_names, params$desc)) {
    message("Channel names were matched to the 'desc' column.")
    matched_flag <- TRUE
    column_to_use <- "desc"
    indices <- match(channel_names, params$desc)
  }
  
  if (.check_match(channel_names, params$name)) {
    message("Channel names were matched to the 'name' column.")
    matched_flag <- TRUE
    column_to_use <- "name"
    indices <- match(channel_names, params$name)
  }
  
  if (!matched_flag) {
    stop("Unable to match values in 'names(map)' to channel names in the data.")
  }
  
  expr_nms <- unname(params$name[indices])
  
  ## Extract the intensities
  message("Extracting cell intensities and thresholding...")
  intensities <- lapply(gslist, function(x) {
    exprs <- exprs( getData(x, path) )[, expr_nms, drop=FALSE]
    for (i in seq_along(node_names)) {
      cNode <- node_names[i]
      cChannel <- channel_names[i]
      gate <- getGate(x, file.path(path, cNode))
      thresh <- gate@min
      exprs[, expr_nms[i]][ exprs[, expr_nms[i]] < thresh ] <- 0
    }
    exprs <- exprs[ rowSums(exprs) > 0, , drop=FALSE]
    return(exprs)
  })
  
  if (!all(names(counts) == names(intensities))) {
    stop("Internal Error: expected `!all(names(counts) == names(intensities))` to be TRUE.", call.=FALSE)
  }
  
  message("Done!")
  return( list(
    data=intensities,
    counts=counts
  ) )
  
}