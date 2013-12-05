## debugging
if (FALSE) {
  library(flowWorkspace)
  gs <- load_gs("../shinyGate/gs/LoveLab")
}

##' Create a COMPASS Container from a GatingSet
##' 
##' This code expects a \code{GatingSet} or \code{GatingSetList}. 
##' It expects a regular expression for the node name 
##' (i.e. '/4\\+$' would match '/4+' in a node name with the plus
##' sign at the end of the string.
##' The user must supply the individual_id, sample_id, and stimulation_id, 
##' but they have default values suitable for the data we commonly see.
##' Sometimes the child node names don't match the marker names exactly.
##' This function will try to make some guesses about how to match these up. 
##' The \code{filter.fun} parameter is a function that does some regular expression string
##' substitution to try and clean up the node names by removing
##' various symobls that are often added to gates, \{+/-\}. The user can provide their
##' own function to do string cleanup.
##' Counts are extracted as well as metadata and single cell data, and these are fed into the
##' COMPASSContainer constructor.
##' 
##' There is likely not sufficient error checking. 
##' 
##' @param gs a \code{GatingSet} or \code{GatingSetList}
##' @param node a \code{regular expression} to match a single node in the gating tree. If more than one node is matched, an error is thrown.
##' @param filter.fun a \code{function} that does string substitution to clean up node names, i.e. turns a 'CD4+' into a 'CD4' to try and
##' match against the \code{parameters} slot of the \code{flowFrames} in \code{gs}
##' @param individual_id a \code{character} identifying the subject id column in the \code{gs} metadata
##' @param sample_id a \code{character} idetifying the sample id column in the \code{gs} metadata.
##' @param stimulation_id a \code{character} identifying the stimulation or treatment columnin the \code{gs} metadata.
##' @param mp a \code{list} mapping node names to markers. This function tries to guess, but may fail. The user can override the guesswork.
##' @param matchmethod a \code{character} either 'regex' or 'Levenshtein' for matching nodes to markers.
##' @param markers a \code{character} vector of marker names to include.
##' @seealso \code{\link{COMPASSContainer}}
##' @examples
##' \dontrun{
##'   COMPASSContainerFromGatingSet(gatingset,'/4\\+$')
##' }
##' @importFrom plyr laply ldply
##' @importFrom knitr kable
##' @importFrom utils adist
##' @importFrom clue solve_LSAP
##' @export
COMPASSContainerFromGatingSet <- function(gs = NULL, node = NULL, filter.fun = NULL, 
                                          individual_id = "PTID", sample_id = "name", stimulation_id = "Stim", 
                                          mp = NULL, countFilterThreshold = 5000, matchmethod = c("regex", "Levenshtein"), 
                                          markers = NA) {
  if (require(flowWorkspace)) {
    
    if (is.null(gs) | is.null(node)) {
      stop("Must specify a gating set and parent node.")
    }
    
    pd <- pData(gs)
    # Do the expected columns exist?
    if (!all(c(sample_id, individual_id, stimulation_id) %in% colnames(pd))) {
      message("Some columns not found in metadata")
      message(sprintf("Expected: %s %s %s", sample_id, individual_id, 
                      stimulation_id))
      message(sprintf("Missing: %s\n", c(sample_id, individual_id, 
                                         stimulation_id)[which(!c(sample_id, individual_id, stimulation_id) %in% 
                                                                 colnames(pd))]))
      stop("Quitting")
    }
    # Can we identify a unique parent node?
    paths <- getNodes(gs[[1]], isPath=TRUE)
    parent.pop <- paths[grepl(node, paths, fixed = FALSE)]
    if (length(parent.pop) > 1) {
      stop(sprintf("The node expression %s is not unique.", node))
    }
    if (length(parent.pop) == 0) {
      stop(sprintf("The node expression %s doesn't identify any nodes.", 
                   node))
    }
    
    # extract all the counts
    message("Extracting cell counts")
    .get_cell_counts <- function(gs, samples, node) {
      nodes <- getNodes(gs[[1]], isPath=FALSE)
      paths <- getNodes(gs[[1]], isPath=TRUE)
      node <- nodes[ node == paths ]
      ind <- .Call("R_getNodeID", gs@pointer, samples[[1]], node)
      output <- unlist( lapply(samples, function(x) {
        .Call("R_getPopStats", gs@pointer, x, ind)$FlowCore["count"]
      }) )
      names(output) <- samples
      return(output)
    }
    
    if (inherits(gs, "GatingSetList")) {
      counts <- unlist( lapply(gs@data, function(x) {
        .get_cell_counts(x, sampleNames(x), parent.pop)
      }))
    } else if (inherits(gs, "GatingSet")) {
      counts <- .get_cell_counts(gs, sampleNames(gs), parent.pop)
    } else {
      stop("Internal error: 'gs' should have either been a GatingSet or a GatingSetList")
    }
    
    # Extract the parent node name from the full population name
    parent.node <- laply(strsplit(parent.pop, "/"), function(x) x[length(x)])
    message(sprintf("Fetching %s", parent.node))
    # Get the children of that parent and filter out boolean gates Test if
    # children exist, and test if non-empty set returned.
    message("Fetching child nodes")
    child.nodes <- getChildren(gs[[1]], parent.node)
    child.nodes <- basename(child.nodes)
    if (length(child.nodes) == 0) {
      stop(sprintf("Population %s has no children! Choose a different parent population.", 
                   parent.node))
    }
    
    child.nodes <- child.nodes[!sapply(child.nodes, function(x) flowWorkspace:::.isBoolGate(gs[[1]], 
                                                                                            x))]
    if (length(child.nodes) == 0) {
      stop(sprintf("All the children of %s are boolean gates. Choose a population with non-boolean child gates.", 
                   parent.node))
    }
    
    # Make sure the child node names are mapped to channel names correctly.
    # This is awful.. we don't have a way to track which dimension of a 2D
    # gate is of importance.. so this code tries to take a guess by matching
    # node names to marker names and doing some deduplication if there's
    # ambiguity.  I cannot even begin to count the number of ways this could
    # fail.  I'll check that the number of mapped nodes at the end matches
    # the expected number of child nodes, and error out if it doesn't. We
    # may also want to let the user pass a map.
    .checkMarkerConsistency <- function(xx) {
      mlist <- lapply(xx, function(x) na.omit(parameters(getData(xx, 
                                                                 use.exprs = FALSE))@data$desc))
      common <- Reduce(intersect, mlist)
      unyn <- Reduce(union, mlist)
      warnflag <- FALSE
      if (!all(common %in% unyn)) {
        warnflag <- TRUE
      }
      message("common markers are: ")
      message(sprintf("%s ", common))
      if (!is.na(markers)) {
        if (all(markers %in% common)) {
          warnflag <- FALSE
        } else {
          warnflag <- TRUE
        }
      }
      if (warnflag) {
        warning("Not all markers are shared across files.")
        message("Disparate markers are:")
        message(sprintf("%s "), setdiff(unyn, common))
      }
    }
    
    .checkMarkerConsistency(gs)
    if (is.null(mp)) {
      params <- parameters(getData(gs[[1]], use.exprs = FALSE))@data
      params <- data.table(params[, c("name", "desc")])
      # make case consistent
      params[, `:=`(desc.upper, toupper(desc))]
      child.nodes.upper <- toupper(child.nodes)
      child.nodes <- data.table(data.frame(child.nodes, child.nodes.upper))
      setkeyv(params, "desc")
      if (class(filter.fun) != "function") {
        filter.fun <- function(x) {
          gsub("-", "", gsub("\\d+\\.", "", gsub("\\\\", "", gsub("/", 
                                                                  "", gsub("\\+", "", x)))))
        }
      }
      child.nodes[, `:=`(child.nodes.upper, filter.fun(child.nodes.upper))]
      
      matchmethod <- match.arg(arg = matchmethod, choices = c("regex", 
                                                              "Levenshtein"))
      if (matchmethod == "Levenshtein") {
        distances <- adist(child.nodes[, child.nodes.upper], na.omit(params[, 
                                                                            desc.upper]))
        matching <- as.vector(solve_LSAP(distances))
        matched <- cbind(na.omit(params)[matching, list(name, desc)], 
                         node = child.nodes[, child.nodes])
        map <- matched
      } else {
        map <- na.omit(unique(ldply(child.nodes[, child.nodes.upper], 
                                    function(x) {
                                      params[desc.upper %like% x, `:=`(child.nodes.upper, 
                                                                       x)]
                                    })))
        map <- data.table(merge(map, child.nodes, by = "child.nodes.upper", 
                                all.x = TRUE))
        map[, `:=`(node, child.nodes)]
        tbl <- table(map$node)
        if (any(tbl > 1)) {
          row.remove <- sapply(which(tbl > 1), function(x) {
            row.keep <- which(map$desc %in% filter.fun(names(tbl)[x]))
            all.row <- which(map$node %in% names(tbl)[x])
            row.remove <- setdiff(all.row, row.keep)
          })
          map <- map[-c(row.remove), ]
        }
        
        # Some error checking
        if (nrow(map) != length(child.nodes[, child.nodes])) {
          message(sprintf("We failed to guess the mapping between the node %s and the markers in the flowFrame\n", 
                          child.nodes))
          message("Our best guess was:")
          kable(map)
          message("Expected nodes:")
          message(sprintf("%s ", child.nodes$child.nodes))
          message("Available dyes:")
          message(sprintf("%s ", na.omit(as.vector(params$desc))))
          message("Try specifying the mapping manually.")
          stop("Quitting")
        }
        map <- map[, c(2, 3, 6), with = FALSE]
      }
      
      # Filter based on selected markers
      if (!is.na(markers)) {
        setkey(map, desc)
        map <- map[markers, ]
      }
      message("We will map the following nodes to markers:")
      kable(map)
      
      
      # construct the map
      mp <- map[, desc]
      names(mp) <- map[, node]
      mp <- as.list(mp)
    }
    # Construct the expression
    expr <- as.name(paste(names(mp), collapse = "|"))
    message(sprintf("Extracting single cell data for %s", as.character(expr)))
    
    # extract the single cell values
    sc_data <- getData(obj = gs, y = expr, pop_marker_list = mp)
    
    message("Filtering low counts")
    filter <- counts > countFilterThreshold
    keep.names <- names(counts)[filter]
    sc_data <- sc_data[keep.names]
    counts <- counts[keep.names]
    pd <- subset(pd, eval(as.name(sample_id)) %in% keep.names)
    message(sprintf("Filtering %s samples due to low counts", length(filter) - 
                      length(keep.names)))
    
    message("Creating COMPASS Container")
    cc <- COMPASSContainer(data = sc_data, counts = counts, meta = pd, 
                           individual_id = individual_id, sample_id = sample_id, stimulation_id = stimulation_id)
    return(cc)
  }
  stop("This function requires flowWorkspace to be installed")
}

