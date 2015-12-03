##' Create a COMPASS Container from a GatingSet
##'
##' This code expects a \code{GatingSet} or \code{GatingSetList}.
##' It expects a regular expression for the node name
##' (i.e. '/4\\+$' would match '/4+' in a node name with the plus
##' sign at the end of the string. Alternatively, you can supply a
##' partial path.
##' The user must supply the individual_id and sample_id,
##' but they have default values suitable for the data we commonly see.
##' Sometimes the child node names don't match the marker names exactly.
##' This function will try to make some guesses about how to match these up.
##' The \code{filter.fun} parameter is a function that does some regular expression string
##' substitution to try and clean up the node names by removing
##' various symobls that are often added to gates, \{+/-\}. The user can provide their
##' own function to do string cleanup.
##' Counts are extracted as well as metadata and single cell data, and these are fed into the
##' \code{COMPASSContainer} constructor.
##'
##' There is likely not sufficient error checking.
##'
##' @param gs a \code{GatingSet} or \code{GatingSetList}
##' @param node a \code{regular expression} to match a single node in the gating tree. If more than one node is matched, an error is thrown.
##' @param filter.fun a \code{function} that does string substitution to clean up node names, i.e. turns a 'CD4+' into a 'CD4' to try and
##' match against the \code{parameters} slot of the \code{flowFrames} in \code{gs}
##' @param individual_id a \code{character} identifying the subject id column in the \code{gs} metadata
##' @param sample_id a \code{character} idetifying the sample id column in the \code{gs} metadata.
##' @param mp a \code{list} mapping node names to markers. This function tries to guess, but may fail. The user can override the guesswork.
##' @param matchmethod a \code{character} either 'regex' or 'Levenshtein' for matching nodes to markers.
##' @param markers a \code{character} vector of marker names to include.
##' @param swap a \code{logical} default FALSE. Set to TRUE if the marker and channel names are swapped.
##' @param countFilterThreshold \code{numeric} threshold. if the number of cells expressing at
##'   least one marker of interest is less than this threshold, we remove that
##'   file. Default is 5000.
##' @seealso \code{\link{COMPASSContainer}}
##' @examples \dontrun{
##' ## gs is a GatingSet from flowWorkspace
##' COMPASSContainerFromGatingSet(gs, "4+")
##' }
##' @importFrom plyr laply ldply
##' @importFrom knitr kable
##' @importFrom utils adist
##' @importFrom clue solve_LSAP
##' @export
COMPASSContainerFromGatingSet<-function(gs = NULL, node = NULL, filter.fun = NULL,
                                        individual_id = "PTID", sample_id = "name",
                                        mp = NULL,
                                        matchmethod = c("Levenshtein","regex"),
                                        markers = NA,swap=FALSE, countFilterThreshold = 5000) {
  if (requireNamespace("flowWorkspace",quietly = TRUE)) {

    ## R CMD check silencing
    desc.upper <- desc <- name <- NULL

    if (is.null(gs) | is.null(node)) {
      stop("Must specify a gating set and parent node.")
    }
    unique.node<-node
    ## Make 'node' act more like a regular expression if it isn't one already
    n <- nchar(node)
    if (!substring(node, 1, 1) == "/") node <- paste0("/", node)
    if (!substring(node, n, n) == "$") node <- paste0(node, "$")
    node <- gsub("(?<!\\\\)\\+", "\\\\+", node, perl=TRUE)

    # extract all the counts
    message("Extracting cell counts")
    .getOneStat<-function(x,y){
      parent.counts<-flowWorkspace::lapply(x,function(xx,yy=y){
        flowWorkspace::getTotal(xx,yy)
      })
      parent.counts <- unlist(parent.counts)
      names(parent.counts) <- flowWorkspace::sampleNames(x)
      parent.counts
    }

    nnames <- flowWorkspace::getNodes(gs[[1]], path="full")
    parent.pop<-nnames[grepl(node, nnames, fixed = FALSE)]
    if (length(parent.pop) > 1) {
      stop(gettextf("The node expression %s is not unique.", node))
    }
    if (length(parent.pop) == 0) {
      stop(gettextf("The node expression %s doesn't identify any nodes.",
                    node))
    }

    # Extract the parent node name from the full population name
    # we can just use the parent.pop
    parent.node <- laply(strsplit(parent.pop, "/"), function(x) x[length(x)])
    message(gettextf("Fetching %s", parent.node))

    counts<-.getOneStat(gs,unique.node)

    #stats <- getPopStats(gs, statistic = "count")

    pd <- pData(gs)
    # Do the expected columns exist?
    if (!all(c(sample_id, individual_id) %in% colnames(pd))) {
      message("Some columns not found in metadata")
      message(gettextf("Expected: %s %s", sample_id, individual_id))
      message(gettextf("Missing: %s\n", c(sample_id, individual_id)[which(!c(sample_id, individual_id) %in%
                                                                            colnames(pd))]))
      stop("Quitting")
    }
    #validity check for name column (since at flowSet level, flowCore now allows name column to be different from row.names)
    if(!isTRUE(all.equal(as.vector(rownames(pd)), as.vector(pd[[sample_id]]))))
      stop("sample names are not consistent with rownames of pData!")

    # Get the children of that parent and filter out boolean gates Test if
    # children exist, and test if non-empty set returned.
    message("Fetching child nodes")
    full.child.nodes<-flowWorkspace::getChildren(gs[[1]], unique.node,path="auto")
    child.nodes <- basename(flowWorkspace::getChildren(gs[[1]], unique.node))

    if (length(child.nodes) == 0) {
      stop(gettextf("Population %s has no children! Choose a different parent population.",
                    parent.node))
    }

    child.nodes <- child.nodes[!sapply(full.child.nodes, function(x) .isBoolGate(gs[[1]],
                                                                                 x))]
    full.child.nodes <- full.child.nodes[!sapply(full.child.nodes, function(x) .isBoolGate(gs[[1]],x))]

    if (length(child.nodes) == 0) {
      stop(gettextf("All the children of %s are boolean gates. Choose a population with non-boolean child gates.",
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
      if (inherits(xx, "GatingSetList")) {
        mlist <- unlist( recursive=FALSE, lapply(xx@data, function(x) {
          dat <- flowWorkspace::getData(x, use.exprs=FALSE)
          lapply( objects(dat@frames), function(obj) {
            fr <- get(obj, envir=dat@frames)
            return(na.omit( parameters(fr)@data$desc ))
          })
        }) )
      } else if (inherits(xx, "GatingSet")) {
        dat <- flowWorkspace::getData(xx, use.exprs=FALSE)
        mlist <- lapply( objects(dat@frames), function(obj) {
          fr <- get(obj, envir=dat@frames)
          return(na.omit( parameters(fr)@data$desc ))
        })
      } else {
        stop("Expected object of type 'GatingSetList' or 'GatingSet'")
      }
      mlist <- flowWorkspace::lapply(xx, function(x) na.omit(parameters(getData(x,
                                                                                use.exprs = FALSE))@data$desc))
      common <- Reduce(intersect, mlist)
      unyn <- Reduce(union, mlist)
      warnflag <- FALSE
      if (!all(common %in% unyn)) {
        warnflag <- TRUE
      }
      message("common markers are: ")
      message(gettextf("%s ", common))
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
        message(gettextf("%s ", setdiff(unyn, common)))
      }
    }

    .checkMarkerConsistency(gs)
    if (is.null(mp)) {
      params <- parameters(flowWorkspace::getData(gs[[1]], use.exprs = FALSE))@data
      params <- data.table(params[, c("name", "desc")])
      if(swap){
        setnames(params,c("name","desc"),c("desc","name"))
      }
      # make case consistent
      params[, `:=`(desc.upper, toupper(desc))]
      child.nodes.upper <- toupper(child.nodes)
      child.nodes <- data.table(data.frame(child.nodes, child.nodes.upper,full.child.nodes))
      setkeyv(params, "desc")
      if (class(filter.fun) != "function") {
        filter.fun <- function(x) {
          gsub("-", "", gsub("\\d+\\.", "", gsub("\\\\", "", gsub("/",
                                                                  "", gsub("\\+", "", x)))))
        }
      }
      child.nodes[, `:=`(child.nodes.upper, filter.fun(child.nodes.upper))]

      matchmethod <- match.arg(arg = matchmethod, choices = c("Levenshtein",
                                                              "regex"))

      if (matchmethod == "Levenshtein") {
        distances <- adist(child.nodes[, child.nodes.upper], na.omit(params[,
                                                                            desc.upper]))
        matching <- as.vector(solve_LSAP(distances))
        matched <- cbind(as.data.frame(params[matching, list(name, desc)]),  data.frame(node=child.nodes[, full.child.nodes]))
        map <- data.table(matched)
      } else {
        if(swap){
          stop("matchmethod regex not supported when swap=TRUE");
        }
        map <- na.omit(unique(ldply(child.nodes[, child.nodes.upper],
                                    function(x) {
                                      params[desc.upper %like% x, `:=`(child.nodes.upper,
                                                                       x)]
                                    })))
        map <- data.table(merge(map, child.nodes, by = "child.nodes.upper",
                                all.x = TRUE))
        map[, `:=`(node, child.nodes)]
        #drop uncompensated markers
        map<-map[name%like%"<"]
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
          message(gettextf("We failed to guess the mapping between the node %s and the markers in the flowFrame\n",
                           child.nodes))
          message("Our best guess was:")
          kable(map)
          message("Expected nodes:")
          message(gettextf("%s ", child.nodes$child.nodes))
          message("Available dyes:")
          message(gettextf("%s ", na.omit(as.vector(params$desc))))
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
      kable(as.data.frame(map))


      # construct the map
      #       if(swap){
      #         mp <- as.character(map[,name])
      #       }else{
      mp <- as.character(map[, desc])
      #       }
      names(mp) <- map[, node]
      mp <- as.list(mp)
    }
    # Construct the expression
    expr <- as.name(paste(names(mp), collapse = "|"))
    message(gettextf("Extracting single cell data for %s", as.character(expr)))

    # extract the single cell values
    #exprs can now be a vector of characters
    expr<-do.call(c,strsplit(as.character(expr),"\\|"))
    sc_data <- try(getSingleCellExpression( x=gs, nodes = expr, map = mp,swap=swap))
    if(inherits(sc_data,"try-error")){
      message("getData failed. Perhaps the marker list is not unique in the flowFrame.")
      message("All markers and channels:")
      kable(na.omit(data.frame(params[,1:2,with=FALSE])))
      stop()
    }


    message("Creating COMPASS Container")
    cc <- COMPASSContainer(data = sc_data, counts = counts, meta = pd,
                           individual_id = individual_id, sample_id = sample_id, countFilterThreshold = countFilterThreshold)
    return(cc)
  }

  .needsFlowWorkspace()

}
