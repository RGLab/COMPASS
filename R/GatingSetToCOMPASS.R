## debugging
if (FALSE) {
  library(flowWorkspace)
  gs <- load_gs("../shinyGate/gs/LoveLab-Gated")
}

##' Extract A GatingSet Node to a COMPASSContainer
##' 
##' This function can be used to take gated data from a \code{flowWorkspace}
##' object (e.g., a \code{GatingSet}), and retrieve it as a 
##' \code{COMPASSContainer}, suitable for ##' modelling with \code{COMPASS}.
##' 
##' @param gs A \code{GatingSet} object.
##' @param node The name of the node. We extract all children of that node,
##'   under the assumption that the markers of interest have been gated at
##'   the chosen node.
##' @param children A character vector, or list, giving the names of the children of
##'   \code{node} that we wish to extract. If unspecified, we take all
##'   children of that node. If the character vector passed is named, we map
##'   from value to name; e.g. if a node is named \code{154} but we prefer the
##'   name \code{CD154+}, we could pass \code{`CD154+`="154"}.
##' @param meta A \code{data.frame} of meta data. If it is missing, we attempt
##'   to pull the metadata from the \code{GatingSet}, and throw an error if it
##'   is unavailable.
##' @param individual_id The name of the vector in \code{meta} that denotes the
##'   individuals from which samples were drawn. 
##' @param sample_id The name of the vector in \code{meta} that denotes the samples.
##'   This vector should contain all of the names in the \code{data} input.
##' @export
##' @seealso \code{\link{COMPASSContainer}}
GatingSetToCOMPASS <- function(gs, node, children, 
                               meta, individual_id, sample_id) {
  
  if (require(flowWorkspace)) {
    
    if (!inherits(gs, "GatingSet")) {
      stop("'gs' must be an object of class 'GatingSet'")
    }
    
    if (missing(children)) {
      children <- getChildren(gs[[1]])
    }
    
    doWarn <- FALSE
    .children <- getChildren(gs[[1]], node, isPath=TRUE)
    for (i in 1:length(gs)) {
      if (!identical(.children, getChildren(gs[[i]], node, isPath=TRUE))) {
        doWarn <- TRUE
      }
    }
    
    if (doWarn) {
      warning("The children of node '", node, "' are not identical across ",
              "each of the samples encapsulated by your GatingSet; be sure that ",
              "your samples have been gated uniformly at each node.")
    }
    
    ## get the data at each node
    data <- vector("list", length(children))
    for (i in 1:length(gs)) {
      data[[i]] <- getData(gs, children[i])
    }
    data <- getData(gs, node)
    intensities <- Map(exprs, data)
    names(intensities) <- sampleNames(gs)
    counts <- integer( length(gs) )
    for (i in seq_along(counts)) {
      counts[[i]] <- getTotal( gs[[i]], node, flowJo=FALSE )
    }
    names(counts) <- sampleNames(gs)
    
    ## try to get the metadata, if available
    if (missing(meta)) {
      tryCatch(meta <- pData(gs),
               error=function(e) {
                 stop("No metadata available!")
               })
    }
    
    CC <- COMPASSContainer(
      data=intensities,
      counts=counts,
      meta=meta,
      individual_id=individual_id,
      sample_id=sample_id
    )
    
    return(CC)
    
  } else {
    stop("This function requires 'flowWorkspace' to be installed; ",
         "try running 'biocLite(\"flowWorkspace\")' to install.")
  }
  
}


##' Create a COMPASS Container from a GatingSet
##'
##'This code expects a \code{GatingSet} or \code{GatingSetList}. 
##'It expects a regular expression for the node name (i.e. "/4\\+$" would match "/4+" in a node name with the plus
##'sign at the end of the string.
##'The user must supply the individual_id, sample_id, and stimulation_id, 
##'but they have default values suitable for the data we commonly see.
##'Sometimes the child node names don't match the marker names exactly.
##'This function will try to make some guesses about how to match these up. 
##'The \code{filter.fun} parameter is a function that does some regular expression string
##'substitution to try and clean up the node names by removing
##'various symobls that are often added to gates, {+/-\}. The user can provide their
##'own function to do string cleanup.
##'Counts are extracted as well as metadata and single cell data, and these are fed into the
##'COMPASSContainer constructor.
##'
##'@note There is likely not sufficient error checking. 
##'
##'@param gs a \code{GatingSet} or \code{GatingSetList}
##'@param node a \code{regular expression} to match a single node in the gating tree. If more than one node is matched, an error is thrown.
##'@param filter.fun a \code{function} that does string substitution to clean up node names, i.e. turns a "CD4+" into a "CD4" to try and
##'match against the \code{parameters} slot of the \code{flowFrames} in \code{gs}
##'@param individual_id a \code{character} identifying the subject id column in the \code{gs} metadata
##'@param sample_id a \code{character} idetifying the sample id column in the \code{gs} metadata.
##'@param stimulation_id a \code{character} identifying the stimulation or treatment columnin the \code{gs} metadata.
##'@param mp a \code{list} mapping node names to markers. This function tries to guess, but may fail. The user can override the guesswork.
##'@seealso \code{\link{COMPASSContainer}}
##'@examples
##'\dontrun{
##' COMPASS:::COMPASSContainerFromGatingSet(gatingset,"/4\\+$")
##'}
##'@export
COMPASSContainerFromGatingSet <- function(gs=NULL,node=NULL,filter.fun=NULL,individual_id="PTID",sample_id="name",stimulation_id="Stim",mp=NULL){
  if(require(flowWorkspace)){
  if(is.null(gs)|is.null(node)){
    stop("Must specify a gating set and parent node.")
  }
  #extract all the counts
  message("Extracting cell counts")
  stats<-getPopStats(aeras,statistic="count")
  
  pd<-pData(gs)
  #Do the expected columns exist?
  if(!all(c(sample_id,individual_id,stimulation_id)%in%colnames(pd))){
    message("Some columns not found in metadata")
    message(sprintf("Expected: %s %s %s",sample_id,individual_id,stimulation_id))
    message(sprintf("Missing: %s\n",c(sample_id,individual_id,stimulation_id)[which(!c(sample_id,individual_id,stimulation_id)%in%colnames(pd))]))
    stop("Quitting")
  }
  #Can we identify a unique parent node?
  parent.pop<-rownames(stats)[grepl(node,rownames(stats),fixed=FALSE)]
  if(length(parent.pop)>1){
    stop(sprintf("The node expression %s is not unique.",node))
  }
  if(length(parent.pop)==0){
    stop(sprintf("The node expression %s doesn't identify any nodes.",node))
  }
  
  # Grab the counts for the parent
  counts<-stats[which(rownames(stats)%in%parent.pop),]
  
  # Extract the parent node name from the full population name
  parent.node <- laply(strsplit(parent.pop,"/"),function(x)x[length(x)])
  message(sprintf("Fetching %s",parent.node))
  # Get the children of that parent and filter out boolean gates
  # Test if children exist, and test if non-empty set returned.
  message("Fetching child nodes")
  child.nodes <- getChildren(gs[[1]],parent.node)
  if(length(child.nodes)==0){
    stop(sprintf("Population %s has no children! Choose a different parent population.",parent.node))
  }
  
  child.nodes <- child.nodes[!sapply(child.nodes,function(x)flowWorkspace:::.isBoolGate(gs[[1]],x))]
  if(length(child.nodes)==0){
    stop(sprintf("All the children of %s are boolean gates. Choose a population with non-boolean child gates.",parent.node))
  }
  
  # Make sure the child node names are mapped to channel names correctly.
  # This is awful.. we don't have a way to track which dimension of a 2D gate is of importance.. so this code tries to take a guess by matching node names to marker names and doing some deduplication if there's ambiguity.
  # I cannot even begin to count the number of ways this could fail.
  # I'll check that the number of mapped nodes at the end matches the expected number of child nodes, and error out if it doesn't. We may also want to let the user pass a map.
  if(is.null(mp)){
    params<-parameters(getData(gs[[1]]))@data
    params <- data.table(params[,c("name","desc")])
    setkey(params,desc)
    if(class(filter.fun)!="function"){
      filter.fun<-function(x){
        gsub("\\\\","",gsub("/","",gsub("\\+","",x)))
      }
    }
    map <- na.omit(unique(ldply(child.nodes,function(x)params[desc%like%filter.fun(x),node:=x])))
    tbl <- table(map$node)
    if(any(tbl>1)){
      row.remove<-sapply(which(tbl>1),function(x){
        row.keep<-which(map$desc%in%filter.fun(names(tbl)[x]))
        all.row<-which(map$node%in%names(tbl)[x])
        row.remove<-setdiff(all.row,row.keep)
      })
    }
    map<-map[-c(row.remove),]
    
    #Some error checking
    if(nrow(map)!=length(child.nodes)){
      message(sprintf("We failed to guess the mapping between the node %s and the markers in the flowFrame\n",child.nodes))
      message("Our best guess was:")
      kable(map)
    }
    message("We will map the following nodes to markers:")
    kable(map)
    
    
    #construct the map
    mp<-map[,"desc"]
    names(mp)<-paste(parent.node,map[,"node"],sep="/")
    mp<-as.list(mp)
  }
  #Construct the expression
  expr<-as.name(paste(paste(parent.node,child.nodes,sep="/"),collapse="|"))
  message(sprintf("Extracting single cell data for %s",as.character(expr)))
  
  #extract the single cell values
  sc_data<-getData(obj=gs,y=expr,pop_marker_list=mp)
  message("Creating COMPASS Container")
  cc<-COMPASSContainer(data=sc_data,counts=counts,meta=pd,individual_id=individual_id,sample_id=sample_id,stimulation_id=stimulation_id)
  return(cc)
  }
  stop("This function requires flowWorkspace to be installed")
}

