##' Fit the COMPASS Model
##' 
##' This function fits the \code{COMPASS} model.
##' 
##' @section Category Filter:
##' The category filter is used to exclude categories (combinations of
##' markers expressed for a particular cell) that are expressed very rarely.
##' It is applied to the \code{treatment} \emph{counts} matrix, which is a 
##' \code{N} samples by \code{K} categories matrix. Those categories which
##' are mostly unexpressed can be excluded here. For example, the default
##' criteria,
##' 
##' \code{category_filter=function(x) colSums(x > 5) > 2}
##' 
##' indicates that we should only retain categories for which at least two samples
##' had at least 5 cells expressing that particular combination of markers.
##' 
##' @param data An object of class \code{COMPASSData}.
##' @param treatment An \R expression, evaluated within the metadata, that
##'   returns \code{TRUE} for those samples that should belong to the
##'   treatment group. For example, if the samples that received a positive
##'   stimulation were named \code{"92TH023 Env"} within a variable in
##'   \code{meta} called \code{Stim}, you could write \code{Stim == "92TH023 Env"}.
##' @param control An \R expression, evaluated within the metadata, that
##'   returns \code{TRUE} for those samples that should belong to the
##'   control group. See above for details.
##' @param subset An expression used to subset the data. We keep only the samples
##'   for which the expression evaluates to \code{TRUE} in the metadata.
##' @param category_filter A filter for the categories that are generated. This is a
##'   function that will be applied to the \emph{treatment counts} matrix generated from
##'   the intensities. Only categories meeting the \code{category_filter} criteria will
##'   be kept.
##' @param filter_lowest_frequency A number specifying how many of the least
##'  expressed markers should be removed.
##' @param filter_specific_markers Similar to \code{filter_lowest_frequency}, 
##'   but lets you explicitly exclude markers.
##' @param model A string denoting which model to fit; \code{"discrete"}
##'   indicates we should fit the discrete \code{COMPASS} model, while
##'   \code{"continuous"} indicates we should fit the continuous
##'   \code{COMPASS} model.
##' @param iterations The number of iterations (per 'replication') to perform.
##' @param replications The number of 'replications' to perform. In order to
##'   conserve memory, we only keep the model estimates from the last replication.
##' @param verbose Boolean; if \code{TRUE} we output progress information.
##' @param ... Other arguments; currently unused.
##' 
##' @seealso \code{\link{COMPASSContainer}}, for constructing the data object
##'   required by \code{COMPASS}.
##' @return A \code{list} with class \code{COMPASSResult} with two components,
##'   the \code{fit} containing parameter estimates and parameter acceptance
##'   rates, and \code{data} containing the generated data used as input for
##'   the model.
##' @export
COMPASS <- function(data, treatment, control, subset=NULL, 
  category_filter=function(x) colSums(x > 5) > 2,
  filter_lowest_frequency=0, filter_specific_markers=NULL, 
  model=c("discrete", "continuous"), 
  iterations=40000, replications=8,
  verbose=TRUE, ...) {
  
  if (class(data) != "COMPASSContainer") {
    stop("'data' must be an object of class 'COMPASSContainer'; see the ",
      "constructor 'COMPASSContainer' for more details.")
  }
  
  ## used for brevity in later parts of code
  sid <- data$sample_id
  iid <- data$individual_id
  
  vmessage <- function(...) if (verbose) message(...) else invisible(NULL)
  
  call <- match.call()
  treatment <- call$treatment
  control <- call$control
  subset <- call$subset
  
  ## subset the data
  if (!is.null(subset)) {
    keep <- data$meta[[sid]][eval(subset, envir=data$meta)]
    vmessage("Subsetting has removed ", length(data$data) - length(keep),
      " of ", length(data$data), " samples.")
    data$data <- data$data[ names(data$data) %in% keep ]
    data$meta <- data$meta[ data$meta[[sid]] %in% keep, ]
  }
  
  .get_data <- function(data, expr, group) {
    which <- eval(expr, envir=data$meta)
    samples <- data$meta[[sid]][which]
    samples <- samples[ samples %in% names(data$data) ]
    individuals <- unique(data$meta[[iid]][ data$meta[[sid]] %in% samples ])
    individuals <- individuals[ !is.na(individuals) ]
    vmessage("There are a total of ", length(samples), " samples from ", 
      length(individuals), " individuals in the '", group, "' group.")
    
    if (length(samples) > length(individuals)) {
      vmessage("There are multiple samples per individual; ",
        "these will be aggregated.")
    }
    
    output <- vector("list", length(individuals))
    for (i in seq_along(output)) {
      indiv <- individuals[i]
      
      ## get the samples belonging to the current individual
      c_samples <- samples[ samples %in% data$meta[[sid]][ data$meta[[iid]] == indiv ] ]
      output[[i]] <- do.call(rbind, data$data[match(c_samples, names(data$data))])
    }
    names(output) <- individuals
    return(output)
  }
  
  y_s <- .get_data(data, treatment, "treatment")
  y_u <- .get_data(data, control, "control")
  
  ## make sure the names match up
  ys_names <- names(y_s)
  yu_names <- names(y_u)
  
  diff <- setdiff(
    union(ys_names, yu_names),
    intersect(ys_names, yu_names)
  )
  
  if (length(diff)) {
    vmessage("The selection criteria for 'treatment' and 'control' do not produce ",
      "paired samples for each individual. The following individual(s) will be dropped:\n\t",
      paste(diff, collapse=", "))
    keep <- intersect(ys_names, yu_names)
    ys_keep <- match(keep, ys_names)
    yu_keep <- match(keep, yu_names)
    y_s <- y_s[ys_keep]
    y_u <- y_u[yu_keep]
  }
  
  if (length(y_s) == 0 || length(y_u) == 0) {
    stop("Filtering has removed all samples.")
  }
  
  ## reorder y_u to match order of y_s
  y_u <- y_u[ match(names(y_s), names(y_u)) ]
  
  ## filter lowest frequency markers (i.e. marginalize over rarely expressed 
  ## markers so we get fewer singleton categories)
  proportions_expressed <- colMeans(do.call(rbind, y_s) > 0)
  
  ## markers we specifically want to drop
  drop_markers <- NULL
  c1 <- 0 < filter_lowest_frequency
  c2 <- filter_lowest_frequency < (length(proportions_expressed) - 2)
  if (c1 && c2) {
    drop_markers <- c(sort(proportions_expressed, decreasing = FALSE)[1:filter_lowest_frequency])
  }
  
  keep_markers <- setdiff(
    names(proportions_expressed), 
    c(names(drop_markers), filter_specific_markers)
  )
  
  ## subset based on the markers we keep
  y_s <- lapply(y_s, function(x) x[, keep_markers, drop=FALSE])
  y_u <- lapply(y_u, function(x) x[, keep_markers, drop=FALSE])
  
  ## remove the null cells
  y_s <- lapply(y_s, function(x) x[rowSums(x) > 0, , drop=FALSE])
  y_u <- lapply(y_u, function(x) x[rowSums(x) > 0, , drop=FALSE])
  
  ## generate the categories matrix here, and with it, the counts
  .generate_categories <- function(data) {
    ## i'm sorry
    tmp <- unique( as.data.table( lapply( lapply( as.data.table( do.call( rbind, data ) ), as.logical ), as.integer ) ) )
    tmp[, c("Counts") := apply(.SD, 1, sum)]
    setkeyv(tmp, c("Counts", rev(names(tmp))))
    output <- as.matrix(tmp)
    output<-output[output[,"Counts"]>0,]
    
    ## the model output requires the last row to be the 'null' category;
    ## ie, number of cells that did not express one of the combinations
    output <- rbind(output, 0)
    return(output)
  }
  
  categories <- .generate_categories( c(y_s, y_u) )
  
  .counts <- function(y, categories, counts) {
    
    ## transform categories matrix into a list suitable for the counts function
    combos <- as.list( as.data.frame( apply(categories[, -ncol(categories)], 1, function(x) {
      tmp <- c( which(x == 1), -which(x == 0) )
      tmp <- tmp[ match(1:length(tmp), abs(tmp)) ]
      return(tmp)
    })))
    
    m <- .Call("COMPASS_cell_counts", y, combos, PACKAGE="COMPASS")
    rownames(m) <- names(y)
    
    ## set the last column to be the 'null'
    m[, ncol(m)] <- counts[ names(y) ] - apply(m[,-ncol(m), drop=FALSE], 1, sum)
    
    return(m)
  }
  
  ## we have to regenerate cell count totals (by individual) to account
  ## for aggregation
  .update_total_cell_counts <- function(counts, individuals, expr) {
    which <- eval(expr, data$meta)
    new_counts <- sapply(individuals, function(ind) {
      ## get the samples corresponding to the current individual, expr
      which2 <- data$meta[[iid]] == ind
      keep <- sapply(as.logical(which * which2), isTRUE)
      smp <- as.character(data$meta[[sid]][keep])
      return( sum(counts[smp]) )
    })
    
    return(new_counts)
    
  }
  
  counts_s <- .update_total_cell_counts(data$counts, names(y_s), treatment)
  counts_u <- .update_total_cell_counts(data$counts, names(y_u), control)
  
  n_s <- .counts(y_s, categories, counts_s)
  n_u <- .counts(y_u, categories, counts_u)
  
  ## check for negative counts
  .check_negative_counts <- function(x) {
    if (any(x < 0)) stop("Internal error: negative counts in '", deparse(substitute(x)), "'")
    return( invisible(NULL) )
  }
  
  .check_negative_counts(n_s)
  .check_negative_counts(n_u)
  
  ## Check for rows in n_s, n_u with zero counts
  .zero_counts <- function(x, group) {
    bad_ids <- character()
    j <- ncol(x)
    for (i in 1:nrow(x)) {
      if (x[i, j] < 1) {
        bad_ids <- c(bad_ids, rownames(x)[i])
      }
    }
    if (length(bad_ids))
      vmessage("The following individual(s) had no cells available for ",
        "the '", group, "' group and will be removed:\n\t", 
        paste(bad_ids, sep=", "))
    
    return(bad_ids)
  }
  
  bad_ids <- c(.zero_counts(n_s, "treatment"), .zero_counts(n_u, "control"))
  if (length(bad_ids)) {
    n_s <- n_s[ !(rownames(n_s) %in% bad_ids), ]
    n_u <- n_u[ !(rownames(n_u) %in% bad_ids), ]
    counts_s <- counts_s[ !(names(counts_s) %in% bad_ids) ]
    counts_u <- counts_u[ !(names(counts_u) %in% bad_ids) ]
    y_s <- y_s[ !(names(y_s) %in% bad_ids) ]
    y_u <- y_u[ !(names(y_u) %in% bad_ids) ]
  }
  
  vmessage("The model will be run on ", length(y_s), " paired samples.")
  
  ## filter the categories matrix
  if (!is.null(category_filter)) {
    category_filter <- match.fun(category_filter)
    del <- !category_filter(n_s)
    if (is.logical(del)) {
      del <- which(del)
    }
    
    ## only do this if del actually has a length, otherwise we get a 1-column vector.
    if (length(del) > 0) { 
      
      vmessage("The category filter has removed ", length(del), " of ", nrow(categories), " categories.")
      
      ## since we're dropping some categories, those cells must be accounted 
      ## for. We add them to the negative cell category.
      n_s[, ncol(n_s)] <- n_s[,ncol(n_s)] + rowSums(n_s[,del, drop = FALSE])
      n_u[, ncol(n_u)] <- n_u[,ncol(n_u)] + rowSums(n_u[,del, drop = FALSE])
      
      ## now we drop them
      n_s <- n_s[, -c(del), drop=FALSE]
      n_u <- n_u[, -c(del), drop=FALSE]
      categories <- categories[ -c(del), , drop=FALSE ]
    } else {
      vmessage("The category filter did not remove any categories.")
    }
  }
  
  if (nrow(categories) < 2) {
    stop("There must be at least 2 categories (including the null categoy) for testing.")
  }
  
  vmessage("There are a total of ", nrow(categories), " categories to be tested.")
  
  model <- match.arg(model)
  
  ## go to the model fitting processes
  switch(model,
    discrete={
      vmessage("Fitting discrete COMPASS model.")
      output <- list(
        fit=.COMPASS.discrete(n_s=n_s, n_u=n_u, categories=categories,
          iterations=iterations, replications=replications, verbose=verbose, ...),
        data=list(n_s=n_s, n_u=n_u, counts_s=counts_s, counts_u=counts_u,
          categories=categories, meta=data$meta, sample_id=data$sample_id,
          individual_id=data$individual_id),
        orig=data
      )
    },
    continuous={
      vmessage("Fitting continuous COMPASS model.")
      output <- list(
        fit=.COMPASS.continuous(y_s=y_s, y_u=y_u, n_s=n_s, n_u=n_u,
          categories=categories, 
          iterations=iterations, replications=replications, verbose=verbose, ...),
        data=list(y_s=y_s, y_u=y_u, n_s=n_s, n_u=n_u, 
          counts_s=counts_s, counts_u=counts_u, categories=categories, 
          meta=data$meta, sample_id=data$sample_id, individual_id=data$individual_id),
        orig=data
      )
    }
  )
  
  class(output) <- "COMPASSResult"
  return(output)
  
}
