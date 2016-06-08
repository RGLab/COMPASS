##' Plot a COMPASSResult
##'
##' This function can be used to visualize the mean probability of response;
##' that is, the probability that there is a difference in response between
##' samples subjected to the 'treatment' condition, and samples subjected
##' to the 'control' condition.
##'
##' @aliases plot
##' @param x An object of class \code{COMPASSResult}.
##' @param y This argument gets passed to \code{row_annotation}, if
##'   \code{row_annotation} is missing. It can be used to group rows (individuals)
##'   by different conditions as defined in the metadata.
##' @param subset An \R expression, evaluated within the metadata, used to
##'   determine which individuals should be kept.
##' @param threshold A numeric threshold for filtering under-expressed
##'   categories. Any categories with mean score < \code{threshold} are
##'   removed.
##' @param minimum_dof The minimum degree of functionality for the categories
##'   to be plotted.
##' @param maximum_dof The maximum degree of functionality for the categories
##'   to be plotted.
##' @param must_express A character vector of markers that should be included
##'   in each subset plotted. For example, \code{must_express=c("TNFa & IFNg")}
##'   says we include only subsets that are positive for both
##'   \code{TNFa} or \code{IFNg}, while \code{must_express=c("TNFa", "IFNg")}
##'   says we should keep subsets which are positive for either \code{TNFa} or
##'   \code{IFNg}.
##' @param row_annotation A vector of names, pulled from the metadata, to be
##'   used for row annotation.
##' @param palette The colour palette to be used.
##' @param show_rownames Boolean; if \code{TRUE} we display row names (ie,
##'   the individual ids).
##' @param show_colnames Boolean; if \code{TRUE} we display column names
##'   (ie, the column name associated with a cytokine; typically not needed)
##' @param measure Optional. By default, we produce a heatmap of the mean
##'   gammas produced in a model fit. We can override this by supplying a
##'   matrix of suitable dimension as well; these can be generated with
##'   the \code{Posterior*} functions -- see \code{\link{Posterior}} for
##'   examples.
##' @param order_by Order rows within a group. This should be a function;
##' e.g. \code{FunctionalityScore}, \code{mean}, \code{median}, and so on.
##' Set this to \code{NULL} to preserve the original ordering of the data.
##' @param ... Optional arguments passed to \code{pheatmap}.
##' @importFrom RColorBrewer brewer.pal
##' @importFrom grDevices colorRampPalette
##' @return The plot as a \code{grid} object (\code{grob}). It can be redrawn
##' with e.g. \code{grid::grid.draw()}.
##' @examples
##' ## visualize the mean probability of reponse
##' plot(CR)
##'
##' ## visualize the proportion of cells belonging to a category
##' plot(CR, measure=PosteriorPs(CR))
##' @export
plot.COMPASSResult <- function(x, y, subset=NULL,
                               threshold=0.01,
                               minimum_dof=1,
                               maximum_dof=Inf,
                               must_express=NULL,
                               row_annotation,
                               palette=colorRampPalette(brewer.pal(10,"Purples"))(20),
                               show_rownames=FALSE,
                               show_colnames=FALSE,
                               measure=NULL,
                               order_by=FunctionalityScore,
										 markers=NULL,
                               ...) {

  call <- match.call()

  ## We do some gymnastics to figure out what the subset expression is
  ## If the caller does something like the following:
  ##
  ## subset_call <- call("%in%", a, b)
  ## plot(CR, subset=subset_call)
  ##
  ## then this function sees 'subset' as a promise that it cannot properly
  ## evaluate, and call$subset just sees 'subset_call'. So we have to explicitly
  ## evaluate the symbol in the parent frame to recover the call.
  if (!is.null(subset)) {
    subset <- call$subset
    n <- 1
    while (is.symbol(subset)) {
      subset <- eval(subset, envir=parent.frame(n))
      n <- n + 1
    }
    if (!is.language(subset)) {
      stop("'subset' should be an expression", call.=FALSE)
    }
  }
  ## Construct an object based on a reduced set of markers for plotting a heatmap
  if(!is.null(markers)){
  	if(!all(markers%in%markers(x))){
  		stop("Invalid marker names")
  	}
  	message("Computing a heatmap based on ",paste(markers,collapse=", "))
  	new_categories = unique(categories(x,FALSE)[,markers,drop=FALSE])
  	all_categories=categories(x,FALSE)[,markers,drop=FALSE]
  	dmat = as.matrix(pdist(new_categories,all_categories))
  	cat_indices = apply(dmat,1,function(y)which(y==0))
  	new_mean_gamma=sapply(cat_indices,function(i)apply(Gamma(x)[,i,],1,mean))
  	new_categories=cbind(new_categories,Counts=rowSums(new_categories))
  	reord=c(setdiff(1:nrow(new_categories),which(new_categories[,"Counts"]==0)),which(new_categories[,"Counts"]==0))
   new_categories=new_categories[reord,]
  	new_mean_gamma=new_mean_gamma[,reord]
   colnames(new_mean_gamma)=apply(new_categories[,-ncol(new_categories)],1,function(x)paste0(x,collapse=""))
  	# copy and reassign
  	X=x
  	X$fit$mean_gamma=new_mean_gamma
  	X$fit$categories = new_categories
  	x=X
  }	
  
  ## Number of markers
  .n <- ncol(x$fit$categories) - 1

  ## If order_by is missing, then the user is expecting default behavior
  ## ie, FunctionalityScore
  if (missing(order_by)) {
    order_fun <- function(x) {
      FunctionalityScore(x, n=.n)
    }
  } else {
    if (!is.null(order_by)) {
      if (is.symbol(call$order_by)) {
        if (call$order_by == "PolyfunctionalityScore") {
          ## .degree is generated downstream
          order_fun <- function(x) {
            PolyfunctionalityScore(x, degree=.degree, n=.n)
          }
        } else if (call$order_by == "FunctionalityScore") {
          order_fun <- function(x) {
            FunctionalityScore(x, n=.n)
          }
        } else {
          FUN <- match.fun(order_by)
          order_fun <- function(x) {
            apply(x, 1, FUN)
          }
        }
      } else {
        FUN <- match.fun(order_by)
        order_fun <- function(x) {
          apply(x, 1, FUN)
        }
      }
    }
  }

  ## try to override mean_gamma with measure
  if (!is.null(measure)) {
    ## normalize measure to have the same order + columns
    m <- as.data.frame(measure)
    missing_names <- colnames(x$fit$mean_gamma)[ !(colnames(x$fit$mean_gamma) %in% names(m)) ]
    m[missing_names] <- 0
    m <- m[ match(names(m), colnames(x$fit$mean_gamma)) ]
    x$fit$mean_gamma <- as.matrix(m)
  }

  ## y is effectively an alias for row_annotation, just for usability
  if (missing(row_annotation)) {
    if (missing(y)) {
      y <- NULL
    }
    row_annotation <- y
  }

  ## Rather than subsetting incrementally, we grow a set of column indices,
  ## which will finally be subsetted at the end. We begin by accepting
  ## everything, and ruling entries out based on the different
  ## subsetting criteria.

  ## Get some useful variables -- TODO: write getters for some of these
  M <- x$fit$mean_gamma
  cats <- x$fit$categories
  ind <- 1:ncol(x$fit$mean_gamma)

  ## We compute the subsets manually from the categories matrix because older
  ## COMPASS fits don't have them :/
  subsets_df <- as.data.frame(cats[, -ncol(cats), drop=FALSE])
  for (i in seq_along(subsets_df)) {
    tmp <- subsets_df[[i]]
    subsets_df[[i]] <- paste0(
      swap(tmp, c(0, 1), c("!", "")),
      colnames(subsets_df)[[i]]
    )
  }
  subsets <- do.call( function(...) paste(..., sep="&"), subsets_df )

  colnames(M) <- subsets
  dof <- cats[, "Counts"]

  if (is.null(rownames(cats))) {
    rownames(cats) <- subsets
  }

  ## Perform 'must_express' subsetting
  if (!is.null(must_express)) {
    must_express_ind <- Reduce(union, lapply(must_express, function(i) {
      which(eval(parse(text=i), envir=as.data.frame(cats)) == 1)
    }))
    ind <- intersect(ind, must_express_ind)
  }

  ## Construct the base of the 'rowann' data.frame -- annotates rows
  rowann <- data.frame(.id=rownames(M))
  
  ## the merge() below expects the COMPASS metadata x$data$meta to be a data.frame
  
  if (is(x$data$meta, "data.table")) {
    x$data$meta <- as.data.frame(x$data$meta)
  }
      
  rowann <- merge(
    rowann,
    x$data$meta[c(x$data$individual_id, row_annotation)],
    by.x=".id",
    by.y=x$data$individual_id
  )
  rowann <- rowann[!duplicated(rowann[[".id"]]), ,drop=FALSE]
  rownames(rowann) <- rowann[[".id"]]
  rowann <- rowann[-c(which(names(rowann)==".id"))]

  ## make sure M, rowann names match up
  rowann <- rowann[ match(rownames(M), rownames(rowann)), , drop=FALSE ]

  ## keep only those meeting the min, max dof criteria
  dof_ind <- which(dof >= minimum_dof & dof <= maximum_dof)
  ind <- intersect(ind, dof_ind)

  ## remove under-expressed categories
  m <- apply(M, 2, function(x) {
    mean(x, na.rm = TRUE)
  })
  keep <- m > threshold
  gone <- m <= threshold
  if (length(gone)) {
    message("The 'threshold' filter has removed ", sum(gone),
            " categories:\n", paste( colnames(M)[gone], collapse=", "))
  }
  ind <- intersect(ind, which(keep))

  ## finally, subset
  if (!length(ind)) {
    stop("no marker subsets available for plotting after subsetting")
  }

  M <- M[, ind, drop=FALSE]
  cats <- cats[ind, , drop=FALSE]
  dof <- dof[ind]

  ## handle subsetting -- this subsets rows, rather than columns, and hence
  ## lives apart of the 'ind' collection
  if (is.call(subset)) {
    keep <- unique(x$data$meta[[x$data$individual_id]][eval(subset, envir=x$data$meta)])
    M <- M[ rownames(M) %in% keep, , drop=FALSE]
    rowann <- rowann[ rownames(rowann) %in% keep, , drop=FALSE]
  }

  ## Generate .degree if using PolyfunctionalityScore
  if (is.symbol(call$order_by) && call$order_by == "PolyfunctionalityScore") {
    cats_mat <- as.matrix(cats)
    mode(cats_mat) <- "integer"
    .degree <- apply(cats_mat, 1, sum)
  }

  ## Reorder within groups based on 'order_by'
  if (!is.null(order_by)) {
    if (!is.null(row_annotation)) {
      grouping <- as.integer( factor( apply( rowann[row_annotation], 1, function(x) {
        paste(x, collapse = " ")
      } ) ) )
    } else {
      grouping <- rep(1, nrow(M))
    }
    groups <- unique(grouping)
    for (i in seq_along(groups)) {
      group <- groups[i]
      ind <- which(group == grouping)
      score <- order_fun(M[ind, , drop=FALSE])
      ord <- ind[ order(score, decreasing=TRUE) ]

      M[ind, ] <- M[ord, ]

      rownames(M)[ind] <- rownames(M)[ord]
      rownames(rowann)[ind] <- rownames(rowann)[ord]

    }
  }

  ## Group the data by row annotation
  if (!is.null(row_annotation)) {
    o <- do.call(order, as.list(rowann[row_annotation]))
    M <- M[o, , drop=FALSE]
  } else {
    rowann <- NA
  }

  ## Make sure the categories matrix is re-set as an appropriate df
  cats_df <- as.data.frame(cats[, -ncol(cats)]) ## drop the "Counts" column
  cats_df[] <- lapply(cats_df, function(x) {
    as.factor(x)
  })
  rownames(cats_df) <- colnames(M)

  ## Reorder data within degrees of functionality
  means <- apply(M, 2, mean)
  ord <- order(dof, means, decreasing = FALSE)
  M <- M[, ord, drop=FALSE]

  pheatmap(M,
           color=palette,
           show_rownames=show_rownames,
           show_colnames=show_colnames,
           row_annotation=rowann,
           cluster_rows=FALSE,
           cluster_cols=FALSE,
           cytokine_annotation=cats_df,
           ...
  )

  return( invisible(grid.grab()) )

}
