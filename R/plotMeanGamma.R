##' Plot a COMPASSResult
##' 
##' This function can be used to visualize the mean probability of response --
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
##' @param remove_unexpressed_categories Boolean, if \code{TRUE} we remove
##'   any unexpressed categories.
##' @param minimum_dof The minimum degree of functionality for the categories
##'   to be plotted.
##' @param maximum_dof The maximum degree of functionality for the categories
##'   to be plotted.
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
##' @examples
##' ## visualize the mean probability of reponse
##' plot(CR)
##' 
##' ## visualize the proportion of cells belonging to a category
##' plot(CR, measure=PosteriorPs(CR))
##' @export
plot.COMPASSResult <- function(x, y, subset, 
  remove_unexpressed_categories=TRUE, 
  minimum_dof=1, 
  maximum_dof=Inf, 
  row_annotation,
  #palette=seq_gradient_pal(low="black", high="red")(seq(0, 1, length=20)),
  palette=colorRampPalette(brewer.pal(10,"Purples"))(20),
  show_rownames=FALSE, 
  show_colnames=FALSE, 
  measure=NULL,
  order_by=FunctionalityScore,
  ...) {
  
  call <- match.call()
  
  ## If order_by is missing, then the user is expecting default behavior
  ## ie, FunctionalityScore
  if (missing(order_by)) {
    order_fun <- FunctionalityScore
  } else {
    if (!is.null(order_by)) {
      
      if (is.symbol(call$order_by)) {
        if (call$order_by == "PolyfunctionalityScore") {
          ## .degree is generated downstream
          order_fun <- function(x) {
            PolyfunctionalityScore(x, degree=.degree)
          }
        } else if (call$order_by == "FunctionalityScore") {
          order_fun <- FunctionalityScore
        } else {
          FUN <- match.fun(order_by)
          order_fun <- function(x, f=FUN) {
            apply(x, 1, f)
          }
        }
      } else {
        FUN <- match.fun(order_by)
        order_fun <- function(x, f=FUN) {
          apply(x, 1, f)
        }
      }
    }
  }
  
  ## try to override mean_gamma with measure
  if (!is.null(measure)) {
    x$fit$mean_gamma <- measure
  }
  
  subset_expr <- match.call()$subset
  
  if (missing(row_annotation)) {
    if (missing(y)) {
      y <- NULL
    } 
      row_annotation <- y
  }
  
  nc <- ncol(x$fit$gamma)
  M <- x$fit$mean_gamma[, -nc]
  
  ## get the dof
  dof <- x$fit$categories[, "Counts"]
  dof <- dof[ -length(dof) ]

  rowann <- data.frame(.id=rownames(M))
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
 
  cats <- x$fit$categories[-nc,]
  cats <- data.frame(cats)
  cats <- cats[,1:(ncol(cats)-1)]
#   cats <- as.data.frame( lapply(cats, function(x) {
#     swap(x, 0, -1)
#   }))
  cats <- as.data.frame( lapply(cats, function(x) {
    factor(x, levels=c(0, 1))
  }))
  
  ## keep only those meeting the min, max dof criteria
  M <- M[, dof >= minimum_dof & dof <= maximum_dof, drop=FALSE]
  cats <- cats[dof >= minimum_dof & dof <= maximum_dof, ]
  if (ncol(M) == 0) {
    stop("No categories left after subsetting for 'minimum_dof', 'maximum_dof'")
  }
  
  ## remove unexpressed categories
  if (remove_unexpressed_categories) {
    m <- apply(M, 2, sum)
    keep <- m != 0
    M <- M[, keep, drop=FALSE]
    cats <- cats[keep, ]
  }
  
  colnames(M) <- rownames(cats)
  
  ## handle subsetting
  if (!missing(subset)) {
    keep <- x$data$meta[[x$data$individual_id]][eval(subset_expr, envir=x$data$meta)]
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
      grouping <- as.integer( factor( apply( rowann[row_annotation], 1, paste ) ) )
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
  
  pheatmap(M,
    color=palette,
    show_rownames=show_rownames,
    show_colnames=show_colnames,
    row_annotation=rowann,
    cluster_rows=FALSE,
    cluster_cols=FALSE,
    cytokine_annotation=cats,
    ...
  )
  
  return (invisible(M))
  
}
