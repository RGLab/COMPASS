##' Plot a Heatmap of the Mean Gammas
##' 
##' This function plots the mean gammas.
##' 
##' @method plot COMPASSResult
##' @S3method plot COMPASSResult
##' @param x An object of class \code{COMPASSResult}.
##' @param y This argument gets passed to \code{row_annotation}, if
##'   \code{row_annotation} is missing.
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
##' @param ... Optional arguments passed to \code{pheatmap}.
##' @importFrom scales seq_gradient_pal
plot.COMPASSResult <- function(x, y, subset, 
  remove_unexpressed_categories=TRUE, minimum_dof=1, maximum_dof=Inf, 
  row_annotation,
  palette=seq_gradient_pal(low="black", high="red")(seq(0, 1, length=20)),
  show_rownames=FALSE, 
  show_colnames=FALSE, ...) {
  
  subset_expr <- match.call()$subset
  
  if (missing(row_annotation)) {
    row_annotation <- y
  }
  
  nc <- ncol(x$fit$gamma)
  M <- x$fit$mean_gamma[, -nc]
  
  ## compute dof from the colnames of M
  dof <- sapply( strsplit( colnames(M), "", fixed=TRUE ), function(x) {
    sum( as.integer(x) )
  })
  
  rowann <- data.frame(.id=rownames(M))
  rowann <- merge(
    rowann, 
    x$data$meta[c(x$data$individual_id, row_annotation)], 
    by.x=".id",
    by.y=x$data$individual_id
  )
  rowann <- rowann[!duplicated(rowann[[".id"]]), ]
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
  
  ## reorder the data
  o <- do.call(order, as.list(rowann[row_annotation]))
  
  pheatmap(M[o,],
    color=palette,
    show_rownames=show_rownames,
    show_colnames=show_colnames,
    row_annotation=rowann,
    cluster_rows=FALSE,
    cluster_cols=FALSE,
    cytokine_annotation=cats,
    ...
  )
  
  return (invisible(M[o, ]))
  
}
