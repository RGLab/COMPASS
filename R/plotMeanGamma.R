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
##' @param row_annotation A vector of names, pulled from the metadata, to be
##'   used for row annotation.
##' @param palette The colour palette to be used.
##' @param ... Optional arguments passed to \code{pheatmap}.
##' @importFrom grid grid.pretty
##' @importFrom RColorBrewer brewer.pal
plot.COMPASSResult <- function(x, y, subset, row_annotation,
  palette=brewer.pal(n=9, "Blues"), show_rownames=FALSE, 
  show_colnames=FALSE, ...) {
  
  subset_expr <- match.call()$subset
  
  if (missing(row_annotation)) {
    row_annotation <- y
  }
  
  nc <- ncol(x$fit$gamma)
  M <- x$fit$mean_gamma[, -nc]
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
  cats <- as.data.frame( lapply(cats, function(x) {
    swap(x, 0, -1)
  }))
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
  
}
