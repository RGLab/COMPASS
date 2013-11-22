##' Plot a Heatmap of the Mean Gammas
##' 
##' This function plots the mean gammas.
##' 
##' @param x An object of class \code{COMPASSResult}.
##' @param y An object of class \code{COMPASSResult}.
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
##' @importFrom scales div_gradient_pal
##' @export
plot2 <- function(x, y, subset, 
  remove_unexpressed_categories=TRUE, minimum_dof=1, maximum_dof=Inf, 
  row_annotation, 
  palette=div_gradient_pal(low="blue", mid="black", high="red")(seq(0, 1, length=20)),
  show_rownames=FALSE, 
  show_colnames=FALSE, ...) {
  
  subset_expr <- match.call()$subset
  
  nc_x <- ncol(x$fit$gamma)
  M_x <- x$fit$mean_gamma[, -nc_x]
  
  nc_y <- ncol(y$fit$gamma)
  M_y <- y$fit$mean_gamma[, -nc_y]
  
  ## make sure the row order is the same
  M_x <- M_x[ order(rownames(M_x)), , drop=FALSE ]
  M_y <- M_y[ order(rownames(M_y)), , drop=FALSE ]
  
  ## find the common PTIDs, incase the fits have different ones
  if (!all( rownames(M_x) %in% rownames(M_y))) {
    warning("Not all individuals are shared in common between the two ",
      "fit objects; some will be dropped.")
    common <- intersect( rownames(M_x), rownames(M_y) )
    M_x <- M_x[ rownames(M_x) %in% common, , drop=FALSE]
    M_y <- M_y[ rownames(M_y) %in% common, , drop=FALSE]
    meta_x <- x$data$meta[ c(x$data$individual_id, row_annotation) ]
    meta_y <- y$data$meta[ c(y$data$individual_id, row_annotation) ]
    meta <- meta_x[ meta_x[[x$data$individual_id]] %in% common, ]
  } else {
    meta <- x$data$meta[ c(x$data$individual_id, row_annotation) ]
  }
  
  rowann <- data.frame(.id=rownames(M_x))
  rowann <- merge(
    rowann, 
    meta[c(x$data$individual_id, row_annotation)], 
    by.x=".id",
    by.y=x$data$individual_id
  )
  rowann <- rowann[!duplicated(rowann[[".id"]]), ]
  rownames(rowann) <- rowann[[".id"]]
  rowann <- rowann[-c(which(names(rowann)==".id"))]
  
  ## make sure M, rowann names match up
  rowann <- rowann[ match(rownames(M_x), rownames(rowann)), , drop=FALSE ]
  
  ## get the common categories
  cats <- unique( rbind( 
    x$fit$categories,
    y$fit$categories
  ) )
  
  ## remove the null category
  cats <- cats[ -nrow(cats), ]
  
  cats <- data.frame(cats)
  cats <- cats[,1:(ncol(cats)-1)]
  cats <- as.data.frame( lapply(cats, function(x) {
    factor(x, levels=c(0, 1))
  }))
  
  cats_str <- apply(cats, 1, function(x) {
    paste0(x, collapse="")
  })
  
  ## for all of the categories not in common between M_x, M_y,
  ## set them to zero
  M_x <- as.data.frame(M_x)
  M_y <- as.data.frame(M_y)
  for (cat in cats_str) {
    if (!(cat %in% names(M_x))) {
      M_x[[cat]] <- 0
    }
    if (!(cat %in% names(M_y))) {
      M_y[[cat]] <- 0
    }
  }
  
  ## reorder M_x, M_y
  M_x <- M_x[, order(colnames(M_x)), drop=FALSE]
  M_y <- M_y[, order(colnames(M_y)), drop=FALSE]
  
  if (!all(colnames(M_x) == colnames(M_y))) {
    stop("Internal error: could not match categories between the matrices ",
      "from 'x' and 'y'")
  }
  
  M <- M_x - M_y
  
  ## compute dof from the colnames of M
  dof <- sapply( strsplit( colnames(M), "", fixed=TRUE ), function(x) {
    sum( as.integer(x) )
  })
  
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
  
  print( pheatmap(M[o,],
    color=palette,
    show_rownames=show_rownames,
    show_colnames=show_colnames,
    row_annotation=rowann,
    cluster_rows=FALSE,
    cluster_cols=FALSE,
    cytokine_annotation=cats,
    ...
  ) )
  
  return(M[o, ])
  
}
