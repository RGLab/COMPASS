# This function returns the formatted category name for each row of the categories matrix
# and the ordered categories matrix
# Note that it ORDERS the categories matrix columns first
# Used by mergeMatricesForPlotCOMPASSResultStack function
getCatsAndSubsetNames <- function(cats) {
  # First, re-arrange the columns of the cats data frame to a fixed order (i.e., sorted)
  # This avoids situations where the same category gets called two different names
  # e.g. "M1&M2" vs "M2&M1"
  fixedCatsColsOrder <- c(sort(colnames(cats[, -ncol(cats), drop=FALSE])), "Counts")
  cats <- cats[, fixedCatsColsOrder]
  # Manually create the rownames for the cats data frames, e.g. "M1&!M2&!M3&!M4&!M5&!M6"
  # We compute the subset names manually from the categories matrix because older
  # COMPASS fits don't have them (according to the original COMPASSResult plot function)
  subsets_df <- as.data.frame(cats[, -ncol(cats), drop=FALSE]) # without the Counts column
  for (j in seq_along(subsets_df)) {
    tmp <- subsets_df[[j]]
    subsets_df[[j]] <- paste0(
      swap(tmp, c(0, 1), c("!", "")),
      colnames(subsets_df)[[j]]
    )
  }
  subsets <- do.call( function(...) paste(..., sep="&"), subsets_df )
  list("cats"=cats, "subsets"=subsets)
}

# This function is called by plotCOMPASSResultStack
##' @import data.table
mergeMatricesForPlotCOMPASSResultStack <- function(x,
                                                   threshold=0.01,
                                                   minimum_dof=1,
                                                   maximum_dof=Inf,
                                                   row_annotation,
                                                   variable) {
  # Make sure the marker names for each COMPASSResult are the same.
  # Different orders are fine and are dealt with later.
  if (length(unique(lapply(x, function(c) { sort(colnames(c$data$categories)) }))) != 1) {
    stop("Marker names of all COMPASSResult objects must be the same.")
  }

  # MList, rowannList, and catsList below are each a list of length length(x)
  MList <- lapply(x, function(c) c$fit$mean_gamma)
  # Due to the way rowann objects are constructed, MList objects and rowannList objects will have rows in the same order
  rowannList <- lapply(MList, function(m) data.frame(.id=rownames(m))) # list of data-frames with a .id column, otherwise empty
  catsList <- lapply(x, function(c) c$fit$categories)

  # Prepare the data from each COMPASSResult for merging later on
  for (i in seq_along(x)) {
    #####################
    # The following chunk of code populates rowannList
    #####################

    # the merge() below expects the COMPASS metadata x$data$meta to be a data.frame
    if (is(x[[i]]$data$meta, "data.table")) {
      x[[i]]$data$meta <- as.data.frame(x[[i]]$data$meta)
    }
    # Add a column to the metadata with the value of `variable` for that specific COMPASS run
    # to be pulled into rowannList
    x[[i]]$data$meta[, variable] <- rep(names(x)[[i]], nrow(x[[i]]$data$meta))
    # Populate rowannList (empty so far except for .id in each data frame), merging by id
    rowannList[[i]] <- merge(
      rowannList[[i]],
      # Add the row_annotation columns to the rowann data
      x[[i]]$data$meta[c(x[[i]]$data$individual_id, row_annotation)],
      by.x=".id",
      by.y=x[[i]]$data$individual_id
    )
    # Subset for rows that aren't duplicated
    rowannList[[i]] <- rowannList[[i]][!duplicated(rowannList[[i]][[".id"]]), ,drop=FALSE]
    # Replace rowannList's rownames with the value in the .id column. Then remove the .id column
    rownames(rowannList[[i]]) <- rowannList[[i]][[".id"]]
    rowannList[[i]] <- rowannList[[i]][-c(which(names(rowannList[[i]])==".id"))] # remove the .id column
    # make sure M and rowann rows are in the same order, and their names match up
    rowannList[[i]] <- rowannList[[i]][ match(rownames(MList[[i]]), rownames(rowannList[[i]])), , drop=FALSE ]
    stopifnot(rownames(rowannList[[i]]) == rownames(MList[[i]]))

    #####################
    # The following chunk of code
    # 1) adds rownames to the categories data frame,
    # 2) adds colnames to the mean_gamma data frame, and
    # 3) adds `variable` values to the rownames of the rowann data frame and the mean_gamma data frame
    #####################

    catsData <- getCatsAndSubsetNames(catsList[[i]])
    catsList[[i]] <- catsData$cats # Columns might be re-ordered, nothing else
    subsets <- catsData$subsets
    # Add rownames to catsList[[i]] for the first time
    rownames(catsList[[i]]) <- subsets
    # Assign column names to MList[[i]]
    # NOTE: Making the dangerous assumption here that the row order of the categories data frame
    # is the same as the column order of the M data frame
    # NOTE 2: MList[[i]] may already have colnames, but we re-assign anyway in case the
    # names are something like "M2&M1" instead of "M1&M2"
    colnames(MList[[i]]) <- subsets

    # Add `variable` value, or the name of `x`, to rownames of M and rowann, so we can distinguish
    # between experiments done under different conditions on the same individuals
    rownames(MList[[i]]) <- lapply(rownames(MList[[i]]), function(n) paste(c(n, ", ", names(x)[[i]]), collapse=""))
    rownames(rowannList[[i]]) <- lapply(rownames(rowannList[[i]]), function(n) paste(c(n, ", ", names(x)[[i]]), collapse=""))
  }

  # Obtain all the unique category names across the mean gamma matrices
  allCatNames <- Reduce(union, lapply(MList, colnames))

  # Make all mean gamma matrices have the same categories
  for (i in seq_along(x)) {
    # Get list of columns/categories missing from MList[[i]]
    colsMissing <- setdiff(allCatNames, colnames(MList[[i]]))
    # Add missing categories as columns of MList[[i]], still out of order
    # Values are all 0
    MList[[i]] <- cbind(MList[[i]], sapply(colsMissing,
                                           function(c) c=rep(0, nrow(MList[[i]])))) # c being the name of the category

    # At this point, catsList[[i]] may be a matrix. Turn it into a data frame.
    catsList[[i]] <- as.data.frame(catsList[[i]])
    # Save the name in a column, to be removed later
    catsList[[i]]$name <- rownames(catsList[[i]])
  }

  #####################
  # Finally, we can begin the process of merging all the data frames together
  # in order to obtain one data frame each for catsMerged, MMerged, and rowannMerged
  #####################

  # Bind together the categories matrices
  # Columns should NOT be out of order at this point.

  stopifnot(length(unique(lapply(catsList, function(c) { colnames(c) }))) == 1)
  catsMerged <- setDT(rbindlist(catsList)) %>% setkey(name) %>% unique
  catsMerged <- catsMerged[order(catsMerged$Counts),]
  # Make all columns a factor, except for `name`
  catsMergedRowNamesTmp <- catsMerged$name
  catsMerged <- as.data.frame(sapply(catsMerged, as.factor))
  catsMerged$name <- catsMergedRowNamesTmp # not a factor
  rownames(catsMerged) <- catsMerged$name

  # Put all the mean gamma matrix column names in the same order as the rows of catsMerged
  # Also order the rows in rowannList using row_annotations
  # And then put the MList matrix rows in the same order as those in rowannList
  for (i in seq_along(MList)) {
    MList[[i]] <- as.data.frame(MList[[i]][, catsMerged$name])
    setorderv(rowannList[[i]], row_annotation)
    MList[[i]] <- MList[[i]][ match(rownames(rowannList[[i]]), rownames(MList[[i]])), , drop=FALSE ]
  }

  # Now bind together the mean gamma matrices
  # use column names to make sure numbers are under the correct categories
  MMerged <- as.matrix.noquote(rbindlist(MList, use.names=TRUE))
  # And the rowann matrices
  rowannMerged <- as.data.frame(rbindlist(rowannList))
  rownames(rowannMerged) <- unname(unlist(lapply(rowannList, rownames)))

  #####################
  # Do some filtering and formatting of the merged data frames
  #####################

  ind <- 1:ncol(MMerged)
  dof <- as.numeric(as.matrix.noquote(catsMerged[, "Counts"]))
  # Keep only those categories meeting the min, max dof criteria
  dof_ind <- which(dof >= minimum_dof & dof <= maximum_dof)
  ind <- intersect(ind, dof_ind)
  # Use ind as filter which gets smaller, containing column numbers corresponding to columns indices of M
  # Remove underexpressed categories based on threshold
  # Obtain the mean for each column (category/cell subset) in the mean_gamma matrix
  m_combined <- apply(MMerged, 2, function(x) {
    mean(unlist(x), na.rm = TRUE)
  })
  keep <- m_combined > threshold
  gone <- m_combined <= threshold
  if (length(gone)) {
    message("The 'threshold' filter has removed ", sum(gone),
            " categories:\n", paste( colnames(MMerged)[gone], collapse=", "))
  }
  ind <- intersect(ind, which(keep))
  # Finally, subset
  if (!length(ind)) {
    stop("no marker subsets available for plotting after subsetting")
  }
  MMerged <- MMerged[, ind, drop=FALSE]
  catsMerged <- catsMerged[ind, , drop=FALSE]

  # rbindlist doesn't copy over rownames, so manually assign them to MMerged
  rownames(MMerged) <- unname(unlist(lapply(MList, rownames)))

  # Finally remove Counts column from catsMerged after assigning rownames
  rownames(catsMerged) <- catsMerged$name
  catsMerged <- as.data.frame(catsMerged[,setdiff(names(catsMerged), c("name", "Counts"))])

  # MMerged may be a noquote matrix object. Convert to a matrix object that has columns of numeric type (not a list)
  MMerged <- apply(as.matrix(unclass(MMerged)), 2, unlist)

  list("MMerged"=MMerged, "rowannMerged"=rowannMerged, "catsMerged"=catsMerged)
}

##' Plot multiple COMPASSResults
##'
##' This function can be used to visualize the mean probability of response;
##' that is, the probability that there is a difference in response between
##' samples subjected to the 'treatment' condition, and samples subjected
##' to the 'control' condition.
##' This version is used for plotting multiple COMPASSResult objects.
##' The COMPASS runs must all use the same markers.
##' The code is heavily based on the plot.COMPASSResult and plot2 functions.
##' Not all options from plot.COMPASSResult are supported yet.
##'
##' @param x A named list of objects of class \code{COMPASSResult}. The names are values of type \code{variable}
##' @param threshold A numeric threshold for filtering under-expressed
##'   categories. Any categories with mean score < \code{threshold} are
##'   removed.
##' @param minimum_dof The minimum degree of functionality for the categories
##'   to be plotted.
##' @param maximum_dof The maximum degree of functionality for the categories
##'   to be plotted.
##' @param row_annotation A vector of names, pulled from the metadata, to be
##'   used for row annotation.
##' @param variable What to call the variable that is different across the objects in x.
##' @param palette The colour palette to be used.
##' @param show_rownames Boolean; if \code{TRUE} we display row names (ie,
##'   the individual ids).
##' @param ... Optional arguments passed to \code{pheatmap}.
##' @importFrom RColorBrewer brewer.pal
##' @importFrom grDevices colorRampPalette
##' @import data.table
##' @import magrittr
##' @return The plot as a \code{grid} object (\code{grob}). It can be redrawn
##' with e.g. \code{grid::grid.draw()}.
##' @export
##' @examples
##' \dontrun{
##' # allCompassResults is a list of 4 COMPASSResults
##' names(allCompassResults) <- c("Antigen 85A", "CFP-10", "CMV", "ESAT-6")
##' plotCOMPASSResultStack(allCompassResults,
##'     row_annotation = c("Antigen", "PATIENT ID", "Time"),
##'     variable = "Antigen", show_rownames = FALSE,
##'     main = "Heatmap of Mean Probability of Response to Antigens, CD8+",
##'     fontsize = 14, fontsize_row = 13, fontsize_col = 11)
##'   }
plotCOMPASSResultStack <- function(x,
                                   threshold=0.01,
                                   minimum_dof=1,
                                   maximum_dof=Inf,
                                   row_annotation,
                                   variable,
                                   palette=colorRampPalette(brewer.pal(9,"Purples"))(20),
                                   show_rownames=FALSE,
                                   ...) {

  data2plot <- mergeMatricesForPlotCOMPASSResultStack(x,
                                                      threshold,
                                                      minimum_dof,
                                                      maximum_dof,
                                                      row_annotation,
                                                      variable)

  # Now we have MMerged, rowannMerged, and catsMerged. Time to plot:
  pheatmap(data2plot$MMerged,
           color=palette,
           show_rownames=show_rownames,
           show_colnames=show_colnames,
           row_annotation=data2plot$rowannMerged,
           cluster_rows=FALSE,
           cluster_cols=FALSE,
           cytokine_annotation=data2plot$catsMerged,
           ...)
}
